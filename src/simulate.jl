#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

using Random

"""
	simulate(JADEmodel::JADEModel, parameters::JADESimulation)

This function carries out a simulation of a JADE model; all the specifications
for this simulation are defined within the `parameters` argument.

The output from this simulation will automatically be written into the
Output/`data_dir`/`policy_dir`/`sim_dir` subdirectory.

### Required Arguments
`JADEmodel` is the model detailing the scenario that we wish to simulate, and
all the corresponding data.

`parameters` contains all the simulation information, including the number of
replications, the type of simulation, the hydrological years to sample from, etc.
"""
function simulate(JADEmodel::JADEModel, parameters::JADESimulation)
    d = JADEmodel.d
    sddpm = JADEmodel.sddpm

    check_settings_compatibility(rundata = d.rundata, simulation = parameters)

    @info("Simulating policy '" * d.rundata.policy_dir * "' for model data in " * joinpath("Input", d.rundata.data_dir))

    cuts_path = joinpath(@JADE_DIR, "Output", d.rundata.data_dir, d.rundata.policy_dir, "cuts.json")

    if length(sddpm.nodes[1].bellman_function.global_theta.cuts) == 0 # If no cuts are loaded in memory, try loading them from disk:
        if isfile(cuts_path)
            @info("Loading existing cuts file: " * joinpath("Output", d.rundata.data_dir, d.rundata.policy_dir, "cuts.json"))
            previous_rundata = load_model_parameters(d.rundata.data_dir, d.rundata.policy_dir)
            check_rundata(d.rundata, previous_rundata, :partial)
            SDDP.read_cuts_from_file(sddpm, cuts_path)
            if has_upper_bound(sddpm.nodes[d.rundata.number_of_wks].bellman_function.global_theta.theta) # Removes any upper bound on the final stage if it exists.
                delete_upper_bound(sddpm.nodes[d.rundata.number_of_wks].bellman_function.global_theta.theta)
            end
        else
            @info("No cuts file found in " * joinpath("Output", d.rundata.data_dir, d.rundata.policy_dir))
        end
    end
    
    # Defines a list of primal variables to extract from the simulation. These are decision variables like generation, reservoir levels, spills, etc.
    get_primal = [
        :thermal_use,
        :hydro_disp,
        :naturalflows,
        :releases,
        :spills,
        :reslevel,
        :inflow,
        :flowover,
        :flowunder,
        :flowpenalties,
        :lostload,
        :terminalcost,
        :transflow,
        :lostloadcosts,
        :contingent_storage_cost,
        :carbon_emissions,
    ]
    if ApplyFuelConstraints
        fuel_primal = [:fuelstoragelevel, :fuel_use_TJ, :fuel_contract, :fuel_injection, :fuel_withdrawal ]
        get_primal = [get_primal; fuel_primal]
    end

    # Defines dual variables to extract, such as:
    get_dual = Dict{Symbol,Function}(
        :prices => (sp) -> d.rundata.scale_objective * JuMP.dual.(sp[:defineShedding]),  # shadow prices of demand.
        :mwv    => (sp) -> -d.rundata.scale_objective * JuMP.dual.(sp[:rbalance]) / 1E3 / d.rundata.scale_reservoirs, # marginal water value, scaled appropriately.
    )
    if ApplyFuelConstraints
        get_dual[:mfsv] = (sp) -> -d.rundata.scale_objective * JuMP.dual.(sp[:fuelStorageBalance]) / 1E3  # marginal fuel value
    end

    # Initial State Handling
    initial_state = SDDP._initial_state(JADEmodel.sddpm) # Gets the default initial state from the model
    if parameters.initial_state == nothing # Then checks if a custom initial state is provided:
        @info("Using default initial reservoir levels")
    else # Validate and use custom initial state
        if length(initial_state) != length(parameters.initial_state)
            @info("An initial value must be specified for each of the following states:")
            for key in keys(initial_state)
                @info(key[10:end-1])
            end
            error("initial_state dictionary has incorrect number of states")
        end
        for key in keys(parameters.initial_state)
            if key ∉ keys(initial_state)
                @info("An initial value must be specified for each of the following states:")
                for key in keys(initial_state)
                    @info(key[10:end-1])
                end
                error("Invalid state: " * key)
            end
        end
        initial_state = parameters.initial_state
        @info("Using simulation-specific initial reservoir levels")
    end

    # Simulation Type Logging
    if parameters.sim_type == :monte_carlo
        @info("Simulating Monte Carlo inflows with " * string(parameters.replications) * " replications")
    elseif parameters.sim_type == :historical
        if parameters.randomize_years == false
            @info("Simulating historical inflows in sequence provided")
        else
            @info("Simulating annual historical inflows randomly with " * string(parameters.replications) * " replications")
        end
    end

    Random.seed!(parameters.random_seed) # Sets the random seed for reproducibility.
    wks = d.rundata.number_of_wks * parameters.number_of_cycles      # Calculates the total number of weeks to simulate.
    lastwk = d.rundata.number_of_wks * parameters.number_of_cycles   # Define last week of simulation
    # Adjusts based on whether the model is in steady state:
    if !d.rundata.steady_state
        wks += 1 - parameters.initial_stage
    else
        lastwk += parameters.initial_stage - 1
    end

    primaldualvariables = vcat(get_primal, [:stage_objective, :bellman_term, :prices, :mwv])
    if ApplyFuelConstraints
        primaldualvariables = [primaldualvariables; :mfsv]
    end

    if parameters.sim_type == :monte_carlo # Monte Carlo Simulation
        if !d.rundata.steady_state || parameters.reset_starting_levels == true  # Case 1: Not Steady State or Resetting Starting Levels
            results = SDDP.simulate(
                sddpm,
                parameters.replications,
                get_primal,
                custom_recorders = get_dual,
                sampling_scheme = InSampleMonteCarlo2(
                    max_depth = wks,
                    initial_stage = parameters.initial_stage,
                    terminate_on_cycle = false,
                    terminate_on_dummy_leaf = false,
                ),
                incoming_state = initial_state,
            )
            # Post-processing Results
            for i in 1:parameters.replications
                for τ in parameters.initial_stage:lastwk
                    t = τ - parameters.initial_stage + 1
                    for key in keys(results[i][t][:reslevel]) # Rescaling Reservoir Levels
                        results[i][t][:reslevel][key[1]] = SDDP.State(
                            results[i][t][:reslevel][key[1]].in  * d.rundata.scale_reservoirs,
                            results[i][t][:reslevel][key[1]].out * d.rundata.scale_reservoirs,
                        )
                    end
                    # Assigning Inflow Year
                    if results[i][t][:noise_term][:scenario] == 0
                        results[i][t][:inflow_year] = d.rundata.start_yr
                    else
                        results[i][t][:inflow_year] = d.rundata.sample_years[round(Int, results[i][t][:noise_term][:scenario])]
                    end
                end
            end
        else # Case 2: Steady State and No Reset
            sequence = SDDP.simulate(
                sddpm,
                1,                                               # Instead of running multiple replications separately, this runs one long simulation 
                get_primal,
                custom_recorders = get_dual,
                sampling_scheme = InSampleMonteCarlo2(
                    max_depth = wks * parameters.replications,   # with wks * replications stages. This is more efficient when the model is in steady state.
                    initial_stage = parameters.initial_stage,
                    terminate_on_cycle = false,
                    terminate_on_dummy_leaf = false,
                ),
                incoming_state = initial_state,
            )
            # Splitting Results into Replications
            results = Vector{Dict{Symbol,Any}}[]    # Creates a nested structure to hold results per replication and per stage.
            for i in 1:parameters.replications
                push!(results, Dict{Symbol,Any}[])
                for τ in parameters.initial_stage:lastwk
                    t = τ - parameters.initial_stage + 1
                    temp = Dict{Symbol,Any}()
                    push!(results[i], temp)

                    for sym in primaldualvariables # Extracting and Rescaling Variables
                        if sym == :reslevel
                            results[i][t][sym] = Dict{Symbol,SDDP.State}()
                            for key in keys(sequence[1][(i-1)*wks+t][sym])
                                results[i][t][sym][key[1]] = SDDP.State(
                                    sequence[1][(i-1)*wks+t][sym][key].in  * d.rundata.scale_reservoirs,
                                    sequence[1][(i-1)*wks+t][sym][key].out * d.rundata.scale_reservoirs,
                                )
                            end
                        else
                            results[i][t][sym] = sequence[1][(i-1)*wks+t][sym]
                        end
                    end
                    #Assigning Inflow Year
                    if sequence[1][(i-1)*wks+t][:noise_term][:scenario] == 0
                        results[i][t][:inflow_year] = d.rundata.start_yr
                    else
                        results[i][t][:inflow_year] = d.rundata.sample_years[round(Int,sequence[1][(i-1)*wks+t][:noise_term][:scenario])]
                    end
                end
            end
        end
    elseif parameters.sim_type == :historical # Historical Simulation
        sequence = nothing                      # Initialises the output container `sequence` to hold the raw simulation output.
        results = Vector{Dict{Symbol,Any}}[]    # Initialises the output container `results` to store processed results per replication and stage.

        if !d.rundata.steady_state || parameters.reset_starting_levels == true # Case 1: Not Steady State or Resetting Starting Levels
            sample_paths = Vector{Tuple{Int,Dict{Symbol,Float64}}}[]   # Creates a list of inflow paths, each being a sequence of (week, inflow_data) tuples.
            push!(sample_paths, Tuple{Int,Dict{Symbol,Float64}}[])
            seq = []
            count = 0
            for year in parameters.sim_years # Loop Over Simulation Years
                years = collect(year:(year-1+ceil(Int, (d.rundata.start_wk - 1 + lastwk) / WEEKSPERYEAR)))  # Compute the list of years needed to cover the simulation horizon.
                inflow_mat = getinflows_for_historical(get_file_directory("inflows.csv", d.rundata),d.rundata,years) # Load inflow data using getinflows_for_historical.
                i = 1
                extrawks = d.rundata.steady_state ? parameters.initial_stage - 1 : 0
                for t in parameters.initial_stage:(d.rundata.number_of_wks+extrawks)  # For each week:
                    s_inflows = Dict{Symbol,Float64}()                                # Create a dictionary s_inflows 
                    for c in d.sets.CATCHMENTS_WITH_INFLOW
                        s_inflows[c] = inflow_mat[(t+d.rundata.start_wk-2)%WEEKSPERYEAR+1][c][i]  # with inflows for each catchment.
                    end
                    # Assign a :scenario tag to indicate the inflow year.
                    if d.rundata.first_week_known && t == 1
                        s_inflows[:scenario] = d.rundata.start_yr 
                    else
                        s_inflows[:scenario] = years[i]
                    end

                    if (t + d.rundata.start_wk - 2) % WEEKSPERYEAR == WEEKSPERYEAR - 1
                        i += 1
                    end
                    push!(sample_paths[end], ((t - 1) % WEEKSPERYEAR + 1, s_inflows)) # Push the inflow data into the current sample path.
                end
                count += 1
                if count == parameters.number_of_cycles
                    count = 0
                    push!(sample_paths, Tuple{Int,Dict{Symbol,Float64}}[])
                end
            end
            
            # Runs the simulation using the historical inflow paths.
            sims = SDDP.simulate(
                sddpm,
                parameters.replications,
                get_primal,
                sampling_scheme = SDDP.Historical(sample_paths),
                custom_recorders = get_dual,
                incoming_state = initial_state,
            )

            # Flattens the simulation results into a single sequence.
            for sim in sims
                append!(seq, sim)
            end
            sequence = [seq]

        else # Case 2: Steady State and No Reset
            sample_path = Tuple{Int,Dict{Symbol,Float64}}[] # Weekly Inflow Construction
            for year in parameters.sim_years                # Similar to Case 1, but without cycle splitting. All inflows are pushed into one path.
                years = collect(year:(year-1+ceil(Int, (d.rundata.start_wk + d.rundata.number_of_wks - 1 + parameters.initial_stage - 1) / WEEKSPERYEAR)))
                inflow_mat = getinflows_for_historical(get_file_directory("inflows.csv", d.rundata),d.rundata,years)
                i = 1

                for t in
                    parameters.initial_stage:(d.rundata.number_of_wks+parameters.initial_stage-1)
                    s_inflows = Dict{Symbol,Float64}()
                    for c in d.sets.CATCHMENTS_WITH_INFLOW
                        s_inflows[c] =
                            inflow_mat[(t+d.rundata.start_wk-2)%WEEKSPERYEAR+1][c][i]
                    end
                    s_inflows[:scenario] = years[i]
                    if (t + d.rundata.start_wk - 2) % WEEKSPERYEAR == WEEKSPERYEAR - 1
                        i += 1
                    end
                    push!(sample_paths[end], ((t - 1) % WEEKSPERYEAR + 1, s_inflows))
                end
            end
            # Runs the simulation using the single historical path.
            sequence = SDDP.simulate(
                sddpm,
                1,
                get_primal,
                sampling_scheme = SDDP.Historical(sample_path),
                custom_recorders = get_dual,
                incoming_state = initial_state,
            )
        end
        # Post-Processing Results
        for i in 1:parameters.replications
            push!(results, Dict{Symbol,Any}[])  # Initialises result containers for each replication.
            for τ in parameters.initial_stage:lastwk # Loop Over Stages. For each stage:
                t = τ - parameters.initial_stage + 1
                temp = Dict{Symbol,Any}()
                push!(results[i], temp)
                for sym in primaldualvariables  # Extract and rescale variables.
                    if sym == :reslevel
                        results[i][t][sym] = Dict{Symbol,SDDP.State}()
                        for key in keys(sequence[1][(i-1)*wks+t][sym])
                            results[i][t][sym][key[1]] = SDDP.State(
                                sequence[1][(i-1)*wks+t][sym][key].in  * d.rundata.scale_reservoirs,
                                sequence[1][(i-1)*wks+t][sym][key].out * d.rundata.scale_reservoirs,
                            )
                        end
                    else
                        results[i][t][sym] = sequence[1][(i-1)*wks+t][sym]
                    end
                end
                results[i][t][:inflow_year] = round(Int, sequence[1][(i-1)*wks+t][:noise_term][:scenario]) # Assign inflow year from :noise_term.
            end
        end
    end

    # Final Post-Processing Loop
    for i in 1:parameters.replications
        for τ in parameters.initial_stage:lastwk
            t = τ - parameters.initial_stage + 1
            # Scaling Stage Results
            results[i][t][:stage_objective] *= d.rundata.scale_objective
            results[i][t][:bellman_term] *= d.rundata.scale_objective
            for r in d.sets.RESERVOIRS
                results[i][t][:mwv][r] /= d.reservoirs[r].sp # Divides the marginal water value for each reservoir by its specific power factor (sp) to normalize it.
            end
            # Computes the cumulative running cost over time by summing up the stage objectives.
            if t == 1
                results[i][t][:running_cost] = results[i][t][:stage_objective]
            else
                results[i][t][:running_cost] =
                    results[i][t-1][:running_cost] + results[i][t][:stage_objective]
            end

            # Computes the total system storage
            results[i][t][:total_storage] = 0
            # by summing the out level of each reservoir,
            for r in keys(d.reservoirs)
                results[i][t][:total_storage] += results[i][t][:reslevel][r].out * d.reservoirs[r].sp * 1000 # Convert Million m3 to GWH
            end
        end
    end

    @info("Saving output in " * joinpath("Output", d.rundata.data_dir, d.rundata.policy_dir, parameters.sim_dir))
    # write_sim_results(results, d, parameters)

    tidy_variables = [:reslevel, :thermal_use, :transflow, :prices, :lostload, :flowover, :flowunder, :flowpenalties, 
                      :contingent_storage_cost, :carbon_emissions, :spills, :total_storage, :inflow_year, :mwv]
    if ApplyFuelConstraints
        tidy_variables = [tidy_variables;[:fuelstoragelevel, :fuel_contract, :fuel_use_TJ,:fuel_injection, :fuel_withdrawal, :mfsv]]
    end 
    output_tidy_results(
        results,
        d,
        parameters,
        variables = tidy_variables,
    )

    @info("Done.")
    return results
end
