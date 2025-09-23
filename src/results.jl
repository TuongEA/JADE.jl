#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

"""
    write_sim_results(results, data::JADEData, parameters::JADESimulation)

    This function writes results from a simulation.

    ### Required arguments
    `results` is a dictionary of results.

    `data` is a JADE data object for the problem simulated.

    `parameters` contains all the simulation settings.
"""
function write_sim_results(results::Array{Array{Dict{Symbol,Any},1},1}, d::JADEData, parameters::JADESimulation)

    # Just aliases
    nsims = parameters.replications
    number_of_wks = d.rundata.number_of_wks * parameters.number_of_cycles

    if d.rundata.discount == 0.0 || parameters.reset_starting_levels == true
        number_of_wks -= parameters.initial_stage - 1
    end

    T = d.durations
    s = d.sets

    data_dir = joinpath(@__JADE_DIR__,"Output",d.rundata.data_dir,d.rundata.policy_dir,parameters.sim_dir)

    # Prepare output directory
    !ispath(data_dir) && mkpath(data_dir)

    # How we will store information from a variable (except inflows)
    function adddata(getvalue::Function, header::String, filename::String)
        open(filename, "w") do f
            # Headers
            write(f, "simulation, stage, ", header)
            # Data values
            for i in 1:nsims
                for j in 1:number_of_wks
                    # Get the value from the simulation
                    value = getvalue(i, j)
                    write(f, "\n", string(i), ",", string(j), ",", string(value))
                end
            end
        end
    end

    # How we will store information if load blocks are the headers
    function adddata(getvalue::Function, filename::String)
        open(filename, "w") do f
            # Headers
            write(f, "simulation, stage")
            for bl in s.BLOCKS
                write(f, ",", string(bl))
            end
            # Data values
            for i in 1:nsims
                for j in 1:number_of_wks
                    write(f, "\n", string(i), ",", string(j))
                    for bl in s.BLOCKS
                        value = getvalue(i, j, bl)
                        write(f, ",", string(value))
                    end
                end
            end
        end
    end


    #------------------------------------
    # Final objective
    #------------------------------------
    open(joinpath(data_dir,"TotalCost.csv"), "w") do io
        for i in 1:nsims
            println(io, results[i][end][:running_cost])
        end
        return
    end

    #------------------------------------
    # Stored energy
    #------------------------------------
    adddata("energy_GWh", joinpath(data_dir,"StoredEnergy.csv")) do i, j
        return sum(d.reservoirs[r].sp * results[i][j][:reslevel][r].out * 10^3 for r in s.RESERVOIRS) # Define function for getvalue(i,j)
    end

    #------------------------------------
    # Arc flows
    #------------------------------------
    !ispath(joinpath(data_dir,"NaturalFlows")) && mkpath(joinpath(data_dir,"NaturalFlows"))

    for a in s.NATURAL_ARCS # Write results for each water arc
        adddata(joinpath(data_dir,string("NaturalFlows/", string(a[1]), "_", string(a[2]), ".csv"))) do i, j, bl
            return results[i][j][:naturalflows][a, bl]
        end
    end

    #------------------------------------
    # Spill flows
    #------------------------------------
    !ispath(joinpath(data_dir,"SpillFlows")) && mkpath(joinpath(data_dir,"SpillFlows"))

    for a in s.STATION_ARCS # Write results for each station's hydro arc
        adddata(joinpath(data_dir,string("SpillFlows/", string(a[1]), "_", string(a[2]), ".csv"))) do i, j, bl
            return results[i][j][:spills][a, bl]
        end
    end

    #------------------------------------
    # Release flows
    #------------------------------------
    !ispath(joinpath(data_dir,"ReleaseFlows")) && mkpath(joinpath(data_dir,"ReleaseFlows"))

    for a in s.STATION_ARCS # Write results for each station's hydro arc
        adddata(joinpath(data_dir,string("ReleaseFlows/", string(a[1]), "_", string(a[2]), ".csv"))) do i, j, bl
            return results[i][j][:releases][a, bl]
        end
    end

    #------------------------------------
    # Output for lower bound penalties
    #------------------------------------
    # Calculate amount of flow below bound, in cumecs * hours, scenario i, stage j
    function flowunder(i, j)
        timenow = TimePoint(d.rundata.start_yr, d.rundata.start_wk) + (j - 1) % d.rundata.number_of_wks
        return sum([
            results[i][j][:flowunder][a, bl] * T[timenow][bl] for
            a in s.NATURAL_ARCS, bl in s.BLOCKS if d.natural_arcs[a].minflow != 0.0
        ])
    end

    adddata("lb_penalty", joinpath(data_dir,"FlowLBCost.csv")) do i, j
        return flowunder(i, j) * SECONDSPERHOUR * d.spMax * d.rundata.penalty_lb
    end

    #------------------------------------
    # Output for upper bound penalties
    #------------------------------------
    # Calculate amount of flow above bound, in cumecs * hours, scenario i, stage j
    function flowover(i, j)
        timenow = TimePoint(d.rundata.start_yr, d.rundata.start_wk) + (j - 1) % d.rundata.number_of_wks
        return sum([
            results[i][j][:flowover][a, bl] * T[timenow][bl] for
            a in s.NATURAL_ARCS, bl in s.BLOCKS if d.natural_arcs[a].maxflow != Inf
        ])
    end

    adddata("ub_penalty", joinpath(data_dir,"FlowUBCost.csv")) do i, j
        return flowover(i, j) * SECONDSPERHOUR * d.spMax * d.rundata.penalty_ub
    end

    #------------------------------------
    # Output for thermal costs
    #------------------------------------
    adddata("thermal_cost", joinpath(data_dir,"ThermalCosts.csv")) do i, j
        timenow = TimePoint(d.rundata.start_yr, d.rundata.start_wk) + (j - 1) % d.rundata.number_of_wks
        return sum([
            results[i][j][:thermal_use][name, bl] * station.heatrate * d.fuel_costs[timenow][station.fuel] * T[timenow][bl] for
            (name, station) in d.thermal_stations, bl in s.BLOCKS
        ])
    end

    #------------------------------------
    # Lost Load
    #------------------------------------
    adddata("cost", joinpath(data_dir,"LostLoadCost.csv")) do i, j
        return results[i][j][:lostloadcosts]
    end

    #------------------------------------
    # Total present cost
    #------------------------------------
    adddata("cost", joinpath(data_dir,"PresentCost.csv")) do i, j
        return results[i][j][:stage_objective]
    end

    #------------------------------------
    # Future costs
    #------------------------------------
    adddata("cost", joinpath(data_dir,"FutureCost.csv")) do i, j
        return results[i][j][:bellman_term]
    end

    adddata("cost", joinpath(data_dir,"SummedCosts.csv")) do i, j
        return results[i][j][:stage_objective] + results[i][j][:bellman_term]
    end

    #------------------------------------
    # Spilled energy
    #------------------------------------
    spillarc(station) = d.hydro_stations[station].arc
    function spillsum(i, j, bl)
        timenow = TimePoint(d.rundata.start_yr, d.rundata.start_wk) + (j - 1) % d.rundata.number_of_wks
        return sum(results[i][j][:spills][spillarc(m), bl] * d.hydro_stations[m].sp * T[timenow][bl] for m in s.HYDROS)
    end

    for bl in s.BLOCKS
        adddata("energy_MWh", joinpath(data_dir,"SpilledEnergy_$(bl).csv")) do i, j
            return spillsum(i, j, bl)
        end
    end

    #------------------------------------
    # Contingent storage cost
    #------------------------------------
    adddata("contingent_storage_cost_\$", joinpath(data_dir,"ContingentStorageCost.csv")) do i, j
        return results[i][j][:contingent_storage_cost]
    end

    #------------------------------------
    # Archive run file
    #------------------------------------
    open(joinpath(joinpath(data_dir,"sim_parameters.json")), "w") do f
        return JSON.print(f, parameters)
    end

    return nothing
end

"""
This function saves output from policy generation.
"""
function write_training_results(sddpm::SDDP.PolicyGraph, d::JADEData, solveoptions::JADESolveOptions)
    data_dir = joinpath(@__JADE_DIR__, "Output", d.rundata.data_dir)

    !ispath(data_dir) && mkpath(data_dir)

    #------------------------------------
    # Archive run and solve related data
    #------------------------------------
    !ispath(joinpath(data_dir, d.rundata.policy_dir)) && mkpath(joinpath(data_dir, d.rundata.policy_dir))

    @info("Archiving training parameters...")
    open(joinpath(data_dir, d.rundata.policy_dir, "rundata.json"), "w") do f
        return JSON.print(f, d.rundata)
    end

    open(joinpath(data_dir, d.rundata.policy_dir, "solveoptions.json"), "w") do f
        return JSON.print(f, solveoptions)
    end

    @info("Writing log to file...")
    open(joinpath(data_dir, d.rundata.policy_dir, "convergence.json"), "w") do f
        return JSON.print(f, sddpm.most_recent_training_results.log)
    end

    @info("Writing cuts to JSON file...")
    SDDP.write_cuts_to_file(sddpm, joinpath(data_dir, d.rundata.policy_dir, "cuts.json"))

    if solveoptions.write_eohcuts
        if d.rundata.number_of_wks == 52 && d.rundata.steady_state == true
            EOH_dir = joinpath(@__JADE_DIR__, "Output", d.rundata.data_dir, "EOH")
            @info("Creating EOH cuts in " * EOH_dir )

            !ispath(EOH_dir) && mkpath(EOH_dir)

            JADE.write_cuts_to_file(sddpm, joinpath(EOH_dir, "$(d.rundata.policy_dir).eoh"), d.rundata.start_wk)

            open(joinpath(EOH_dir, "$(d.rundata.policy_dir).rdt"), "w") do f
                return JSON.print(f, d.rundata)
            end
            open(joinpath(EOH_dir, "$(d.rundata.policy_dir).slv"), "w") do f
                return JSON.print(f, solveoptions)
            end
        else
            @warn(
                "End-of-horizon cuts not written, since these must be generated\nfrom a 52-week infinite horizon model."
            )
        end
    end

    @info("Writing cuts to DOASA-compatible files...")
    write_DOASA_cuts(sddpm, d, joinpath(data_dir, d.rundata.policy_dir, "Cuts"))

    return nothing
end

function output_tidy_results(results::Vector{Vector{Dict{Symbol,Any}}}, d::JADEData, parameters::JADESimulation; variables::Array{Symbol,1})
    data_dir = joinpath(@__JADE_DIR__, "Output", d.rundata.data_dir, d.rundata.policy_dir, parameters.sim_dir)

    # Prepare output directory
    !ispath(data_dir) && mkpath(data_dir)

    open(joinpath(data_dir,"tidy_results.csv"), "w") do f
        println(f, "simulation,stage,year,week,type,name,demand,node,value")
        for i in 1:length(results)
            time = TimePoint(d.rundata.start_yr, d.rundata.start_wk) + (parameters.initial_stage - 1)
            for j in 1:length(results[i])
                for s in variables
                    if s ∉ keys(results[i][j])
                        error(string(s) * " not found")
                    end
                    data = results[i][j][s]
     
                    temp = "$i,$((j + parameters.initial_stage - 2) % WEEKSPERYEAR + 1),$(time.year),$(time.week),$s"

                    temp2 = ""
                    if typeof(data) <: AbstractArray || typeof(data) <: Dict
                        if typeof(data) <: JuMP.Containers.SparseAxisArray
                            data = data.data
                        end
                        for k in keys(data)
                            if typeof(k) <: JuMP.Containers.DenseAxisArrayKey
                                if length(getfield(k, 1)) == 1
                                    temp2 = replace(string(getfield(k, 1)[1]), "," => ";") * ",-,-"

                                elseif length(getfield(k, 1)) == 2
                                    if s == :prices
                                        temp2 =  "-," * replace(string(getfield(k, 1)[2]), "," => ";") * "," * string(getfield(k, 1)[1])
                                    else
                                        temp2 = replace(string(getfield(k, 1)[1]), "," => ";") * "," * string(getfield(k, 1)[2]) * ",-"
                                    end
                                
                                elseif length(getfield(k, 1)) == 3
                                    temp2 = replace(string(getfield(k, 1)[3]), "," => ";") * "," * string(getfield(k, 1)[2]) * "," * string(getfield(k, 1)[1])
                                end

                                if typeof(data[k]) <: SDDP.State
                                    if j == 1
                                        time2 = time - 1
                                        println(f, join([i, (j + parameters.initial_stage - 2) % WEEKSPERYEAR, time2.year, time2.week, s, temp2, data[k].in], ","))
                                    end

                                    println(f,join([temp, temp2 ,string(data[k].out)], ","))
                                else
                                    println(f, temp * "," * temp2 * "," * string(data[k]))
                                end

                            elseif typeof(data[k]) <: SDDP.State
                                if j == 1
                                    time2 = time - 1
                                    println(f, join([i, (j + parameters.initial_stage - 2) % WEEKSPERYEAR, time2.year, time2.week, s, k, "-,-",data[k].in], ","))
                                end
                                println(f, join([temp,k,"-,-",data[k].out], "," ))

                            elseif typeof(k) <: Tuple
                                if length(k) == 1
                                    println(f, join([i, (j + parameters.initial_stage - 2) % WEEKSPERYEAR + 1, time.year, time.week, s, k[1], "-,-", data[k]], ","))

                                elseif length(k) == 2
                                    println(f, join([i, (j + parameters.initial_stage - 2) % WEEKSPERYEAR + 1, time.year, time.week, s, k[1], k[2], "-",data[k]], ","))

                                elseif length(k) == 3
                                    println(f, join([i, (j + parameters.initial_stage - 2) % WEEKSPERYEAR + 1, time.year, time.week, s, 
                                                     replace(string(k[3]), "," => ";"), k[2], k[1], data[k]], ","))
                                end
                            else
                                println(f, join([temp,k,"-,-", data[k]], ","))
                            end
                        end
                    else
                        println(f, join([temp, "-,-,-", data], ","))
                    end
                end
                time += 1
            end
        end
    end
    return nothing
end

function write_DOASA_cuts(sddpm::SDDP.PolicyGraph, d::JADEData, path::String)
    !ispath(path) && mkpath(path)
    
    for t in keys(sddpm.nodes)
        oracle = sddpm.nodes[t].bellman_function.global_theta
        states = sddpm.nodes[t].bellman_function.global_theta.states
        states = collect(keys(states))

        states_ordered = Symbol[]
        for i in 1:length(states)
            push!(states_ordered, :a)
        end

        for s in states
            
            if findfirst("reslevel[", string(s)) != nothing
                s_sym = Symbol(string(s)[10:end-1])
                states_ordered[d.reservoirs[s_sym].index] = s

            elseif findfirst("fuelstoragelevel[", string(s)) != nothing
                s_sym = string(s)[20:end-2]                 # Get substring from fuelstoragelevel[(:GENE, :GAS)] => GENE, :GAS
                s_sym = replace(s_sym," :" =>"")            # Repalce " :" by "" --> GENE, :GAS => GEN,GAS
                parts = split(s_sym, ",")                   # Split GEN,GAS => parts[1] = GEN and parts[2] = GAS
                s_sym = Symbol(parts[1]), Symbol(parts[2])  # Define Tuple{Symbol,Symbol} => (:GENE, :GAS)
               
                states_ordered[length(d.sets.RESERVOIRS) + findfirst(isequal(s_sym), d.sets.STORED_FUELS)] = s
            
            else # Not sure the use of this
                s_sym = Symbol(string(s)[14:end-1])
                states_ordered[length(d.sets.RESERVOIRS,) + findfirst(isequal(s_sym), d.sets.STORED_FUELS)] = s
            end
        end

        open(joinpath(path, "BendersCuts_" * string(t) * "_1.csv"), "w") do f
            for k in 1:length(oracle.cuts)
                print(f, oracle.cuts[k].intercept)
                print(f, ",")

                for s in states_ordered
                    print(f, "," * string(-oracle.cuts[k].coefficients[s]))
                end
                println(f, "")
            end
        end
    end
    return nothing
end
