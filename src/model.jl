#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

"""
    JADEsddp(d::JADEData, optimizer = nothing)

    This function builds a JADE model using SDDP.jl from the JADEData object `d`.
    Returns an `SDDPModel` object.

    ### Required Arguments
    `d` is a `JADEData` object.
    `optimizer` is a `JuMP` optimizer.
"""
function JADEsddp(d::JADEData, optimizer = nothing)

    # Constants
    penalty_ub = d.rundata.penalty_ub
    penalty_lb = d.rundata.penalty_lb
    number_of_wks = d.rundata.number_of_wks
    nscenarios = d.rundata.nscenarios
    nmargins = length(d.terminal_eqns)
    scale_factor = d.rundata.scale_reservoirs
    scale_obj = d.rundata.scale_objective

    @assert nmargins > 0                             # Ensures that terminal conditions are defined.

    # Some more convenient names
    LOOPS = d.loops                                  # Transmission network loops (used for Kirchhoff’s laws or loss modeling) - seems to be redundant
    s = d.sets                                       # model sets

    if optimizer == nothing                          # If no optimizer is provided
        error("No solver specified")                 # throws an error.
    elseif typeof(optimizer) <: Function             # If a function is passed (e.g., ()->Gurobi.Optimizer)
        d.parallel_optimizer = optimizer             # it is stored 
        optimizer = d.parallel_optimizer()           # and then called to create an optimizer instance.
    end

    #------------------------------------------------------------------------
    graph = SDDP.LinearGraph(number_of_wks)          # Creates a linear-stage graph with number_of_wks stages (e.g., 52 for a year)

    if d.rundata.steady_state                        # If the model is steady-state
        SDDP.add_edge(graph, number_of_wks => 1, d.rundata.discount)  # Connects the last stage back to the first, forming a loop.
    end

    if d.rundata.weekly_discounting && d.rundata.discount != 0     # Applies weekly discounting to future costs if enabled. 
        wk_disc = d.rundata.discount^(1 / 52)                      # Adjusts the discount factor to reflect weekly compounding.
        for (index, obj) in graph.nodes
            graph.nodes[index] = [(obj[1][1], wk_disc)]            # Modifying the discount factors associated with each node (stage) in the SDDP graph
        end
    end

    sddpm = SDDP.PolicyGraph(
        graph,                                        # Initializes the SDDP model using the "graph".
        sense = :Min,                                 # Sets the objective to minimization.
        lower_bound = 0,                              # lower_bound = 0 ensures non-negative costs.
        optimizer = optimizer,
    ) do md, stage                                    # The do block defines the model for each stage (this part continues below).

#########-------------------------------------------------------------------------
######### Year and week for the current stage
#########-------------------------------------------------------------------------
        timenow = TimePoint(d.rundata.start_yr, d.rundata.start_wk) + stage - 1  # Constructs a TimePoint object for the current stage by 
                                                                                 # adding (stage - 1) weeks to the starting year and week.
        CONTINGENT = [
            r for r in s.RESERVOIRS if sum(                                      # Loops through all reservoirs in s.RESERVOIRS. For each reservoir r:
                d.reservoirs[r].contingent[timenow][j].level for                 # checks its contingent storage levels at timenow.
                j in 1:length(d.reservoirs[r].contingent[timenow])               # If the sum of levels is greater than zero, it’s considered active.
            ) > 0.0
        ]                                                                        # CONTINGENT becomes a list of reservoirs that have non-zero contingent water.rundata

        dr_keys = [
            (n, bl) for n in keys(d.dr_tranches[timenow]), bl in s.BLOCKS if     # Builds a list of (node, block) pairs where demand response tranches are defined.
            bl in keys(d.dr_tranches[timenow][n])
        ]

        en_keys = unique([
            (n, sector, lb) for n in keys(d.en_tranches[timenow]) for            # Builds a list of (node, sector, loadblock) tuples for energy tranches.
            (sector, lb) in keys(d.en_tranches[timenow][n])
        ])
#########------------------------------------------------------------------------

#########------------------------------------------------------------------------
######### State variable: water in reservoirs
#########------------------------------------------------------------------------
        JuMP.@variable(
            md,
            -sum( d.reservoirs[r].contingent[timenow][j].level / scale_factor for       # Calculates the total contingent water in reservoir r at timenow.
                  j in 1:length(d.reservoirs[r].contingent[timenow]) ) / scale_factor   # Scales it down using scale_factor.
            <= reslevel[r in s.RESERVOIRS] <=                                           # a state variable representing the water level in reservoir r at the current stage.
            d.reservoirs[r].capacity[timenow] / scale_factor,                           # The maximum water level allowed in the reservoir at timenow, scaled appropriately.
            SDDP.State,                                                                 # keyword SDDP.State provided by the SDDP.jl package. It tells the model that the variable being declared is a state variable in a multi-stage stochastic optimization problem.
            initial_value = d.reservoirs[r].initial / scale_factor                      # Sets the starting water level for each reservoir.
        )
#########------------------------------------------------------------------------

#########------------------------------------------------------------------------
######### Other variables
#########------------------------------------------------------------------------
        FLOWOVER = [a for a in s.NATURAL_ARCS if d.natural_arcs[a].maxflow != Inf]      # Get list of natural (hydro) arcs that have maximum flow limits.
        FLOWUNDER = [a for a in s.NATURAL_ARCS if d.natural_arcs[a].minflow != 0.0]     # Get list of natural (hydro) arcs that have minimum flow limits > 0.
        SPILLOVER = [a for a in s.STATION_ARCS if d.station_arcs[a].maxflow != Inf]     # Get list of station arcs that have maximum spill flow limits.

        JuMP.@variables(
            md,
            begin
                hydro_disp[s.HYDROS, s.BLOCKS] >= 0              # Dispatch of energy in MW from hydro stations
                
                thermal_use[s.THERMALS, s.BLOCKS] >= 0           # Amount of thermal energy used, in MW

                transflow[s.TRANS_ARCS, s.BLOCKS]                # Transmission flows between nodes in MW
                
                naturalflows[s.NATURAL_ARCS, s.BLOCKS] >= 0      # Water flows in cumecs
                
                releases[s.STATION_ARCS, s.BLOCKS] >= 0          # Water going through hydro stations
                
                spills[s.STATION_ARCS, s.BLOCKS] >= 0            # Water spilled from hydro stations
                
                terminalcost >= 0                                # Residual cost the final time period
                
                0 <= lostload[                                   # Lost load amount in MW
                    n in keys(d.dr_tranches[timenow]),
                    bl in keys(d.dr_tranches[timenow][n]),
                    k in keys(d.dr_tranches[timenow][n][bl]),
                ] <= d.dr_tranches[timenow][n][bl][k].q
                
                0 <= energyshedding[                              # Energy shedding
                    n in keys(d.en_tranches[timenow]),
                    s in keys(d.en_tranches[timenow][n]),
                    k in 1:length(d.en_tranches[timenow][n][s]),
                ] <= d.en_tranches[timenow][n][s][k].q * 1E3

                flowover[FLOWOVER, s.BLOCKS] >= 0                 # Any flows over upper bound (cumecs)
                
                flowunder[FLOWUNDER, s.BLOCKS] >= 0               # Any flows below lower bound (cumecs)
                
                spillover[SPILLOVER, s.BLOCKS] >= 0               # Spill flows over upper bound (cumecs)
                
                contingent[                                       # Contingent storage tranche
                    r in CONTINGENT,
                    1:length(d.reservoirs[r].contingent[timenow]),
                ] >= 0
                
                inflow[[s.CATCHMENTS_WITH_INFLOW; [:scenario]]]   # To track inflow levels seen
            end
        )

        if d.rundata.losses != :none  # This loss mode is not currently used by EA                   
            JuMP.@variables(
                md,
                begin
                    postransflow[s.TRANS_ARCS, s.BLOCKS] >= 0      # Positive partof transflow
                    
                    negtransflow[s.TRANS_ARCS, s.BLOCKS] >= 0      # Negative part of transflow
                    
                    postransflowtranche[                           # tranches of transmission positive flow for losses
                        (i, j) in s.TRANS_ARCS,
                        1:length(d.transmission[(i, j)].Poslosses),
                        s.BLOCKS ] >= 0

                    negtransflowtranche[                            # tranches of transmission negative flow for losses
                        (i, j) in s.TRANS_ARCS,
                        1:length(d.transmission[(i, j)].Neglosses),
                        s.BLOCKS ] >= 0
                    
                    losses[s.TRANS_ARCS, s.BLOCKS] >= 0             # losses in MW on trnamission line

                    node_losses[s.NODES, s.BLOCKS] >= 0             # Transmission loss allocated to node
                end
            )
        end
#########------------------------------------------------------------------------

#########------------------------------------------------------------------------
######### Define handy expressions
#########------------------------------------------------------------------------

        JuMP.@expressions(   # JuMP.@expressions(md, begin ... end)   A macro that defines multiple expressions at once in the JuMP model md.
            md,
            begin
                
                totHours, sum(d.durations[timenow][bl] for bl in s.BLOCKS)          # Number of hours in a week
             
                transmission[n in s.NODES, bl in s.BLOCKS],                         # Net transmission to any node: transmission to, minus transmission away
                sum(transflow[(i, j), bl] for (i, j) in s.TRANS_ARCS if j == n) -
                sum(transflow[(i, j), bl] for (i, j) in s.TRANS_ARCS if i == n)

                
                flowpenalties,   # Penalties for going over/under flow bounds. Note spMax is in MWh/m^3. 
                SECONDSPERHOUR * sum(
                    (   # Using general UB/LB flow penalty from run.csv if no upper/lower-bound penalty is provided for natural hydro arc and/or hydro station
                        sum( penalty_lb * d.spMax * flowunder[a, bl] for a in FLOWUNDER if d.natural_arcs[a].lb_penalty == -1.0 ) +
                        sum( penalty_ub * d.spMax * flowover[a, bl]  for a in FLOWOVER  if d.natural_arcs[a].ub_penalty == -1.0 ) +
                        sum( penalty_ub * d.spMax * spillover[a, bl] for a in SPILLOVER if d.station_arcs[a].penalty    == -1.0 ) +
                        
                        # Using specific UB/LB flow penalty provided for natural hydro arc and/or hydro station
                        sum( d.natural_arcs[a].lb_penalty * flowunder[a, bl] for a in FLOWUNDER if d.natural_arcs[a].lb_penalty > 0 )/1000 +
                        sum( d.natural_arcs[a].ub_penalty * flowover[a, bl]  for a in FLOWOVER  if d.natural_arcs[a].ub_penalty > 0 )/1000 +
                        sum( d.station_arcs[a].penalty    * spillover[a, bl] for a in SPILLOVER if d.station_arcs[a].penalty    > 0 )/1000    
                    ) * d.durations[timenow][bl] for bl in s.BLOCKS
                )

                # Water flow in minus flow out at any hydro hub
                netflow[n in s.CATCHMENTS, bl in s.BLOCKS],
                sum(naturalflows[(i, j), bl] for (i, j) in s.NATURAL_ARCS if j == n) - sum(naturalflows[(i, j), bl] for (i, j) in s.NATURAL_ARCS if i == n) +
                sum(releases[(i, j), bl] for (i, j) in s.STATION_ARCS if j == n) - sum(releases[(i, j), bl] for (i, j) in s.STATION_ARCS if i == n) +
                sum(spills[(i, j), bl] for (i, j) in s.STATION_ARCS if j == n) - sum(spills[(i, j), bl] for (i, j) in s.STATION_ARCS if i == n)
            end
        )

        # Total supply of electricity at any node and block
        if d.rundata.losses != :none
            JuMP.@expression(
                md,
                supply[n in s.NODES, bl in s.BLOCKS],
                d.durations[timenow][bl] * (
                    transmission[n, bl] - node_losses[n, bl] +
                    sum(thermal_use[m, bl] for m in d.nodehas[n].thermal) +
                    sum(hydro_disp[m, bl]  for m in d.nodehas[n].hydro)
                )
            )
        else
            JuMP.@expression(
                md,
                supply[n in s.NODES, bl in s.BLOCKS],
                d.durations[timenow][bl] * (
                    transmission[n, bl] +
                    sum(thermal_use[m, bl] for m in d.nodehas[n].thermal) +
                    sum(hydro_disp[m, bl]  for m in d.nodehas[n].hydro)
                )
            )
        end
#########------------------------------------------------------------------------

#########------------------------------------------------------------------------
######### Define constraints
#########------------------------------------------------------------------------
        JuMP.@constraints(
            md,
            begin
                # hydro/water flow lower and upper bounds constraint
                natOver[a in FLOWOVER, bl in s.BLOCKS]   , flowover[a, bl]  >= naturalflows[a, bl] - d.natural_arcs[a].maxflow
                natUnder[a in FLOWUNDER, bl in s.BLOCKS] , flowunder[a, bl] >= d.natural_arcs[a].minflow - naturalflows[a, bl]
                spillOver[a in SPILLOVER, bl in s.BLOCKS], spillover[a, bl] >= spills[a, bl] - d.station_arcs[a].maxflow

                # Capacity constraints

                # Hydro plant capacities
                maxHydroGeneration[m in s.HYDROS, bl in s.BLOCKS],
                hydro_disp[m, bl] <= d.hydro_stations[m].capacity - get(d.outage[timenow], (m, bl), 0.0)

                # Thermal plant capacities
                maxThermalGeneration[m in s.THERMALS, bl in s.BLOCKS],
                thermal_use[m, bl] <= d.thermal_stations[m].capacity - get(d.outage[timenow], (m, bl), 0.0)

                # Transmission line capacities
                transUpper[(n, m) in s.TRANS_ARCS, bl in s.BLOCKS],
                transflow[(n, m), bl] <= d.transmission[(n, m)].poscapacity - d.transmission[(n, m)].posoutage[timenow][bl]

                transLower[(n, m) in s.TRANS_ARCS, bl in s.BLOCKS],
                -d.transmission[(n, m)].negcapacity + d.transmission[(n, m)].negoutage[timenow][bl] <= transflow[(n, m), bl]

                # Set thermal station generation to zero if the station has not yet been commissioned, or is decommissioned
                notCommissioned[t in s.THERMALS, bl in s.BLOCKS; timenow < d.thermal_stations[t].commission], thermal_use[t, bl] <= 0
                decommissioned[t in s.THERMALS, bl in s.BLOCKS;  timenow > d.thermal_stations[t].decommission 
                                                              && d.thermal_stations[t].decommission != TimePoint(0, 0)  # Ignore this constraint if thermal station has no decommissioned date
                                                              ], thermal_use[t, bl] <= 0

                # All flow through hydro-stations is used in hydro generation
                hydroDispatch[m in s.HYDROS, bl in s.BLOCKS], hydro_disp[m, bl] == d.hydro_stations[m].sp * releases[d.hydro_stations[m].arc, bl]

                # Define shedding: the load shed over all sectors has to equal the shortage. In MWh.
                # Total energy shed (in MWh) at node n and block bl is at least equal to the shortage between demand and supply.
                defineShedding[n in s.NODES, bl in s.BLOCKS],
                sum(lostload[n, bl, k] for k in keys(d.dr_tranches[timenow][n][bl])) * d.durations[timenow][bl] >=
                d.demand[timenow][(n, bl)] - d.fixed[timenow][(n, bl)] - supply[n, bl]

                # The total energy shed across all relevant blocks and tranches for a sector must be 
                # less than or equal to the energy shedding capacity available for that sector and node.
                energyShedding[(n, sector, loadblocks) in en_keys],
                sum(lostload[n, bl, (s, name)] * d.durations[timenow][bl] for bl in loadblocks, (s, name) in keys(d.dr_tranches[timenow][n][bl]) if s == sector
                ) <= sum(energyshedding[n, (sector, loadblocks), k] for k in 1:length(d.en_tranches[timenow][n][(sector, loadblocks)]))

            end
        )

        if length(CONTINGENT) != 0
            JuMP.@constraints(
                md,
                begin
                    # Output level of the reservoir is at least the negative of the total contingent usage, scaled appropriately ???
                    contingentstorage[r in CONTINGENT],
                    reslevel[r].out >= -sum( contingent[r, j] / scale_factor for j in 1:length(d.reservoirs[r].contingent[timenow]) )

                    # The amount used from each contingent tranche does not exceed its available level.
                    maxcontingenttranche[r in CONTINGENT, j in 1:(length(d.reservoirs[r].contingent[timenow])-1)],
                    contingent[r, j] <= d.reservoirs[r].contingent[timenow][j].level
                end
            )
        end
#########------------------------------------------------------------------------

#########------------------------------------------------------------------------
######### ABP losses code
#########------------------------------------------------------------------------
        if d.rundata.losses != :none
            JuMP.@constraints(
                md,
                begin
                    # Define total positive flow on an arc.
                    definePosFlowTotal[(i, j) in s.TRANS_ARCS, bl in s.BLOCKS],
                    postransflow[(i, j), bl] == sum(postransflowtranche[(i, j), k, bl] for k in 1:length(d.transmission[(i, j)].Poslosses))

                    # Define total negative flow on an arc.
                    defineNegFlowTotal[(i, j) in s.TRANS_ARCS, bl in s.BLOCKS],
                    negtransflow[(i, j), bl] == sum(negtransflowtranche[(i, j), k, bl] for k in 1:length(d.transmission[(i, j)].Neglosses))
                    
                    # Tranche Flow Limits (Positive)
                    definePosLossTranche[(i, j) in s.TRANS_ARCS, bl in s.BLOCKS, k in 1:length(d.transmission[(i, j)].Poslosses)],
                    postransflowtranche[(i, j), k, bl] <= d.transmission[(i, j)].Poslosses[k][1]

                    # Tranche Flow Limits (Negative)
                    defineNegLossTranche[(i, j) in s.TRANS_ARCS, bl in s.BLOCKS, k in 1:length(d.transmission[(i, j)].Neglosses)],
                    negtransflowtranche[(i, j), k, bl] <= d.transmission[(i, j)].Neglosses[k][1]
                    
                    # Total losses on transmission line (Arc)
                    defineArcLosses[(i, j) in s.TRANS_ARCS, bl in s.BLOCKS],
                    losses[(i, j), bl] >=
                    sum(postransflowtranche[(i, j), k, bl] * d.transmission[(i, j)].Poslosses[k][2] for k in 1:length(d.transmission[(i, j)].Poslosses)) + 
                    sum(negtransflowtranche[(i, j), k, bl] * d.transmission[(i, j)].Neglosses[k][2] for k in 1:length(d.transmission[(i, j)].Neglosses))

                    # Net Flow Calculation
                    definePosFlow[(i, j) in s.TRANS_ARCS, bl in s.BLOCKS],
                    transflow[(i, j), bl] == postransflow[(i, j), bl] - negtransflow[(i, j), bl]
                end
            )

            # Half line losses are to be added to load at each end of the line
            JuMP.@constraints(
                md,
                begin
                    defineNodalLosses[n in s.NODES, bl in s.BLOCKS],
                    node_losses[n, bl] ==
                    0.5 * sum(losses[(i, j), bl] for (i, j) in s.TRANS_ARCS if j == n) +
                    0.5 * sum(losses[(i, j), bl] for (i, j) in s.TRANS_ARCS if i == n)
                end
            )
        end
#########------------------------------------------------------------------------

#########------------------------------------------------------------------------
######### Fuel variables and constraints
#########------------------------------------------------------------------------
    if ApplyFuelConstraints

        # ------------------------------------------------------------------------
        # State variable: fuel storage
        # ------------------------------------------------------------------------
        JuMP.@variable(
            md,
            0 <= fuelstoragelevel[stg in s.STORED_FUELS] <=                             # a state variable representing the fuel (PJ) in storage stg at the current stage.
            d.fuel_storages[stg].capacity,                  
            SDDP.State,                                                                 # keyword SDDP.State provided by the SDDP.jl package. It tells the model that the variable being declared is a state variable in a multi-stage stochastic optimization problem.
            initial_value = d.fuel_storages[stg].initial                                # Sets the starting fuel level for each storage.
        )

        # 
        JuMP.@variables(
            md,
            begin                
                fuel_contract[s.STORED_FUELS] >= 0               # Amount of fuel delivered according to fuel contract, in TJ (ignore fuel contract for now)
                fuel_injection[s.STORED_FUELS] >= 0              # Amount of fuel injected into storage, in TJ (ignore max injection constraints)
                fuel_withdrawal[s.STORED_FUELS] >= 0             # Amount of fuel withdrawn from storage, in TJ (ignore max withdrawal constraints) 
            end
        )

        # Fuel (TJ) used for thermal generation
        JuMP.@expressions(
            md,
            begin
                fuel_use_TJ[sf in s.STORED_FUELS],
                sum(thermal_use[m,bl] * d.thermal_stations[m].heatrate * d.durations[timenow][bl] * 1e-3
                for bl in s.BLOCKS, m in keys(d.thermal_to_storage) if d.thermal_to_storage[m]==sf)
            end
        ) 

        JuMP.@constraints(
            md,
            begin
                # Fuel use/storage balance constraints (TJ)
                fuelUseBalance[sf in s.STORED_FUELS],
                fuel_contract[sf] + fuel_withdrawal[sf] - fuel_injection[sf] >= fuel_use_TJ[sf]

                # Fuel contracts(TJ)
                fuelContractDelivery[sf in s.STORED_FUELS],
                d.fuel_contracts[timenow][sf].mindelivery <= fuel_contract[sf] <= d.fuel_contracts[timenow][sf].maxdelivery

                # Conservation for fuel storages (TJ)
                fuelStorageBalance[sf in s.STORED_FUELS],
                fuelstoragelevel[sf].out - fuelstoragelevel[sf].in == (fuel_injection[sf] - fuel_withdrawal[sf])*1e-3
            end
        )
    end
#########------------------------------------------------------------------------

#########------------------------------------------------------------------------
######### Flow-conservation-related calculations
#########------------------------------------------------------------------------
        # Initialize Scenario List - Creates an empty array of dictionaries. 
        # Each dictionary will represent one inflow scenario, mapping catchments (and a scenario ID) to inflow values.
        inflow_uncertainty = Array{Dict{Symbol,Float64},1}()

        for scenario in 1:length(d.inflow_mat[1][s.CATCHMENTS_WITH_INFLOW[1]])   # Loops over all inflow scenarios available for the current week and catchments.

            s_inflows = Dict{Symbol,Float64}()
            for c in s.CATCHMENTS_WITH_INFLOW           # For each catchment c,
                s_inflows[c] = d.inflow_mat[timenow.week][c][scenario] # stores inflow value for the current week and scenario into s_inflows
            end

            if stage == 1 && d.rundata.first_week_known # If it's the first stage and the first week's inflow is known
                s_inflows[:scenario] = 0                # it sets the scenario ID to 0 (i.e., deterministic).
            else
                s_inflows[:scenario] = scenario         # Otherwise, it uses the actual scenario number. 
            end

            push!(inflow_uncertainty, s_inflows)        # Adds the scenario dictionary to the list of inflow scenarios.

        end

        # Parameterize the Model - This tells SDDP to:
        SDDP.parameterize(md, inflow_uncertainty) do ϕ  # Use inflow_uncertainty as the set of stochastic scenarios.
            for (c, value) in ϕ                         # For each scenario ϕ, fix the inflow values in the model to those specified in the dictionary.
                JuMP.fix(inflow[c], value)              
            end
        end

        JuMP.@constraints(
            md,
            begin
                # Conservation for reservoirs
                rbalance[r in s.RESERVOIRS],
                (reslevel[r].out - reslevel[r].in) * 1E3 * scale_factor ==
                SECONDSPERHOUR * 1E-3 * (sum(d.durations[timenow][bl] * (netflow[r, bl]) for bl in s.BLOCKS) + totHours * inflow[r])

                # Conservation for junction points with inflow
                jbalance[c in s.CATCHMENTS_WITH_INFLOW, bl in s.BLOCKS; c in s.JUNCTIONS],
                netflow[c, bl] + inflow[c] == 0

                # Flow conservation for junctions without an inflow
                conserveFlow[c in s.JUNCTIONS_WITHOUT_INFLOW, bl in s.BLOCKS],
                netflow[c, bl] == 0
            end
        )

        for dr in d.rundata.decision_rules  # Iterates over each decision rule dr defined in the model data (d.rundata.decision_rules) - not used for EA data
            if timenow.week ∉ dr.weeks
                continue    # Skips the rule if the current week (timenow.week) is not included in the rule’s applicable weeks (dr.weeks).
            end

            # This kind of logic is used in Hydro Scheduling Models to
            # 1. Ensuring hydro stations operate within environmental or operational limits.
            # 2. Balancing generation vs. reservoir levels and inflows.
            LHS = 0.0
            if dr.flowtype == :generation # Calculates energy generated by water releases.
                LHS = d.hydro_stations[dr.station].sp * sum( releases[d.hydro_stations[dr.station].arc, bl] * d.durations[timenow][bl] for bl in s.BLOCKS )
            
            elseif dr.flowtype == :spill # Calculates energy equivalent of water spilled (not used for generation).
                LHS = d.hydro_stations[dr.station].sp * sum( spills[d.hydro_stations[dr.station].arc, bl] * d.durations[timenow][bl] for bl in s.BLOCKS )
            
            elseif dr.flowtype == :combined # Calculates energy equivalent of both water releases and spills combined.
                LHS = d.hydro_stations[dr.station].sp * sum(( releases[d.hydro_stations[dr.station].arc, bl] + 
                                                              spills[d.hydro_stations[dr.station].arc, bl]) * d.durations[timenow][bl] for bl in s.BLOCKS )
            else
                error("Invalid flow type: " * string(dr.flowtype))
            end

            # Linear decision rule.
            RHS = dr.intercept + dr.slope * ( reslevel[dr.reservoir].in * scale_factor + SECONDSPERHOUR * totHours * inflow[dr.reservoir] / 1E6 )

            # Adds a constraint to the optimization model md.
            if dr.boundtype == :upper
                JuMP.@constraint(md, LHS <= RHS)
            elseif dr.boundtype == :lower
                JuMP.@constraint(md, LHS >= RHS)
            elseif dr.boundtype == :equality
                JuMP.@constraint(md, LHS == RHS)
            else
                error("Invalid bound type: " * string(dr.boundtype))
            end
        end
#########------------------------------------------------------------------------

#########------------------------------------------------------------------------
######### Objective-related calculations
#########------------------------------------------------------------------------
        JuMP.@expression( # Calulate total cost of lost load and energy shedding
            md,
            lostloadcosts,
            sum(
                sum(
                    lostload[n, bl, k] * d.durations[timenow][bl] * d.dr_tranches[timenow][n][bl][k].p
                    for k in keys(d.dr_tranches[timenow][n][bl])
                ) for (n, bl) in dr_keys
            ) + 
            sum(
                sum(
                    energyshedding[n, (sector, loadblocks), k] * d.en_tranches[timenow][n][(sector, loadblocks)][k].p
                    for k in 1:length(d.en_tranches[timenow][n][(sector, loadblocks)])
                ) for (n, sector, loadblocks) in en_keys
            )
        )

        JuMP.@expression( # Total cost incurred due to contingent storage requirements across multiple reservoirs and time slices
            md,           # To penalize deviations from minimum storage levels and incentivize maintaining reliability reserves in reservoirs.
            contingent_storage_cost,
            sum(
                contingent[r, j] / scale_factor *
                d.reservoirs[r].contingent[timenow][j].penalty for r in CONTINGENT,
                j in 1:length(d.reservoirs[r].contingent[timenow])
            )
        )

        JuMP.@expression( # Calculate total CO2 emissions for each thermal station t and time block bl
            md,
            carbon_emissions[t in s.THERMALS, bl in s.BLOCKS],
            d.carbon_content[d.thermal_stations[t].fuel] *
            d.thermal_stations[t].heatrate *
            thermal_use[t, bl] *
            d.durations[timenow][bl]
        )

        JuMP.@expression( # Calculate total thermal/hydro generation cost and pentalies
            md,
            immediate_cost,
            sum( # Thermal generation cost
                (station.omcost + d.fuel_costs[timenow][station.fuel] * station.heatrate) * thermal_use[name, bl] * d.durations[timenow][bl] +
                carbon_emissions[name, bl] * d.fuel_costs[timenow][:CO2] 
                for (name, station) in d.thermal_stations, bl in s.BLOCKS
            ) +
            sum( # Hydro generation cost
                station.omcost * hydro_disp[name, bl] * d.durations[timenow][bl] 
                for (name, station) in d.hydro_stations, bl in s.BLOCKS
            ) +
            flowpenalties +
            lostloadcosts +
            contingent_storage_cost
        )

        if stage < number_of_wks || !d.rundata.use_terminal_mwvs
            # Stage cost function not including terminal water value
            SDDP.@stageobjective(md, immediate_cost / scale_obj)
        else
            # Convert stored water in Mm³ to MWh
            JuMP.@expression(
                md,
                storedenergy,
                1E6 * scale_factor * sum(d.reservoirs[r].sp * reslevel[r].out for r in s.RESERVOIRS)
            )

            for cut in d.terminal_eqns # -terminalcost for value
                JuMP.@constraint(
                    md,
                    terminalcost >= -(cut.intercept + cut.coefficient * storedenergy) / scale_obj
                )
            end
            # Cost function includes terminal values added
            SDDP.@stageobjective(md, immediate_cost / scale_obj + terminalcost)
        end
#########------------------------------------------------------------------------

    end

    return sddpm
end
