#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

#--------------------------------------------------
# Prepare reservoir-related data
#--------------------------------------------------

struct ContingentTranche
    level::Float64
    penalty::Float64
end

mutable struct Reservoir
    capacity::TimeSeries{Float64}                      # consented capacity (Mil m3) of the reservoir may vary throughout the year.
    initial::Float64                                   # reservoir's storage level (Mil m3) at the start of study period
    sp::Float64                                        # reservoir's specific power factor (kWh/m3) is determined based on the SP values of downstream hydroelectric stations.
    contingent::TimeSeries{Vector{ContingentTranche}}  # avaialble contingent storage and "penalty" cost to use it may vary throughtout the year
    index::Int                                         # index assigned to a reservoir
end

"""
    initialisereservoirs(file::String, limits::String)

    Description:
    Initializes a dictionary of Reservoir objects using metadata and time series data from two CSV files:
    1. One for basic reservoir info (reservoirs.csv)
    2. One for time-varying limits and contingent storage (reservoir_limits.csv)

    Inputs:
    - reservoirs_filename::String: Path to the CSV file containing reservoir metadata.
    - reservoir_limits_filename::String: Path to the CSV file containing time series data for reservoir limits and contingent storage.

    Outputs:
    reservoirs::Dict{Symbol, Reservoir}: A dictionary mapping reservoir names to their corresponding Reservoir objects.
"""
function initialisereservoirs(reservoirs_filename::String,reservoir_limits_filename::String)
    reservoirs = Dict{Symbol,Reservoir}()  # Create an empty dictionary to store reservoirs

        # Load reservoir metadata from "reservoirs.csv"
    for row in CSV.Rows(reservoirs_filename; missingstring = "NA", stripwhitespace = true, comment = "%")

        # Check if the format of "reservoir.csv" is valid
        row = _validate_and_strip_trailing_comment(row, [:RESERVOIR, :INFLOW_REGION, :INI_STATE], [:CAPACITY])

        # Converts the reservoir name to a Symbol and checks for duplicates.
        reservoir = str2sym(row.RESERVOIR)
        if haskey(reservoirs, reservoir)
            error("Reservoir $(reservoir) given twice.")
        end

        reservoirs[reservoir] = Reservoir(                               # Initializes a Reservoir object with:
            TimeSeries{Float64}(TimePoint(0, 0), []),                    # placeholder for capacity time series.
            parse(Float64, row.INI_STATE),                               # get initial storage level.
            0.0,                                                         # specific power set to 0 for now to be replaced later.
            TimeSeries{Vector{ContingentTranche}}(TimePoint(0, 0), []),  # placeholder for contingent time series.
            length(reservoirs) + 1,                                      # index of reservoir, used to create compatible DOASA cut files.
        )
    end

    # Load reservoir limits time series
    limits, column_names = gettimeseries(reservoir_limits_filename)

    # For each reservoir
    for (name, res) in reservoirs 
        # Set capacity time series
        res.capacity = TimeSeries{Float64}(limits.startpoint, [
            limits[i][Symbol("$name MAX_LEVEL")] for i in 1:length(limits)
        ])

        # Initialise contingent storage tranches
        tranches = [ContingentTranche[] for i in 1:length(limits)]

        # Initialise total contingent limit for reservoir for later validation check
        total = zeros(Float64, length(limits))

        # If no contingent levels 1, add default tranche
        if !(Symbol("$name MIN_1_LEVEL") in column_names)
            for i in 1:length(limits)
                push!(tranches[i], ContingentTranche(0.0, 0.0))
            end
        end

        # Process contingent levels
        j = 1
        while Symbol("$name MIN_$(j)_LEVEL") in column_names
            penalty_col = Symbol("$name MIN_$(j)_PENALTY")
            level_col = Symbol("$name MIN_$(j)_LEVEL")
            prev_level_col = Symbol("$name MIN_$(j-1)_LEVEL")

            if !(penalty_col in column_names)
                error("$name contingent storage missing MIN_$(j)_PENALTY")
            end
            # The logic in this loop is questionable - Tuong Nguyen
            for i in 1:length(limits)
                total[i] = -limits[i][level_col] + (j > 1 ? limits[i][prev_level_col] : 0.0)
                push!(tranches[i], ContingentTranche(total[i], limits[i][penalty_col]))
            end
            j += 1
        end

        if maximum(total) - minimum(total) >= 1e-8
            @warn("The maximum contingent storage is not constant for reservoir $name.")
        end

        res.contingent = TimeSeries{Vector{ContingentTranche}}(limits.startpoint, tranches)
    end

    return reservoirs
end

#------------------------------------------------------
# Flow arcs and hydro generator data
#------------------------------------------------------
# Arcs in water network
struct NaturalArc
    minflow::Float64        # flow lower bound
    maxflow::Float64        # flow upper bound
    lb_penalty::Float64     # flow lower bound violation penalty
    ub_penalty::Float64     # flow upper bound violation penalty
end

"""
    getnaturalarcs(file::String)

    # Description:
    This function reads a CSV file (hydro_arcs.csv) containing natural water flow connections (like rivers or canals) between hydro stations.
    It returns a dictionary mapping each origin-destination pair to a NaturalArc object, which holds flow constraints and penalties (optional).

    # Inputs:
    filename::String: The path to the CSV file that contains the natural arc data.

    # Outputs:
    - arcs::Dict{NTuple{2,Symbol}, NaturalArc}: A dictionary where:
        Key: A tuple of two Symbols representing the origin and destination nodes.
        Value: A NaturalArc object containing flow constraints and penalties.

    # Example File (The columns MUST be ordered as shown below.)
    ORIG,DEST,MIN_FLOW,MAX_FLOW
    Clyde_tail,Lake_Roxburgh,66.25,na
    Karapiro_tail,SEA,140.0,550.0

"""
function getnaturalarcs(filename::String)
    arcs = Dict{NTuple{2,Symbol},NaturalArc}()

    for row in CSV.Rows(filename; missingstring = ["NA", "na", "default"], stripwhitespace = true, comment = "%")

        row = _validate_and_strip_trailing_comment(row, [:ORIG, :DEST, :MIN_FLOW, :MAX_FLOW], [:LB_PENALTY, :UB_PENALTY])

        key = (str2sym(row.ORIG), str2sym(row.DEST))
        if haskey(arcs, key)
            error("Arc $(key) given twice.")
        end

        minflow = parse(Float64, coalesce(row.MIN_FLOW, "0.0"))
        maxflow = parse(Float64, coalesce(row.MAX_FLOW, "Inf"))
        lb_penalty = parse(Float64, coalesce(get(row, :LB_PENALTY, missing), "-1"))
        ub_penalty = parse(Float64, coalesce(get(row, :UB_PENALTY, missing), "-1"))

        arcs[key] = NaturalArc(minflow, maxflow, lb_penalty, ub_penalty)
    end

    return arcs
end

struct StationArc
    maxflow::Float64    # Spillway max flow in cumecs (can be zero if no spillway exists, or "na" for unlimited spill).
    penalty::Float64    # Penalty apply to over-spill
    station::Symbol     # Name of the hydro station
end

# Hydro generators
struct HydroStation
    node::Symbol
    capacity::Float64       # Capacity in MW
    sp::Float64             # Specific power in MW/cumec;
    omcost::Float64         # Operation and maintenance cost $/MWh
    arc::NTuple{2,Symbol}   # Mapping water flow arcs (from one node to another).
end

"""
    gethydros(filename::String, nodes::Vector{Symbol})

    # Description:
    This function reads hydroelectric station data from a CSV file (hydro_stations.csv) and constructs two dictionaries:
    1. hydros: Maps generator symbols to HydroStation objects.
    2. station_arcs: Maps water flow arcs (from one node to another) to StationArc objects.

    # Inputs:
    - filename::String: The path to the CSV file containing hydro station data.
    - nodes::Vector{Symbol}: A list of valid power system nodes used to validate entries in the CSV.

    # Outputs:
    - hydros::Dict{Symbol, HydroStation}: A dictionary mapping generator names to their corresponding hydro station data.
    - station_arcs::Dict{NTuple{2,Symbol}, StationArc}: A dictionary mapping water flow arcs to their associated station arc data.

    # Example File:
    GENERATOR,HEAD_WATER_FROM,TAIL_WATER_TO,POWER_SYSTEM_NODE,CAPACITY,SPECIFIC_POWER,SPILLWAY_MAX_FLOW
    Arapuni,Lake_Arapuni,Lake_Karapiro,NI,196.7,0.439847649,99999
"""
function gethydros(filename::String, nodes::Vector{Symbol})
    hydros = Dict{Symbol,HydroStation}()
    station_arcs = Dict{NTuple{2,Symbol},StationArc}()

    for row in CSV.Rows(filename; missingstring=["NA", "na", "default"], stripwhitespace=true, comment="%")
        row = _validate_and_strip_trailing_comment(
            row,
            [:GENERATOR, :HEAD_WATER_FROM, :TAIL_WATER_TO, :POWER_SYSTEM_NODE, :CAPACITY, :SPECIFIC_POWER, :SPILLWAY_MAX_FLOW],
            [:OM_COST, :OVERSPILL_PENALTY]
        )

        gen               = str2sym(row.GENERATOR)
        arc               = (str2sym(row.HEAD_WATER_FROM), str2sym(row.TAIL_WATER_TO))
        node              = str2sym(row.POWER_SYSTEM_NODE)
        capacity          = parse(Float64, row.CAPACITY)
        specific_power    = parse(Float64, row.SPECIFIC_POWER)
        max_spill_flow    = parse(Float64, coalesce(row.SPILLWAY_MAX_FLOW, "Inf"))
        om_cost           = parse(Float64, get(row, :OM_COST, "0.0"))
        overspill_penalty = parse(Float64, get(row, :OVERSPILL_PENALTY, "-1"))

        if haskey(hydros, gen)
            error("Generator $gen given twice.")
        elseif haskey(station_arcs, arc)
            error("Station arc $arc already given.")
        elseif !(node in nodes)
            error("Node $node for generator $gen not found.")
        end

        hydros[gen] = HydroStation(node,capacity,specific_power,om_cost,arc)

        station_arcs[arc] = StationArc(max_spill_flow,overspill_penalty,gen)
            
    end

    return hydros, station_arcs
end

#---------------------------------------------------
# Prepare inflow-related data
#---------------------------------------------------
"""
    adjustinflows(inflow_file::String, rundata::RunData)

    # Description:
    The adjustinflows function processes historical inflow data for multiple catchments 
    and applies a statistical transformation (DIA) to produce adjusted inflow time series. 
    This transformation is inspired by the methodology described in the DOASA model.

    # Inputs:
    - filename::String: The path to the CSV file containing inflow data (inflows.csv).
    - rundata:: An object containing the JADE model settings, required to define the stage problems in SDDP.jl.

    # Outputs:
    - inflows::Dict{Symbol,TimeSeries{Float64}}: A dictionary mapping catchment names to their corresponding "adjusted" inflow data.
    - first_inflows::Dict{Symbol,Float64}: A dictionary mapping catchment name to inflow of the start year and start week

"""
function adjustinflows(inflow_file::String, rundata::RunData)
    if rundata.dialength >= 52
        error("Correlation length too large.")
    end
    # correlation length ω
    ω = if rundata.dialength <= 0
        1
    else
        rundata.dialength
    end

    catchments = Symbol[]
    regions = Symbol[]
    tmp_inflows = Dict{Symbol,Dict{TimePoint,Float64}}()
    # wrangle data
    parsefile(inflow_file, true) do items
        if lowercase(items[1]) == "catchment"
            for it in items[3:end]
                push!(catchments, str2sym(it))
                tmp_inflows[str2sym(it)] = Dict{TimePoint,Float64}()
            end
        elseif lowercase(items[1]) == "inflow_region"
            for it in items[3:end]
                push!(regions, str2sym(it))
            end
            @assert length(regions) == length(catchments)
        elseif lowercase(items[1]) == "year"
            return
        else
            year = parse(Int, items[1])
            week = parse(Int, items[2])
            tp = TimePoint(year, week)
            for (i, catchment) in enumerate(catchments)
                tmp_inflows[catchment][tp] = parse(Float64, items[i+2])
            end
        end
    end
    inflows = Dict{Symbol,TimeSeries{Float64}}()
    first_inflows = Dict{Symbol,Float64}()

    for (location, val) in tmp_inflows
        tps, infl = collect(keys(val)), collect(values(val))
        # sort all our data
        p = sortperm(tps)
        permute!(tps, p)
        permute!(infl, p)

        first_inflows[location] = val[TimePoint(rundata.start_yr, rundata.start_wk)]

        weeks = sort(unique([tp.week for tp in tps]))
        years = sort(unique([tp.year for tp in tps]))
        N = length(years)
        # Follow the steps in DOASA paper
        α = [  # mean weekly historical inflow
            Statistics.mean(inflow for (tp, inflow) in zip(tps, infl) if tp.week == t) for t in weeks
        ]
        W = [] # rolling total inflow starting from week ty
        for ty in 1:length(infl)
            last = ((ty + ω - 2) % length(infl) + 1)
            if last < ty
                push!(W, sum(infl[ty:length(infl)]) + sum(infl[1:last]))
            else
                push!(W, sum(infl[ty:last]))
            end
        end
        m = [  # mean of the rolling inflow
            Statistics.mean(W[i] for (i, tp) in enumerate(tps) if tp.week == t) for t in weeks
        ]
        d = [  # mean + the bigger deviation
            max(0.0, α[tp.week] + (W[i] - m[tp.week]) / √ω) for (i, tp) in enumerate(tps)
        ]
        d_mean = [  # mean of adjusted values
            Statistics.mean(d[i] for (i, tp) in enumerate(tps) if tp.week == t) for t in weeks
        ]
        k = [  # adjusted inflows
            isnan(d[i] * α[tp.week] / d_mean[tp.week]) ? α[tp.week] :
            d[i] * α[tp.week] / d_mean[tp.week] for (i, tp) in enumerate(tps)
        ]
        inflows[location] = TimeSeries(minimum(tps), k)
    end
    return inflows, first_inflows
end

"""
    diatofile(adjusted::Dict{Symbol,TimeSeries{Float64}}, outpath::String, policy_dir::String)

    # Description:
    This function saves DIA adjusted inflows (created by the "adjustinflows" function) to `adjusted_inflows.csv`.

    # Required Arguments
    `adjusted` is the dictionary of DIA-adjusted inflows
    `outpath` is the filename for adjusted inflows
    `policy_dir` is the subdirectory where the outputs are stored
"""
function diatofile(adjusted::Dict{Symbol,TimeSeries{Float64}}, outpath::String, policy_dir::String)
    outpath = joinpath(@__JADE_DIR__, "Output", outpath)
    if !ispath(joinpath(outpath, policy_dir))
        mkpath(joinpath(outpath, policy_dir))
    end
    open(joinpath(outpath, policy_dir * "/adjusted_inflows.csv"), "w") do io
        # print header
        print(io, "YEAR, WEEK")
        for (name, location) in adjusted
            print(io, ", $(name)")
        end
        print(io, "\n")
        # print data
        l = first(keys(adjusted))
        starttime = adjusted[l].startpoint
        endtime = starttime + length(adjusted[l]) - 1
        for t in starttime:endtime
            print(io, t.year, ", ", t.week)
            for (name, location) in adjusted
                print(io, ", ", location[t])
            end
            print(io, "\n")
        end
    end
end

"""
    getinflows(adjusted::Dict{Symbol,TimeSeries{Float64}},firstweekinflows::Dict{Symbol,Float64},rundata::RunData)

    # Description:
    Constructs a weekly inflow dataset for multiple catchments across multiple scenarios, using historical data and optionally overriding the first week's inflow.

    # Input:
    - adjusted        :: A dictionary mapping catchment names (Symbol) to their adjusted inflow time series (TimeSeries{Float64}).
    - firstweekinflows:: A dictionary mapping catchments to their inflow values at the start of the simulation.
    - rundata         :: A struct containing metadata like sample years, number of scenarios, and whether the first week inflow is known.

    # Output:
    - An array of dictionaries. It's a 1D array with 52 elements (one for each week of the year).
    - Each element is a dictionary mapping catchment names (Symbol) to a vector of inflow values (Vector{Float64}), one per scenario.
"""
function getinflows(adjusted::Dict{Symbol,TimeSeries{Float64}}, firstweekinflows::Dict{Symbol,Float64}, rundata::RunData)
    # Creates an empty array of dictionaries. 
    # Each element of the array will correspond to a week of the year.
    # Each dictionary maps catchment names to vectors of inflow values (one per scenario).
    inflows = Dict{Symbol,Vector{Float64}}[]

    # Allocate storage
    for week in 1:WEEKSPERYEAR  # For each week (1 to 52), creates a dictionary.
        push!(inflows, Dict{Symbol,Vector{Float64}}())
        for l in keys(adjusted) # For each catchment l, initializes an empty vector to store inflows for that week.
            inflows[week][l] = Float64[]
        end
    end

    # Populate Inflows from Sample Years
    for year in rundata.sample_years
        # Defines the start and stop TimePoint for that year (week 1 to 52)
        samplestart = TimePoint(year, 1)
        samplestop = TimePoint(year, WEEKSPERYEAR)

        # Push data (in order)
        for (name, location) in adjusted
            for t in samplestart:samplestop
                push!(inflows[t.week][name], location[t])
            end
        end
    end

    # locations =
    @assert length(inflows[1][first(keys(adjusted))]) == rundata.nscenarios

    if rundata.first_week_known
        # In start week, only consider inflows from "problem start year"
        #   (DOASA feature)
        for (name, location) in adjusted
            inflows[rundata.start_wk][name] =
                repeat([firstweekinflows[name]], outer = rundata.nscenarios)
        end
    end

    return inflows
end

"""
    getinflows_for_historical(inflowsfile::String, rundata::RunData, year::Union{Int,Vector{Int}})
    This function is is a wrap around the "getinflows" function to get the actual histrorical inflow by setting dialenght = 1 (no DIA adjustment)
"""
function getinflows_for_historical(inflowsfile::String, rundata::RunData, year::Union{Int,Vector{Int}})
    rd = deepcopy(rundata)
    rd.dialength = 1
    if typeof(year) == Int
        rd.sample_years = [year]
        rd.nscenarios = 1
    else
        rd.sample_years = year
        rd.nscenarios = length(year)
    end
    adjusted_inflows, firstweekinflows = adjustinflows(inflowsfile, rd)
    return getinflows(adjusted_inflows, firstweekinflows, rd)
end


#---------------------------------------------------
# Function to get reservoirs' specific power factors
#---------------------------------------------------
"""
	out_neighbors(vertex::Symbol, edges::Vector{NTuple{2,Symbol}})
    
    # Purpose: 
    This function finds all outgoing neighbors of a given vertex in a directed graph.

    # Inputs:
    - vertex ::Symbol                  : The node whose outgoing neighbors you want to find.
    - edges  ::Vector{NTuple{2,Symbol}}: A list of directed edges, where each edge is a tuple (from, to).

    # Output:
    - Returns a list of all nodes that are directly reachable from vertex.
"""
function out_neighbors(vertex::Symbol, edges::Vector{NTuple{2,Symbol}})
    neighbors = Symbol[]
    for e in edges
        if e[1] == vertex
            push!(neighbors, e[2])
        end
    end
    return neighbors
end

"""
    hasdownstream(sets::Sets, station_arcs::Dict{NTuple{2,Symbol},StationArc})

    # Purpose:
    This function determines which hydro stations are downstream of each reservoir, based on a network of catchments and arcs.

    # Inputs:
    sets        ::Sets                              : A structure containing sets of catchments, reservoirs, arcs, etc.
    station_arcs::Dict{(Symbol, Symbol), StationArc}: Maps arcs to hydro station metadata.  

    # Output:
    Returns a dictionary ::Dict{Symbol, Vector{Symbol}} where each key is a reservoir, and the value is a list of hydro stations downstream from it.

"""
function hasdownstream(sets::Sets, station_arcs::Dict{NTuple{2,Symbol},StationArc})

    @assert !isempty(sets.RESERVOIRS)
    @assert !isempty(sets.STATION_ARCS)
    
    reverseflow = Symbol[] # Currently unused

    # Filter out reverse flow station arcs (no effect currently)
    station_arcs_filtered = Dict(filter(((pair, arc),) -> arc.station ∉ reverseflow, station_arcs))

    # All arcs put together
    arcs = union(sets.NATURAL_ARCS, sets.STATION_ARCS)

    # Remove arcs associated with reverse flow stations (no effect currently)
    for (pair, arc) in station_arcs
        if arc.station in reverseflow
            filter!(x -> x != pair, arcs)
        end
    end
    
    # Build downstream neighbor map
    neighbors = Dict(c => out_neighbors(c, arcs) for c in sets.CATCHMENTS)
    
    # Identify valid catchments
    catchments = [c for c in sets.CATCHMENTS if !isempty(neighbors[c]) || c == :SEA]
  
    # Identify valid reservoirs
    reservoirs = [r for r in sets.RESERVOIRS if (r in sets.CATCHMENTS && !isempty(neighbors[r])) || r ∉ sets.CATCHMENTS]

    # Dictionary structure to store which stations each reservoir has downstream
    reservoir_has_downstream = Dict{Symbol,Array{Symbol}}()

    # For every reservoir, we do a depth-first search (DFS) to look for hydro-stations downstream.
    for r in reservoirs

        # Creates a dictionary node_col to track the state of each catchment during DFS:
        # 0 = white (unvisited)
        # 1 = grey (visited but not fully explored)
        # 2 = black (fully explored)
        node_col = Dict(c => 0 for c in catchments)
        current = r            # Sets the current node to the reservoir r. This is the starting point for DFS.
        node_col[current] = 1  # Marks the starting reservoir as grey, indicating it's being visited.
        theList = Symbol[]     # Initializes an empty list to store hydro stations found downstream of this reservoir.
        grey_stack = [r]       # Initializes a stack for DFS traversal. Starts with the reservoir r.

        # Start/continue a depth first search (DFS) as long as there are nodes in the grey_stack.
        while !isempty(grey_stack)
            current = grey_stack[end]    # Gets the top of the stack (last element) — the current node being explored.

            if current == :SEA           # If the current node is the special terminal node :SEA:
                node_col[current] = 2    # Mark it as black (fully explored).
                pop!(grey_stack)         # Remove it from the stack.
                continue                 # Skip to the next "while" iteration.
            end

            moved_down = false           # Flag to track whether we moved deeper into the graph during this iteration.
            for n in neighbors[current]  # Loops over all downstream neighbors of the current catchment.
                thearc = (current, n)    # Defines the arc (edge) from current to neighbor n.

                if haskey(station_arcs_filtered, thearc)              # Checks if this arc corresponds to a hydro station.
                    station = station_arcs_filtered[thearc].station   # Gets the station associated with this arc.
                    if station ∉ theList && station ∉ reverseflow     # If the station hasn't been added yet and isn't in the reverse flow list (which is empty in this case)
                        push!(theList, station)                       # add it to theList.
                    end
                end

                if node_col[n] == 0       # If the neighbor n hasn't been visited yet:
                    node_col[n] = 1       # Mark it as grey.
                    push!(grey_stack, n)  # Push it onto the grey_stack.
                    moved_down = true     # Set moved_down = true to indicate we went deeper.
                    break                 # break out of for loop to explore this new node next.
                end
            end

            if !moved_down                # If we didn’t go deeper (i.e., all neighbors are visited):
                node_col[current] = 2     # Mark current as black (fully explored).
                pop!(grey_stack)          # Pop it from the stack to backtrack.
            end
        end
        reservoir_has_downstream[r] = theList # After DFS completes for reservoir r, store the list of downstream stations in the result dictionary.

    end # next reservoir

    return reservoir_has_downstream
end


"""
    set_reservoir_sp!(reservoirs::Dict{Symbol,Reservoir}, hydros::Dict{Symbol,HydroStation},
                      sets::Sets, station_arcs::Dict{NTuple{2,Symbol},StationArc})
    # Purpose:
    This function calculate the specific factors of reservoirs as the sum of downstrean stations' specific power factors

    # Inputs:
    reservoirs::Dict{Symbol,Reservoir}              : A dictionary mapping reservoir names to their corresponding Reservoir objects.
    hydros::Dict{Symbol,HydroStation}               : A dictionary mapping generator names to their corresponding hydro station data.
    sets        ::Sets                              : A structure containing sets of catchments, reservoirs, arcs, etc.
    station_arcs::Dict{(Symbol, Symbol), StationArc}: Maps arcs to hydro station metadata.  

    # Output:
    Returns the dictionary reservoirs::Dict{Symbol,Reservoir} with updated specific power factor
    
"""
function set_reservoir_sp!(reservoirs::Dict{Symbol,Reservoir}, hydros::Dict{Symbol,HydroStation},
                           sets::Sets, station_arcs::Dict{NTuple{2,Symbol},StationArc})
                           
    # Get downstream hydro stations for each reservoir
    downstream_stations = hasdownstream(sets, station_arcs)

    # Update the specific power of each reservoir:
    for (res_name, reservoir) in reservoirs
        if haskey(downstream_stations, res_name)
            total_sp = sum(hydros[station].sp for station in downstream_stations[res_name])
            reservoir.sp = total_sp / 3600   # Convert from MW/cumec to MWh/m3
        else
            reservoir.sp = 0
        end
    end
end


#---------------------------------------------------
# Get terminal water value
#---------------------------------------------------

# Defines a simple structure to represent a linear equation of the form: y = intercept + coefficient * x
struct LinearEquation
    intercept::Float64
    coefficient::Float64
end

"""
    getterminalvalue(file::String)

    # Description:
    This function reads data from a CSV file (terminal_water_value.csv) and returns a list of LinearEquation objects represent linear equations.
    - Value of stored energy in hydro lakes at the end of the time horizon.
    - STORED_ENERGY is the cumulative stored energy in GWh.
    - VALUE is the marginal value in \$/MWh. This column should be a decreasing sequence.

    The columns MUST be ordered as shown below.
        STORED_ENERGY, VALUE
        1000, 137.218398630355   % first 1000 GWh is worth ~ \$137/MWh
        1500, 85.9718321526058

    # Inputs: 
    filename::String: The path to the CSV file containing terminal water value data.
    
    # Outputs: 
    A list of LinearEquation objects.
"""
function getterminalvalue(file::String)
    equations = LinearEquation[]  # stores the resulting linear segments.
    cumulative_value = 0.0        # tracks the accumulated value across energy segments.
    previous_energy  = 0.0        # stores the energy value from the previous row.
    parsefile(file, true) do items
        @assert length(items) == 2                # Ensures each row has two columns.

        # Header Check
        if lowercase(items[1]) == "stored_energy" # Skips the header row.
            return
        end

        # Data Conversion
        energy = 1000.0 * parse(Float64, items[1]) # GWh -> MWh conversion
        value = parse(Float64, items[2])           # $/MWh

        # Equation Construction: y = cumulative_value + (x - previous_energy) * value
        intercept = cumulative_value - previous_energy * value  # Calculates the intercept for the linear segment 
        push!(equations, LinearEquation(intercept, value))      # and stores it.

        # Update Accumulated Value
        cumulative_value += (energy - previous_energy ) * value
        previous_energy  = energy
    end
    return [LinearEquation(eqn.intercept - cumulative_value, eqn.coefficient) for eqn in equations]
end
