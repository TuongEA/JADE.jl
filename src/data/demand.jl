#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

#----------------------------------------------------
# Prepare demand-related data
#----------------------------------------------------
"""
    getdemand(file::String, durations::JADE.TimeSeries{Dict{Symbol,Float64}})
    This function reads a demand data file and constructs a time series of demand values for each (node, block) pair. It scales these values by the corresponding durations provided in a separate time series.

    Inputs:
    file::String: Path to a CSV or text file containing demand data.
        # Example
            % Demand (MW):
            NODE,YEAR,WEEK,peak,shoulder,offpeak
            NI,2005,1,2132.15,1616.51,1116.58
    durations::JADE.TimeSeries{Dict{Symbol,Float64}}: A time series mapping each time point to a dictionary of durations for each demand block.

    Outputs:
    demands: A time series of scaled demand values for each (node, block) pair.
    nodes: A list of all nodes found in the data.

    Behavior:
    Parses the file to extract demand blocks and node-specific demand values.
    Constructs a time series of demand values indexed by timepoints.
    Scales each demand value by its corresponding duration.
    Returns the final demand time series and the list of nodes.
"""
function getdemand(file::String, durations::JADE.TimeSeries{Dict{Symbol,Float64}})
    
    data = Dict{TimePoint,Dict{Tuple{Symbol,Symbol},Float64}}()    # A dictionary mapping TimePoint → Dict{(node, block) => demand value}.
    demand_blocks = Symbol[]                                       # List of demand block names (e.g., :peak, :shoulder, :offpeak).
    nodes = Symbol[]                                               # List of node names (e.g., :NI, :HAY, :SI).
    start_year = 0
    start_week = 0

    # Reads the file line by line. Each line is split into items.
    parsefile(file, true) do items
        @assert length(items) >= 4       # At least one demand block
        if lowercase(items[1]) == "node" # Identifies the header row. 
            for it in items[4:end]       # Extracts demand block names from columns 4 onward.
                if it != ""
                    push!(demand_blocks, str2sym(it))
                end
            end
        else # Data row
            timepoint = TimePoint(parse(Int, items[2]), parse(Int, items[3]))  # Parses year and week to create a TimePoint.
            if length(data) == 0              # First data row will contain start year and start week 
                start_year = timepoint.year
                start_week = timepoint.week
            end
            node = str2sym(items[1])          # Converts node name to a Symbol.
            if !(node in nodes)               # Add node to list of nodes if the node is not yet in the list
                push!(nodes, node)
            end

            # checks whether the dictionary data already contains an entry for the key timepoint. 
            if !haskey(data, timepoint) # If it does not, it creates a new entry with that key 
                data[timepoint] = Dict{Tuple{Symbol,Symbol},Float64}() # and assigns it an empty dictionary of type Dict{Tuple{Symbol,Symbol},Float64}.
            end

            # Stores demand values for each (node, block) pair at the given timepoint.
            for (i, block) in enumerate(demand_blocks) # iterates over each demand block (like :peak, :shoulder, :offpeak) 
                data[timepoint][(node, block)] = parse(Float64, items[i+3]) # and assigns the corresponding demand value from the items array
            end
        end
    end

    # Collects and sorts timepoints.
    tps, vals = collect(keys(data)), collect(values(data))
    p = sortperm(tps)

    # Constructs a TimeSeries object starting from the first timepoint.
    demands = TimeSeries{Dict{Tuple{Symbol,Symbol},Float64}}( TimePoint(start_year, start_week), vals[p],
    )

    # Multiplies each demand value by its corresponding duration.
    # k[2] extracts the block name from the (node, block) tuple.
    for tp in tps
        for k in keys(demands[tp])
            demands[tp][k] *= durations[tp][k[2]]
        end
    end

    return demands, nodes
end

#---------------------------------------------
# Load shedding
#---------------------------------------------

"""
    getdemandresponse(file::String, demand, durations, nodes, years, weeks, loadblocks)
    This function reads a demand response configuration file and constructs time series of demand response tranches for both power-based and energy-based programs.

    Inputs:
    file: Path to the demand response tranche definition file.
    demand: A time series of demand values by node and load block.
    durations: A time series of durations for each load block.
    nodes, years, weeks, loadblocks: Lists of valid values used to validate and filter tranche applicability.

    Outputs:
    ptranches: A time series of power-based demand response tranches.
    etranches: A time series of energy-based demand response tranches.
    sectors: A list of all sectors defined in the file.

    Behavior:
    Parses each row of the input file to define a demand response tranche.
    Supports both absolute and proportional tranches.
    Applies tranches to the appropriate timepoints, nodes, and load blocks.
    Ensures no duplicate tranche names exist for the same timepoint and node.
"""
function getdemandresponse(file::String, demand, durations, nodes, years, weeks, loadblocks)
    ptranches = Dict{Symbol,Dict{Symbol,Dict{Tuple{Symbol,Symbol},Tranche}}}[]
    etranches = Dict{Symbol,Dict{Tuple{Symbol,Vector{Symbol}},Vector{Tranche}}}[]

    for i in 1:length(demand.data)
        push!(ptranches, Dict{Symbol,Dict{Symbol,Dict{Tuple{Symbol,Symbol},Tranche}}}())
        push!(etranches, Dict{Tuple{Symbol,Vector{Symbol}},Vector{Tranche}}())
    end

    ptranches = TimeSeries(demand.startpoint, ptranches)
    etranches = TimeSeries(demand.startpoint, etranches)

    sectors = Symbol[]
    parsefile(file, true) do items
        if length(items) == 10
            if lowercase(items[1]) == "sector"
                return
            end

            sector = str2sym(items[1])

            if !(sector in sectors)
                push!(sectors, sector)
            end

            name = str2sym(items[2])
            nodes_ = process_set_items(items[3], nodes)
            years_ = process_set_items(items[4], years)
            weeks_ = process_set_items(items[5], weeks)
            loadblocks_ = process_set_items(items[6], loadblocks)
            mode = lowercase(items[7])
            type_ = lowercase(items[8])
            bound = parse(Float64, items[9])
            bidprice = parse(Float64, items[10])
        else # Tuong: This is for EA Data processing
            @assert length(items) == 7

            if lowercase(items[1]) == "node"
                return
            end

            sector = str2sym(items[3])

            if !(sector in sectors)
                push!(sectors, sector)
            end

            name = str2sym(items[4])
            nodes_ = [str2sym(items[1])]
            years_ = years
            weeks_ = weeks
            loadblocks_ = loadblocks
            mode = "power"
            type_ = "proportional"
            bound = parse(Float64, items[5]) * parse(Float64, items[6])
            bidprice = parse(Float64, items[7])
        end

        if mode == "power" # Tuong: This is for EA Data processing
            for y in years_
                for w in weeks_
                    t = TimePoint(y, w)
                    for n in nodes_
                        if !haskey(ptranches[t], n)
                            ptranches[t][n] =
                                Dict{Symbol,Dict{Tuple{Symbol,Symbol},Tranche}}()
                        end
                        for l in loadblocks_
                            if !haskey(ptranches[t][n], l)
                                ptranches[t][n][l] = Dict{Tuple{Symbol,Symbol},Tranche}()
                            end
                            if haskey(ptranches[t][n][l], (sector, name))
                                error("Two demand response tranches with same name found.")
                            end
                            if type_ == "absolute"
                                ptranches[t][n][l][(sector, name)] =
                                    Tranche(bound, bidprice)
                            elseif type_ == "proportional"
                                if durations[t][l] == 0
                                    ptranches[t][n][l][(sector, name)] =
                                        Tranche(0.0, bidprice)
                                else
                                    ptranches[t][n][l][(sector, name)] = Tranche(
                                        max(
                                            0.0,
                                            bound * demand[t][(n, l)] / durations[t][l],
                                        ),
                                        bidprice,
                                    )
                                end
                            else
                                error(
                                    "Unrecognised demand tranche type " *
                                    type_ *
                                    " use 'absolute' or 'proportional'",
                                )
                            end
                        end
                    end
                end
            end
        elseif mode == "energy"
            for y in years_
                for w in weeks_
                    t = TimePoint(y, w)
                    for n in nodes_
                        if !haskey(etranches[t], n)
                            etranches[t][n] =
                                Dict{Tuple{Symbol,Vector{Symbol}},Vector{Tranche}}()
                        end
                        if !haskey(etranches[t][n], (sector, loadblocks_))
                            etranches[t][n][(sector, loadblocks_)] = Tranche[]
                        end

                        if type_ == "absolute"
                            push!(
                                etranches[t][n][(sector, loadblocks_)],
                                Tranche(bound, bidprice),
                            )
                        elseif type_ == "proportional"
                            push!(
                                etranches[t][n][(sector, loadblocks_)],
                                Tranche(
                                    max(
                                        0.0,
                                        bound * sum(
                                            demand[t][(n, l)] * durations[t][l] for
                                            l in loadblocks_
                                        ),
                                    ) / 1000,
                                    bidprice,
                                ),
                            )
                        else
                            error(
                                "Unrecognised demand tranche type " *
                                type_ *
                                " use 'absolute' or 'proportional'",
                            )
                        end
                    end
                end
            end
        else
            error("Demand-response mode must be either 'power' or 'energy'")
        end
    end
    return ptranches, etranches, sectors
end
