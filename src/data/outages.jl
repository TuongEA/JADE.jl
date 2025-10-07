#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.


# Tuong: This function is obsolete, it is used to read the old format files "generator_outages.csv" and "line_outages.csv"
# These files are now replaced by "station_outages.csv" and "transmission_outages.csv"
function getoutages(outage_file::String, years::Vector{Int}, weeks::Vector{Int}, loadblocks::Vector{Symbol},
    demand::TimeSeries{Dict{Tuple{Symbol,Symbol},Float64}}, T::TimeSeries{Dict{Symbol,Float64}},
)
    
    # Initialize an empty TimeSeries with the same time points as demand
    empty_data = [Dict{Tuple{Symbol, Symbol}, Float64}() for _ in keys(demand)]
    outages = TimeSeries(demand.startpoint, empty_data)
    
    # Parse the outage file line by line
    parsefile(outage_file, true) do items
        @assert length(items) == 5

        # Skip header row
        if lowercase(items[1]) == "year"
            return
        end

        years_ = process_set_items(items[1], years)
        weeks_ = process_set_items(items[2], weeks)
        loadblocks_ = process_set_items(items[3], loadblocks)
        name = str2sym(items[4])
        outage = parse(Float64, items[5])

        for blk in loadblocks_, y in years_, w in weeks_
            tp = TimePoint(y, w)
            key = (name, blk)

            if haskey(outages[tp], key)
                error("Duplicate outage entry for $(key) in week $(tp)")
            end
            outages[tp][key] = outage
                
        end
    end

    return outages
end
