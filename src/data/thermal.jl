#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

#-----------------------------------------------------
# Parameters relating to thermal stations
#-----------------------------------------------------
struct ThermalStation
    node::Symbol
    fuel::Symbol
    heatrate::Float64
    capacity::Float64
    omcost::Float64
    commission::TimePoint
    decommission::TimePoint
end

struct FuelStorage
    initial::Float64       # initial storage (TJ)
    capacity::Float64      # maximum storage (TJ)
    maxinjection::Float64  # maximum weekly injection (TJ/week)
    maxwithdrawal::Float64 # maximum weekly withdrawal (TJ/week)
end

"""
    getthermalstations(file::String, nodes::Vector{Symbol})

    # Description

    Read thermal station list from `file`.
    The columns MUST be ordered as shown below.

    # Example File

        % Thermal station data,,,,,,,,,
        % Heat rate in GJ/MWh,,,,,,,,,
        % Capacity in MW,,,,,,,,,
        % O&M cost in \$/MWh,,,,,,,,,
        GENERATOR,NODE,FUEL,HEAT_RATE,CAPACITY,OMCOST,START_YEAR,START_WEEK,END_YEAR,END_WEEK
        Stratford_220KV,NI,gas,7.6,377,0,0,0,0,0
        Huntly_e3p,NI,gas,7.2,403,0,2007,23,0,0
        Huntly_main_g1,NI,coal,10.3,250,0,0,0,0,0
        Huntly_main_g2,NI,coal,10.3,250,0,0,0,0,0
"""
function getthermalstations(filename::String, nodes::Vector{Symbol})
    thermal_stations = Dict{Symbol,ThermalStation}()
    for row in CSV.Rows(filename; missingstring = ["NA", "na", "default"], stripwhitespace = true, comment = "%")
        row = _validate_and_strip_trailing_comment(
            row,  [:GENERATOR,:NODE,:FUEL,:HEAT_RATE,:CAPACITY,:OMCOST,:START_YEAR,:START_WEEK,:END_YEAR,:END_WEEK])
        station = str2sym(row.GENERATOR)
        if haskey(thermal_stations, station)
            error("Thermal Station ($station) already given.")
        end
        node = str2sym(row.NODE)
        if !(node in nodes)
            error("Node $node for generator $station not found")
        end
        thermal_stations[station] = ThermalStation(
            node,
            str2sym(row.FUEL),
            parse(Float64, row.HEAT_RATE),
            parse(Float64, row.CAPACITY),
            parse(Float64, row.OMCOST),
            TimePoint(parse(Int, row.START_YEAR), parse(Int, row.START_WEEK)),
            TimePoint(parse(Int, row.END_YEAR), parse(Int, row.END_WEEK)),
        )
    end
    return thermal_stations
end

"""
    getfuelstorage(file::String, nodes::Vector{Symbol})

    # Description

    Read fuel storage data from `file`.
    The columns MUST be ordered as shown below.

    # Example File

        % Initial storage in (TJ)
        % Capacity in (TJ) - ignore this if we only model the total avaialble fuel.
        % Max Injection in (TJ/week) - ignore this if we only model the total avaialble fuel.
        % Max Withdraw in (TJ/week) - ignore this if we only model the total avaialble fuel.
        TRADER,FUEL,INITIAL,CAPACITY,MAX_INJECTION,MAX_WITHDRAW
        GENE,coal,23500,na,na,na
        GENE,gas,14389,na,na,na
        CTCT,gas,13303,na,na,na
        TODD,gas,8787,na,na,na
"""
function getfuelstorage(filename::String)
    storages = Dict{NTuple{2,Symbol},FuelStorage}()

    for row in CSV.Rows(filename; missingstring = ["NA", "na", "default"], stripwhitespace = true, comment = "%")

        row = _validate_and_strip_trailing_comment(row, [:TRADER,:FUEL,:INITIAL,:CAPACITY,:MAX_INJECTION,:MAX_WITHDRAW])

        key = (str2sym(row.TRADER), str2sym(row.FUEL))
        if haskey(storages, key)
            error("Storage $(key) given twice.")
        end

        initial = parse(Float64, coalesce(row.INITIAL, "0.0"))
        capacity = parse(Float64, coalesce(row.CAPACITY, "Inf"))
        maxinjection = parse(Float64, coalesce(row.MAX_INJECTION, "Inf"))
        maxwithdrawal = parse(Float64, coalesce(row.MAX_WITHDRAW, "Inf"))

        storages[key] = FuelStorage(initial,capacity,maxinjection,maxwithdrawal)
    end

    return storages
end

"""
    get_station_fuelstorage_mapping(file::String, fuel_storages::Vector{NTuple{2,Symbol}})

    # Description

    Read station - fuel storage mapping data from `file`.
    The columns MUST be ordered as shown below.

    # Example File

        % Thermal station - fuel storage mapping (can be integrated into thermal_station.csv data)
        GENERATOR,TRADER,FUEL
        Huntly_main_g1,GENE,coal
        Huntly_main_g2,GENE,coal
        Huntly_main_g4,GENE,coal
        Huntly_e3p,GENE,gas
        Huntly_peaker,GENE,gas
        Stratford_220KV,CTCT,gas
        Stratford_peakers,CTCT,gas
        McKee_peakers,TODD,gas
        Junction_Road,TODD,gas
"""
function get_station_fuelstorage_mapping(filename::String, fuel_storages::Vector{NTuple{2,Symbol}})
    stations_storage = Dict{Symbol,NTuple{2,Symbol}}()

    for row in CSV.Rows(filename; missingstring = ["NA", "na", "default"], stripwhitespace = true, comment = "%")
        row = _validate_and_strip_trailing_comment(row,  [:GENERATOR,:TRADER,:FUEL])

        station = str2sym(row.GENERATOR)
        if haskey(stations_storage, station)
            error("Thermal Station ($station) already given.")
        end

        fuelstorage = (str2sym(row.TRADER), str2sym(row.FUEL))
        if !(fuelstorage in fuel_storages)
            error("Fuel storage $fuelstorage for generator $station not found")
        end

        stations_storage[station] = (str2sym(row.TRADER), str2sym(row.FUEL))
    end
    return stations_storage
end
"""
    getfuelcosts(file::String)

    # Description

    Read costs and carbon content for fuels of thermal plant.

    # Example File

        % Fuel costs (\$/GJ except CO2 in \$/tonne) and carbon content (tonnes CO2/GJ)
        ,,coal,diesel,gas,CO2
        CO2_CONTENT,,0.0912,0,0.0528,
        YEAR,WEEK,,,,
        2008,1,4,33.11,5.57,0
        2008,2,4,33.11,5.57,0
"""
# Tuong comments: We should split this thermal_fuel_costs.csv into two files to make it less complicated. 
function getfuelcosts(filename::String)
    start_time, data = nothing, Dict{Symbol,Float64}[]
    rows = CSV.Rows(filename; missingstring = ["NA", "na", "default"],stripwhitespace = true, comment = "%" )
    row, row_state = iterate(rows)
    carbon_content = Dict(
        str2sym("$k") => parse(Float64, row[k]) for
        k in CSV.getnames(row) if !(k in (:Column1, :Column2, :CO2))
    )
    # Skip YEAR,WEEK,... row
    _, row_state = iterate(rows, row_state)
    while (ret = iterate(rows, row_state)) !== nothing
        row, row_state = ret
        time = TimePoint(parse(Int, row.Column1), parse(Int, row.Column2))
        if isempty(data)
            start_time = time
        elseif time != start_time + length(data)
            error("Weeks in $filename must be contiguous")
        end
        d = Dict{Symbol,Float64}(
            str2sym("$k") => parse(Float64, row[k]) for
            k in CSV.getnames(row) if !(k in (:Column1, :Column2))
        )
        push!(data, d)
    end
    return TimeSeries{Dict{Symbol,Float64}}(start_time, data), carbon_content
end

function getfuelandcarboncosts(filename::String)
    start_time = nothing
    data = Dict{Symbol, Float64}[]
    rows = CSV.Rows(filename; missingstring = ["NA", "na", "default"],stripwhitespace = true, comment = "%" )
    for row in rows
        time = TimePoint(parse(Int, row.YEAR), parse(Int, row.WEEK))
        if isempty(data)
            start_time = time
        elseif time != start_time + length(data)
            error("Weeks in $filename must be contiguous")
        end
        d = Dict{Symbol,Float64}(str2sym("$k") => parse(Float64, row[k]) for k in CSV.getnames(row) if !(k in (:YEAR, :WEEK)))
        push!(data, d)
    end
    return TimeSeries{Dict{Symbol,Float64}}(start_time, data)
end

function getfuelcarbonccontents(filename::String)
    rows = CSV.Rows(filename; missingstring = ["NA", "na", "default"],stripwhitespace = true, comment = "%")
    row = first(rows)
    carbon_content = Dict(str2sym("$k") => parse(Float64, row[k]) for k in CSV.getnames(row))
    return carbon_content
end