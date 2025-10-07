#  This file is part of JADE source code.
#  This file defines struct types and function applied to time data
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

const SECONDSPERHOUR = 3600
const WEEKSPERYEAR = 52

struct TimePoint
    year::Int
    week::Int
end

# Function to convert TimePoint to total weeks
total_weeks(tp::TimePoint) = tp.year * WEEKSPERYEAR + tp.week
time_point(nwks::Int) = TimePoint(div(nwks-1, WEEKSPERYEAR), (nwks-1) % WEEKSPERYEAR + 1)

# Functions for TimePoint comparison and equality
Base.:<(a::TimePoint, b::TimePoint) = total_weeks(a) < total_weeks(b)
Base.:>(a::TimePoint, b::TimePoint) = b < a
Base.:<=(a::TimePoint, b::TimePoint) = !(b < a)
Base.:>=(a::TimePoint, b::TimePoint) = !(a < b)
Base.:(==)(a::TimePoint, b::TimePoint) = a.year == b.year && a.week == b.week
Base.isless(a::TimePoint, b::TimePoint) = a < b

# Functions for TimePoint (-) and (+) operation
Base.:-(a::TimePoint, b::TimePoint) = total_weeks(a) - total_weeks(b)
Base.:+(a::TimePoint, b::Int) = time_point(total_weeks(a) + b)
Base.:-(a::TimePoint, b::Int) = a + (-b)


struct TimeSeries{T}
    startpoint::TimePoint
    data::Vector{T}
end

Base.length(ts::TimeSeries) = length(ts.data)
Base.isvalid(ts::TimeSeries, tp::TimePoint) = 0 < tp - ts.startpoint + 1 <= length(ts)
Base.getindex(ts::TimeSeries, tp::TimePoint) = isvalid(ts, tp) ? ts.data[tp - ts.startpoint + 1] : throw(KeyError(tp))
Base.getindex(ts::TimeSeries, t::Int) = ts.data[t]

function Base.iterate(ts::TimeSeries, state = (ts.startpoint, ts.data[1]))
    tp, _ = state
    if tp == nothing return nothing end
    isvalid(ts, tp + 1) ? ((tp, ts[tp]), (tp + 1, ts[tp + 1])) : ((tp, ts[tp]), (nothing,nothing))
end

struct TimePointIterator
    startpoint::TimePoint
    endpoint::TimePoint
end

Base.length(tpi::TimePointIterator) = tpi.endpoint - tpi.startpoint + 1
Base.iterate(tpi::TimePointIterator, tp::TimePoint = tpi.startpoint) = tp > tpi.endpoint ? nothing : (tp, tp + 1)
Base.keys(ts::TimeSeries) = TimePointIterator(ts.startpoint, ts.startpoint + length(ts) - 1)
