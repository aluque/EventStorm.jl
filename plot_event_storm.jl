# plot_event_storm.jl : Wed Dec 20 15:55:18 2023
"""
# PlotEventStorm
Plot utils for event_storm.jl

## Running the code

```julia
julia> includet("plot_event_storm.jl")     # Needs Revise.jl
julia> PlotEventStorm.main()
```
"""
module PlotEventStorm

using CSV
using DataFrames
using PyPlot
using Printf

function main()
    
    return NamedTuple(Base.@locals)
end

function plot_profiles(folder, s=:e; kw...)
    times = CSV.read(joinpath(folder, "times.csv"), DataFrame)
    for i in eachindex(times.t)
        fname = joinpath(folder, @sprintf("n-%04d.csv", i))
        df = CSV.read(fname, DataFrame)
        plt.plot(df[!, Symbol(s)], df.z; kw...)
    end
end


function pcolor_profiles(folder, s=:e; log=true, cmap="gnuplot2", kw...)
    times = CSV.read(joinpath(folder, "times.csv"), DataFrame)
    z = nothing
    list = map(eachindex(times.t)) do i
        fname = joinpath(folder, @sprintf("n-%04d.csv", i))
        df = CSV.read(fname, DataFrame)
        z = isnothing(z) ? df.z : z
        return df[!, Symbol(s)]
    end
    v = hcat(list...)
    if log
        kw = (;norm=plt.matplotlib.colors.LogNorm(), kw...)
    end
    plt.pcolormesh(times.t, z, v; cmap, kw...)
end

function plot_slice(folder, zslice, s::Symbol=:e; kw...)
    times = CSV.read(joinpath(folder, "times.csv"), DataFrame)
    z = nothing
    list = map(eachindex(times.t)) do i
        fname = joinpath(folder, @sprintf("n-%04d.csv", i))
        df = CSV.read(fname, DataFrame)
        z = isnothing(z) ? df.z : z
        return df[!, Symbol(s)]
    end
    iz = searchsortedfirst(z, zslice)
    n = [item[iz] for item in list]
    plt.plot(times.t, n; label=string(s), kw...)
end

plot_slice(folder, zslice, s::String; kw...) = plot_slice(folder, zslice, Symbol(s); kw...)
plot_slice(folder, zslice, v::AbstractVector; kw...) = foreach(s -> plot_slice(folder, zslice, s; kw...), v)


end

if abspath(PROGRAM_FILE) == @__FILE__
    PlotEventStorm.main()
end
