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
using ColorSchemes
using Colors
using HDF5

function main()
    
    return NamedTuple(Base.@locals)
end

function plot_profiles(folder, s=:e, scheme = ColorSchemes.isoluminant_cm_70_c39_n256; kw...)
    (times, z, n) = loadsim(folder)
    tnorm = times.t ./ maximum(times.t)
    for i in eachindex(times.t)
        df = n[i]
        rgb = get(scheme, tnorm[i])
        color = (red(rgb), green(rgb), blue(rgb))
        plt.plot(df[!, Symbol(s)], df.z; color, kw...)
    end
end


function pcolor_profiles(folder, s=:e; log=true, cmap="gnuplot2", kw...)
    (times, z, n) = loadsim(folder)
    list = map(eachindex(times.t)) do i
        df = n[i]
        return df[!, Symbol(s)]
    end
    v = hcat(list...)
    if log
        kw = (;norm=plt.matplotlib.colors.LogNorm(), kw...)
    end
    plt.pcolormesh(times.t, z, v; cmap, kw...)
end


function plot_slice(folder, zslice, s::Symbol=:e; kw...)
    (times, z, n) = loadsim(folder)

    list = map(eachindex(times.t)) do i
        df = n[i]
        return df[!, Symbol(s)]
    end

    iz = searchsortedfirst(z, zslice)
    n = [item[iz] for item in list]
    plt.plot(times.t, n; label=string(s), kw...)
end

plot_slice(folder, zslice, s::String; kw...) = plot_slice(folder, zslice, Symbol(s); kw...)
plot_slice(folder, zslice, v::AbstractVector; kw...) = foreach(s -> plot_slice(folder, zslice, s; kw...), v)

function loadsim(folder)
    try
        return loadsim(folder, Val(:hdf5))
    catch e
        return loadsim(folder, Val(:csv))
    end    
end

function loadsim(folder, ::Val{:csv})
    times = CSV.read(joinpath(folder, "times.csv"), DataFrame)
    z = nothing
    n = map(eachindex(times.t)) do i 
        # The conversion to Int is because CSV sometimes returns ChainedVectors for long columns.
        # The indices to the ChainedVector are not integeres but can be converted into them.
        fname = joinpath(folder, @sprintf("n-%04d.csv", convert(Int, i)))
        df = CSV.read(fname, DataFrame)
        z = isnothing(z) ? df.z : z
        return df
    end
    
    return (times, z, n)
end


function loadsim(folder, ::Val{:hdf5})
    local times1, z, species, n1
    h5open(joinpath(folder, "output.hdf5")) do fp
        times1 = Array(fp["times"])
        species = Array(fp["species"])
        z = Array(fp["z"])
        n1 = Array(fp["n"])
    end
    times = DataFrame(t=times1)
    n = map(eachindex(times.t)) do i
        DataFrame("z" => z, (s => n1[j, :, i] for (j, s) in enumerate(species))...)        
    end

    return (times, z, n)
end



end

if abspath(PROGRAM_FILE) == @__FILE__
    PlotEventStorm.main()
end
