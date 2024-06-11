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

function plot_diff_profiles(folder, s=:e, scheme = ColorSchemes.isoluminant_cm_70_c39_n256; kw...)
    (times, z, n) = loadsim(folder)
    tnorm = times.t ./ maximum(times.t)
    df1 = nothing
    
    for i in eachindex(times.t)
        df = n[i]
        if isnothing(df1)
            df1 = df
            continue
        end
        for col in names(df)
            df1[!, col] .= df[!, col] .- df1[!, col]
        end

        rgb = get(scheme, tnorm[i])
        color = (red(rgb), green(rgb), blue(rgb))
        plt.plot(df1[!, Symbol(s)], df.z; color, kw...)
        df1 = df
    end
end


function pcolor_profiles(folder, s=:e; log=true, cmap="gnuplot2", baseline=nothing, kw...)
    (times, z, n) = loadsim(folder)
    if !isnothing(baseline)
        (btimes, bz, bn) = loadsim(baseline)
        subtract_baseline!(n, bn)
    end
    
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


function pcolor_diff_profiles(folder, s=:e; log=true, cmap="Reds", baseline=nothing, kw...)
    (times, z, n) = loadsim(folder)
    xprev = nothing

    if !isnothing(baseline)
        (btimes, bz, bn) = loadsim(baseline)
        subtract_baseline!(n, bn)
    end
    
    list = map(eachindex(times.t)) do i
        df = n[i]
        if isnothing(xprev)
            xprev = df[!, Symbol(s)]
            return xprev
        else
            x = df[!, Symbol(s)] .- xprev
            xprev = df[!, Symbol(s)]
            return x
        end

    end
    v = hcat(list[begin+1:end]...)
    if log
        vmin = get(kw, :vmin, nothing)
        vmax = get(kw, :vmax, nothing)
        kw = (;norm=plt.matplotlib.colors.LogNorm(;vmin, vmax))
    end
    t1 = @.(0.5 * (times.t[begin:end-1] + times.t[begin+1:end]))

    plt.pcolormesh(t1, z, v; cmap, kw...)
end


function plot_diff_max(folder, s=:e; baseline=nothing, kw...)
    (times, z, n) = loadsim(folder)
    xprev = nothing

    if !isnothing(baseline)
        (btimes, bz, bn) = loadsim(baseline)
        subtract_baseline!(n, bn)
    end
    
    list = map(eachindex(times.t)) do i
        df = n[i]
        if isnothing(xprev)
            xprev = df[!, Symbol(s)]
            return xprev
        else
            x = df[!, Symbol(s)] .- xprev
            xprev = df[!, Symbol(s)]
            return x
        end

    end
    v = hcat(list[begin+1:end]...)
    t1 = @.(0.5 * (times.t[begin:end-1] + times.t[begin+1:end]))
    imax = dropdims(map(x->x[1], argmax(v, dims=1)), dims=1)
    #plt.plot(t1, z[imax], "o", ms=0.5, c="r")
    vmax = [v[imax[i], i] for i in axes(v, 2)]
    
    plt.plot(t1, vmax, "o", ms=0.5, c="r")
end


function plot_slice(folder, zslice, s::Symbol=:e; baseline=nothing, kw...)
    (times, z, n) = loadsim(folder)
    if !isnothing(baseline)
        (btimes, bz, bn) = loadsim(baseline)
        subtract_baseline!(n, bn)
    end

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


function subtract_baseline!(df, bl)
    @assert length(df) == length(bl)
    for i in eachindex(df)
        for col in names(df[i])
            col == :z && continue
            df[i][!, col] .-= bl[i][!, col]
        end
    end

    return df
end


end

if abspath(PROGRAM_FILE) == @__FILE__
    PlotEventStorm.main()
end
