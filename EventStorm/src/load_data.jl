#= 
Loading data for atmospheric profiles.
=#

"""
Load gas density data
"""
function load_gas_density(T, fname)
    f = CSV.File(fname, header=[:z, :ngas], types=T)
    z1::Vector{T} = f.z
    z = approxrange(z1)
    
    ngas::Vector{T} = f.ngas
    #interp = LinearInterpolation(z, log.(ngas))
    interp1 = interpolate(log.(ngas), BSpline(Cubic(Line(OnGrid()))))
    interp = Interpolations.scale(interp1, z)
    
    wrap = let interp=interp
        z -> exp(interp(z))
    end
    
    return wrap
end

load_gas_density(fname) = load_gas_density(Float64, fname)


"""
Load background gas profiles
"""
function load_gas_profiles(fname)
    stdatm = DataFrame(CSV.File(fname, header=3, skipto=5,
                                delim=' ', ignorerepeated=true))
    rename!(stdatm, :Z => :z)

    c = Dict("z" => co.kilo,
             "T" => 1.0,
             "P" => co.milli * co.bar)

    for name in names(stdatm)
        if name in keys(c)
            stdatm[!, name] .*= c[name]
        else
            stdatm[!, name] .*= @. co.atm * 1e-5 / (co.k * stdatm.T)
        end
    end
    stdatm[!, :nair] = @. stdatm.P / (co.k * stdatm.T)    

    return stdatm
end

function load_waccm(fname)
    fields = "z p T nair O O3 N NO NO2 NO3 N2O N2O5 H OH H2 HO2 H2O2 HNO3 CO"
    f = CSV.read(fname, DataFrame;
                 delim=" ", ignorerepeated=true,
                 skipto=10,
                 header=String.(split(fields)))
    reverse!(f)

    c = Dict("z" => co.kilo,
             "T" => 1.0,
             "p" => co.atm)

    for name in names(f)
        if name in keys(c)
            f[!, name] .*= c[name]
        else
            f[!, name] .*= co.centi^-3
        end
    end

    return f
end


function interp_gas_profile(T, f, col, zto)
    z1::Vector{T} = f.z
    y::Vector{T} = f[!, col]
    
    interp = LinearInterpolation(z1, log.(y))
    # interp1 = interpolate(log.(y), BSpline(Cubic(Line(OnGrid()))))
    # interp = scale(interp1, z)

    return @.(exp(interp(zto)))
end

interp_gas_profile(f, col, zto) = interp_gas_profile(Float64, f, col, zto)


"""
Load background ionization profiles (mostly CR). 
"""
function load_cr_profile(T, fname)
    f = DataFrame(CSV.File(fname, header=[:z, :q], delim=' ', ignorerepeated=true))

    z1::Vector{T} = f.z * co.kilo
    z = approxrange(z1)
    
    q::Vector{T} = f.q * co.centi^-3
    #interp = LinearInterpolation(z, log.(ngas))
    interp1 = interpolate(log.(q), BSpline(Cubic(Line(OnGrid()))))
    interp = Interpolations.scale(interp1, z)
    
    wrap = let interp=interp
        z -> exp(interp(z))
    end
    
    return wrap
end

load_cr_profile(fname) = load_cr_profile(Float64, fname)


"""
Load reaction rate coefficient data for effective ionization.
"""
function load_effective_ionization(T::Type; comp)
    (k025, k042, k027) = [CSV.File(joinpath(DATA_DIR, "swarm", "$(k).dat"), delim=" ", comment="#",
                                   header=[:en, :k], types=T)
                          for k in ["k025", "k042", "k027"]]
    keff::Vector{T} = @.(comp["N2"] * k025.k + comp["O2"] * k042.k - comp["O2"] * k027.k)
    en::Vector{T} = k025.en

    # Some care is required here: I allowing the reduced to go beyond 1000 Td. This is bc higher
    # quadrature orders in Ipeak involve very intense flashes that otherwise would crash the computation
    # however they typically have very small weights assigned.
    interp = linear_interpolation([0.0; en], [0.0; keff], extrapolation_bc=Flat())

    return interp
end

load_effective_ionization(;kw...) = load_effective_ionization(Float64; kw...)


"""
Load NRLMSIS profiles.
"""
function load_nrlmsis(T::Type, fname::String)
    df = CSV.read(fname, delim=' ', ignorerepeated=true, DataFrame)
    
    z::Vector{T} = df[!, "Heit(km)"] .* co.kilo
    n_o2::Vector{T} = df[!, "O2den(cm-3)"] .* co.centi^-3    
    n_cum_o2 = cumlogint(z, n_o2)
    n_cum_o2 .= n_cum_o2[end] .- n_cum_o2
    
    o2 = linear_interpolation(z, n_o2)
    cum_o2 = linear_interpolation(z, n_cum_o2)
    
    return (;o2, cum_o2)
end

load_nrlmsis(fname::String) = load_nrlmsis(Float64, fname)

"""
Cumulative integral of `y` evaluated at `x` assuming that between gridpoints `log(y)` is linear.
""" 
function cumlogint(x, y)
    @assert axes(x) == axes(y)
    
    c = similar(y)
    c[begin] = 0
    
    for i in eachindex(y)[begin: end - 1]
        c[i + 1] = c[i] + (y[i] - y[i + 1]) * (x[i + 1] - x[i]) / log(y[i] / y[i + 1])        
    end

    return c
end

"""    
Guess a range from a vector `v`.  Raises a warning if the values are too far 
from uniform.
"""
function approxrange(v::AbstractVector; atol::Real=0, rtol::Real=0.001)
    @assert issorted(v)
    
    h = diff(v)
    hmean = sum(h) / length(h)
    
    if !all(x -> isapprox(x, hmean; atol, rtol), h)
        @warn "The supplied vector is not approximately uniformly spaced."
        return v
    end

    l = length(v)
    
    return LinRange(v[begin], v[end], l)
end

