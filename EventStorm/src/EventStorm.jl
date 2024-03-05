"""
# EventStorm
Simulate the effects of a thunderstorm on the ionosphere by solving a PDE including flashes as
discontinuous events.

## Running the code

```julia
julia> using EventStorm
julia> EventStorm.main()
```
"""
module EventStorm

using StaticArrays
using CSV
using DataFrames
using Interpolations
using LinearAlgebra
using Polyester
using DifferentialEquations
using DiffEqCallbacks
using Printf
using Distributions
using Dates

using Constants: co
using DipoleRadiators
using DipoleRadiators: FieldComponents, Propagator, image, pos, remotefield
using Chemise

const DATA_DIR = joinpath(@__DIR__, "..", "data")

include("softstep.jl")
include("electrons.jl")
include("load_data.jl")
include("lxcat.jl")

function main(;
              T = 200,
              Te = T,

              h=0 * co.kilo,
              l=8 * co.kilo,

              rise=2 * co.micro,
              decay=10 * co.micro,
              plateau=0.0,
              
              # number of dipoles
              ndipoles=100,

              # log-normal distribution for the peak currents [Nag 2016], according to Slyunyaev2018
              # Ipeak_median = 13 * co.kilo,
              # Ipeak_log_std = 0.77,

              # log-normal distribution for the peak currents [Berger 1975], according to Slyunyaev2018
              # Ipeak_median = 32 * co.kilo,
              # Ipeak_log_std = 0.56,

              # log-normal distribution for the peak currents from Ingrid's paper (see fit_ingrid.jl)
              Ipeak_median = 20.36 * co.kilo,
              Ipeak_log_std = 1.14,
              
              # Duration of the storm
              storm_duration = 1 * co.hour,

              # Pre and post relaxation times
              pre_relax = 20000,

              post_relax = 10000,
              
              # number of quadrature points
              nquad = 16,
              
              # Minimum radius to consider influence
              ρmin = 70 * co.kilo,

              # Maximum radius to consider effects
              ρmax = 100 * co.kilo,

              # boundaries for integration
              hmin = 70 * co.kilo,
              hmax = 95 * co.kilo,
              
              # Storm rate in flashes / m^2 / s
              storm_rate = 1000 / (3600 * 1e5^2),

              # Air composition
              comp=Dict("N2" => 0.8, "O2" => 0.2),

              # Output folder
              outfolder = expanduser("~/data/storm/" * String(rand('A':'Z', 3)) * "-" * Dates.format(now(), "yyyymmdd-HHMMss")),
              
              # Time between outputs
              output_dt = 300,
              
              # resolution
              points_per_km = 2)
    
    Polyester.reset_threads!()

    r1 = @SVector([0, 0, h])
    r2 = r1 + @SVector([0, 0, l])
    
    # Normalized source with an Ipeak=1A
    source = SoftStep(r2, r1, 1.0, rise, decay, plateau)

    # A time duration that includes the full source
    source_duration = rise + decay + plateau
    
    tl = trans_line(source, ndipoles)

    # No absorption below this line
    hcut = 50 * co.kilo
    z = LinRange(hcut, hmax, ceil(Int, (hmax - hcut) * points_per_km / co.kilo))

    gas_profiles = load_gas_profiles(joinpath(DATA_DIR, "US_Std_1976.dat"))
    waccm_profiles = load_waccm(joinpath(DATA_DIR, "waccm_fg_l38.dat"))
    cr_profile = load_cr_profile(joinpath(DATA_DIR, "Thomas1974_Production.dat"))    
    electron_density = LogInterpolatedElectronDensity(joinpath(DATA_DIR, "earth", "electrons.dat"))
    nrlmsis = load_nrlmsis(joinpath(DATA_DIR, "nrlmsis2.dat"))
    
    keff = load_effective_ionization(;comp)
    
    ne = electron_density.(z)
    q = cr_profile.(z)
    ngas = interp_gas_profile(waccm_profiles, :nair, z)
    n_o3 = interp_gas_profile(gas_profiles, :O3, z)
    n_o = interp_gas_profile(waccm_profiles, :O, z)
    n_no = interp_gas_profile(waccm_profiles, :NO, z)

    # O2 absorption of ionizing Lyman α and β. Flux and cross sections from Kotovsky 2016b
    R = 1e6 * co.centi^-2
    # At 1025.7 A
    lyman_β = @. 30 * R * exp(-nrlmsis.cum_o2(z) * 1.55e-18 * co.centi^2)

    # Kotovsky and Moore give 1.04e-22 cm^2 but that seems to be a typo, as Watanabe (their ref) gives
    # 1.04e-20 cm^2 citing Preston.
    # At 1215.6 A
    lyman_α = @. 5e3 * R * exp(-nrlmsis.cum_o2(z) * 1.04e-20 * co.centi^2)
    
    ne1 = copy(ne)
    
    kmin = searchsortedfirst(z, hmin)
    kmax = searchsortedlast(z, hmax)
    krange = kmin:kmax
    
    latt = zeros(length(z))

    ##
    ## SLOW, FIELD-INDEPENDENT REACTION SET
    ##
    
    # From Kotovsky & Moore 2016
    F = 1.74e-18 + 1.93e-17 * sind(35)^4
    local rs, frs

    ion_ion_attachment = 10 * 6e-8 * co.centi^3 * sqrt(300 / T)
    let compN2 = comp["N2"], compO2 = comp["O2"], ngas=ngas, n_o=n_o, n_o3=n_o3, n_no=n_no, lyman_β=lyman_β,
        lyman_α=lyman_α, q=q
        
        rs = ReactionSet(
            ["e + O3 -> O- + O2" => Biblio(1e-17, "ref"),
             
             "e + O3 -> O2- + O" => Biblio(1e-15, "ref"),
             
             "e + O2 + O2 -> O2- + O2" =>
             Biblio(1.4e-41 * (300 / Te) * exp(-600 / T) * exp(700 * (Te - T) / Te / T), "ref"),
             
             "e + O2 + N2 -> O2- + N2" =>
             Biblio(1.07e-43 * (300 / Te)^2 * exp(-70 / T) * exp(1500 * (Te - T) / Te / T), "ref"),
             
             "O- + O2 -> e + O3" => Biblio(5.0e-21, "ref"),
             
             "O- + O -> e + O2" => Biblio(2.3e-16, "ref"),
             
             "O2- + O -> e + O3" => Biblio(3.3e-16, "ref"),
             
             "O2- + O -> O- + O2" => Biblio(3.310e-16, "ref"),

             "O- + O3 -> O3- + O" => 5.30e-16,
             "O- + O2 + N2 -> O3- + N2" => 1.10e-42 * (300 / T),
             "O- + O2 + O2 -> O3- + O2" => 1.10e-42 * (300 / T),
             "O2- + O3 -> O3- + O2" => 4e-17,
             
             # Lyman-β
             "O2 + β -> e + O2+" => 0.90e-18 * co.centi^2,

             # Lyman-α
             "NO + α -> e + NO+" => 2.02e-18 * co.centi^2,
             
             # Constant rate of e/pos-ion recombination.  Must be improved
             # "e + N2+ -> " => Biblio(3e-13 * 300 / Te + 4e-13 * (300 / Te)^1.5, "Gordillo-Vazquez2016/JGR"),
             # "e + O2+ -> " => Biblio(3e-13 * 300 / Te + 4e-13 * (300 / Te)^1.5, "Gordillo-Vazquez2016/JGR"),
             # "e + NO+ -> " => Biblio(3e-13 * 300 / Te + 4e-13 * (300 / Te)^1.5, "Gordillo-Vazquez2016/JGR"),
             "e + N2+ -> " => Biblio(4.2e-12, "Gordillo-Vazquez2016/JGR"),
             "e + O2+ -> " => Biblio(4.2e-12, "Gordillo-Vazquez2016/JGR"),
             "e + NO+ -> " => Biblio(4.2e-12, "Gordillo-Vazquez2016/JGR"),
             
             # This was for cluster recombination
             # "e + e -> e" => Biblio(1e-11, "ref"),
             
             # Production rate of electrons from cosmic rays. 
             "O2 -> e + O2+" =>  Biblio(F, "ref"),
             "N2 -> e + N2+" =>  Biblio(F, "ref"),

             "O- + N2+ ->" => ion_ion_attachment,
             "O- + O2+ ->" => ion_ion_attachment,
             "O- + NO+ ->" => ion_ion_attachment,
             "O2- + N2+ ->" => ion_ion_attachment,
             "O2- + O2+ ->" => ion_ion_attachment,
             "O2- + NO+ ->" => ion_ion_attachment,
             "O3- + N2+ ->" => ion_ion_attachment,
             "O3- + O2+ ->" => ion_ion_attachment,
             "O3- + NO+ ->" => ion_ion_attachment,
            
             # Phony reaction to keep track of second positive emissions
             "SPS -> SPS" => 1.0
             ];
            
            fix = [
                :N2 => j -> compN2 * ngas[j],
                :O2 => j -> compO2 * ngas[j],
                :O => j -> n_o[j],
                :O3 => j -> n_o3[j],
                :NO => j -> n_no[j],
                :Q => j -> q[j],
                :β => j -> lyman_β[j],
                :α => j -> lyman_α[j],
            ]
        )
    end
    
    ##
    ## FAST, FIELD-DEPENDENT REACTION SET
    ##
    lx = LxCatSwarmData.load(joinpath(DATA_DIR, "swarm/LxCat_Phelps_20230914.txt"))
    tbl = loadtable(eachcol(lx.data), xcol=:en)
    
    let compN2 = comp["N2"], compO2 = comp["O2"], ngas=ngas
        frs = ReactionSet(["e + N2 -> 2 * e + N2+" => RateLookup(tbl, :C25),
                           "e + N2 -> 2 * e + N2+" => RateLookup(tbl, :C26),
                           "e + O2 -> 2 * e + O2+" => RateLookup(tbl, :C43),
                           "e + O2 -> O + O-" => RateLookup(tbl, :C28),
                           "e + O2 + M -> O2- + M" => RateLookup(tbl, :C27),
                           # Instantaneous emission
                           "e + N2 -> e + SPS" => RateLookup(tbl, :C10)
                           ];
                          fix = [:N2 => j -> compN2 * ngas[j],
                                 :O2 => j -> compO2 * ngas[j],
                                 
                                 # The 1e6 comes from Bolsig+'s normalization of 3-body
                                 :M => j -> ngas[j] / 1e6])
    end
        
    u1 = zeros(nspecies(rs), length(krange))
    
    conf = Config(;z, ngas, krange, ρmin, ρmax, r1, r2, source_duration, Ipeak_median, Ipeak_log_std,
                  tl, storm_rate, keff, ne1, latt, rs, frs, u1)
    
    ws = [Workspace(Float64, length(tl)) for _ in krange]
    
    u0 = zeros(nspecies(rs), length(krange))
    for i in axes(u0, 2)
        j = krange[i]
        u0[:, i] = Chemise.init(rs, Dict(:e => ne[j], :O3 => 1e9), 0.0)
    end

    # Events per unit time
    ν = storm_rate * π * (ρmax^2 - ρmin^2)
    nevents = convert(Int, rand(Poisson(ν * storm_duration)))
    event_times = pre_relax .+ storm_duration .* rand(nevents)
    sort!(event_times)
    
    cb = DiffEqCallbacks.PresetTimeCallback(event_times, flash!, save_positions=(false, false))
    
    tspan = (0.0, pre_relax + storm_duration + post_relax)
    prob = ODEProblem(derivs!, u0, tspan, (;conf, ws))
    
    #integrator = init(prob, Rodas4P())
    i = 0
    times = Float64[]

    if !isdir(outfolder)
        mkdir(outfolder)
        @info "$(outfolder) created"
    else
        @warn "$(outfolder) already exists and output may overwrite exisiting files."
    end

    # integrator = init(prob, Rodas4P(), save_everystep=false, callback=cb)
    # flash!(integrator)
    
    # return NamedTuple(Base.@locals)
    sol = solve(prob, Rodas4P(), saveat=output_dt, callback=cb)

    for (i, u) in enumerate(sol.u)
        fname = joinpath(outfolder, @sprintf("n-%04d.csv", i))
        
        CSV.write(fname, DataFrame(Dict(:z => z[krange],
                                        map(s -> s => u[idx(rs, s), :], species(rs))...)))
    end
    CSV.write(joinpath(outfolder, "times.csv"), DataFrame(t=sol.t))
        
    return NamedTuple(Base.@locals)
end

@kwdef mutable struct Config{Z<:AbstractVector, T, Dip, Kinterp, RS, FRS}
    "Altitudes (m)"
    z::Z
    
    "Gas density"
    ngas::Vector{T}

    "Range kmin:kmax where electron density is allowed to change"
    krange::UnitRange{Int}
    
    "Minimum radius to consider influence"
    ρmin::T
    
    "Maximum radius to consider influence"
    ρmax::T

    "Start point of the source discharge"
    r1::SVector{3, T}

    "End point of the source discharge"
    r2::SVector{3, T}

    "Duration of the source current pulse"
    source_duration::T

    "Median of the current peak distribution"
    Ipeak_median::T

    "log-standard deviation of the current peak distribution"
    Ipeak_log_std::T

    "Transmission line dipoles"
    tl::Vector{Dip}
        
    "Number of flashes per unit time and area"
    storm_rate::T
    
    "Interpolation function for the effective ionization rate coefficient"
    keff::Kinterp

    "Space to store a temporary electron density"
    ne1::Vector{T}

    "Space for updating u"
    u1::Array{T, 2}
    
    "Space to store the log-attenuation."
    latt::Vector{T}

    "Reaction set"
    rs::RS

    "Fast reaction set"
    frs::FRS
    
    "Number of gridpoints when integrating in time"
    time_gridpoints::Int = 256
end

"""
Thread-local workspace.
"""
struct Workspace{T, Prop}
    "Space to store propagators for each dipole"
    props::Vector{Prop}

    "Space to store attenuations."
    c::Vector{T}

    function Workspace(T, n)
        return new{T, DipoleRadiators.Propagator{T}}(
            Vector{DipoleRadiators.Propagator{T}}(undef, n),
            Vector{Float64}(undef, n))
    end

end

function derivs!(du, u, p, t)
    (;conf, ws) = p
    (;rs, krange) = conf
    
    for i in axes(u, 2)
        j = krange[i]
        du[:, i] = derivs(@view(u[:, i]), rs; x=j)
    end
end


"""
Sample the next event.
"""
function sample_next(conf)
    (;storm_rate, ρmin, ρmax) = conf

    # Events per unit time
    ν = storm_rate * π * (ρmax^2 - ρmin^2)

    # Next time
    tnext = -log(rand()) / ν

    # Radius
    ρ = sqrt((ρmax^2 - ρmin^2) * rand() + ρmin^2)

    return (tnext, ρ)
end


function flash!(integrator)
    (;conf, ws) = integrator.p
    (;latt, z, ne1, ngas, frs, rs, Ipeak_median, Ipeak_log_std, u1, ρmin, ρmax, krange) = conf
    u = integrator.u
    
    # Sample from distributions
    ρ = sqrt((ρmax^2 - ρmin^2) * rand() + ρmin^2)
    Ipeak = exp(log(Ipeak_median) + randn() * Ipeak_log_std)
    
    electrons = idx(rs, :e)
    for k in axes(u, 2)
        k1 = krange[k]
        ne1[k1] = u[electrons, k]
    end

    logattenuation!(latt, z, ne1, ngas)
    @batch for k in axes(u, 2)        
        u1start = mapspecies(@view(u[:, k]), frs, rs)

        u1end = flash1!(u1start, k, ρ, Ipeak, conf, ws)

        u1[:, k] .= @view(u[:, k])

        mapspecies!(@view(u1[:, k]), u1end, rs, frs)
        
        #u1[electrons, k] *= exp(flash1!(k, ρ, Ipeak, conf, ws))
    end
    set_u!(integrator, u1)
end


"""
Discrete change to the electron density due to a flash.
"""
function flash1!(u1start, k, ρ, Ipeak, conf, ws)
    (;z, tl, r1, r2, source_duration, latt, ngas, krange) = conf

    k1 = krange[k]
    (;props, c) = ws[k]

    r = SA[ρ, 0.0, z[k1]]
    propagators!(props, tl, r)
    attenuation!(c, latt[k1], tl, r)

    (t1, t2) = tminmax(r, r1, r2)
    t2 += source_duration

    
    prob = ODEProblem{false}(fast_derivs, u1start, (t1, t2), (conf, k1, r, Ipeak, ngas[k1], latt[k1], ws[k]))
    
    # This is type-unstable; don't know if it can be solved.
    integrator = init(prob, Tsit5())
    solve!(integrator)
    
    return integrator.u::typeof(u1start)
    #dlogne1(r, t1, t2, Ipeak, ngas[k1], latt[k1], conf, ws[k])
end

function fast_derivs(u, p, t)
    (conf, k1, r, Ipeak, ngas1, latt1, ws1) = p
    (;tl, frs) = conf
    (;props, c) = ws1

    en = norm(electric_field(r, t, Ipeak, tl, latt1, props, c)) / ngas1 / co.Td
                 
    return derivs(u, frs, en, x=k1)
end


"""
Map the densities `n2` of a `ReactionSet` `rs2` into a `SVector` for the `ReactionSet` `rs1`.
Species not in `rs2` are left as `default`.
"""
function mapspecies(n2, rs1, rs2; default=zero(eltype(n2)))
    s1 = species(rs1)
    tpl = map(s1) do spec
        i = idx(rs2, spec)
        isnothing(i) ? default : n2[i]
    end
    return SVector(tpl)
end


"""
Map the densities `n2` of a `ReactionSet` `rs2` into `rs1` overwriting the vector `n1`.
Species not in `rs2` are left unchanged.
"""
function mapspecies!(n1, n2, rs1, rs2)
    s1 = species(rs1)
    for (i, spec) in enumerate(s1)
        j = idx(rs2, spec)
        isnothing(j) || (n1[i] = n2[j])
    end
end


"""
Integrate the effect a single pulse at `r` from a discharge with peak current `Ipeak` lasting
`source_duration` from a transmission line `tl` with pre-computed propagators `props` and 
log-attenuation `lat`.  Uses trapezoidal integration with n nodes.
See `electric_field` for `c`.
"""
function dlogne1(r, t1, t2, Ipeak, gas_dens, latt, conf, ws1)
    (;tl, time_gridpoints, keff) = conf
    (;props, c) = ws1
    
    n = time_gridpoints
    
    @assert t2 > t1

    h = (t2 - t1) / n
    
    function f(t)
        return gas_dens * keff(norm(electric_field(r, t, Ipeak, tl, latt, props, c)) / gas_dens / co.Td)
    end
        
    int = h * (f(t1) + f(t2)) / 2
    for k in 1:n - 1
        tk = h * k + t1
        int += h * f(tk)
    end

    return int
end


"""
Compute the electric field at `z[i]`, `ρ` produced by a peak current `Ipeak` when the log-attenuation
is `latt`. `props` contains the propagators and `c` the attenuation factor for each dipole. 
"""
function electric_field(r, t, Ipeak, tl, latt, props, c)
    
    attenuation!(c, latt, tl, r)
    field = zero(SVector{3, typeof(t)})
    
    for j in eachindex(tl)
        field += total(remotefield(tl[j], props[j], t) * c[j])
    end

    return field * Ipeak
end



"""
Compute the propagators of the electric field from each dipole in the transmission line `tl` to the 
point `r`.
"""
function propagators!(props, tl, r)
    @assert axes(tl) == axes(props)

    for i in eachindex(tl)
        props[i] = Propagator(tl[i], r)
    end

end


"""
Compute the log-attenuation factor at `z` from electron density `ne` and gas density `ngas`.
`mun` is the reduced mobility of electrons. The integral is performed assuming a log-linear
dependence of the conductivity on `z`.
"""
function logattenuation!(latt, z, ne, ngas, mun=8.985943e23)
    T = promote_type(eltype(latt), eltype(z), eltype(ne), eltype(ngas), typeof(mun))
    
    sigma_pre = zero(T)
    z_pre = 2 * z[begin] - z[begin + 1]
    latt_pre = zero(T)
    
    for i in eachindex(z)
        @assert ne[i] > 0
        sigma = convert(T, co.elementary_charge) * ne[i] * mun / ngas[i]
        dlatt = (sigma - sigma_pre) * (z[i] - z_pre) / log(sigma / sigma_pre)
        latt[i] = latt_pre - sqrt(co.mu_0 / co.epsilon_0) * dlatt / 2
        
        (sigma_pre, z_pre, latt_pre) = (sigma, z[i], latt[i])
    end

    return nothing
    # sigma = @. co.elementary_charge * ne * mun / ngas
    # csigma = cumlogint(z, sigma)
    
    #return @. -sqrt(co.mu_0 / co.epsilon_0) * csigma / 2
end


"""
Set the attenuation factor in `f` for all dipoles in the transmission line `tl` as
exp(a / cos(θ)) where `a` is constant θ is the incidence angle for rays propagating 
from each dipole towards `r`.
"""
function attenuation!(f, a, tl, r)
    for i in eachindex(tl)
        f[i] = exp(a / cosinc(pos(tl[i]), r))
    end

    return f
end


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
Cosine of the incidence angle for ray propagating from `r1` to `r2` and inciding into a
 horizontal plane.
"""
function cosinc(r1, r2)
    dr = r2 - r1
    return dr[3] / norm(dr)
end


"""
Construct the transmission-line representation given a current `source` with `n` dipoles.
"""
function trans_line(source, n=100)
    # Min time considered
    tmin = 0.0
    
    # Max time considered
    tmax = 1e-4
    
    # Time interval for the discretization
    dt = 1e-9
    
    pulse = CurrentPulse(t -> current(source, t), tmin, tmax, dt)

    v = 2.0e8 #1.5e8

    # No decay
    λ = 1e6

    tl = mtle(pulse, source.r1, source.r2, v, λ, n; mirror=true, w0=1.0, t0=0.0)
end


""" 
Compute the max and min time retardations at `rf` caused by a source at `r1`-`r2`.  Considers 
the delays of the mirror images.
"""
function tminmax(rf, r1, r2)
    r1m = @SVector [r1[1], r1[2], -r1[3]]
    r2m = @SVector [r2[1], r2[2], -r2[3]]

    minl2 = maxl2 = dot(rf - r1, rf - r1)

    l2 = dot(rf - r2, rf - r2)
    minl2 = min(minl2, l2)
    maxl2 = max(maxl2, l2)

    l2 = dot(rf - r1m, rf - r1m)
    minl2 = min(minl2, l2)
    maxl2 = max(maxl2, l2)

    l2 = dot(rf - r2m, rf - r2m)
    minl2 = min(minl2, l2)
    maxl2 = max(maxl2, l2)
    
    return sqrt(minl2) / co.c, sqrt(maxl2) / co.c
end


end

if abspath(PROGRAM_FILE) == @__FILE__
    EventStorm.main()
end

