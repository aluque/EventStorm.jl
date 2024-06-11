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
using RecursiveArrayTools
using Printf
using Distributions
using Dates
using Random
using ProgressMeter
using HDF5

using DocStringExtensions

@template DEFAULT =
    """
    $(SIGNATURES)
    $(DOCSTRING)
    """

using Constants: co
using DipoleRadiators
using DipoleRadiators: FieldComponents, Propagator, image, pos, remotefield
using Chemise

const DATA_DIR = normpath(joinpath(@__DIR__, "..", "data"))
const IPEAK_SAMPLES = Vector{Float64}()
const RHO_SAMPLES = Vector{Float64}()
const PROGRESS_REF = Ref{Progress}()

include("softstep.jl")
include("electrons.jl")
include("load_data.jl")
include("lxcat.jl")
include("chemistry.jl")
include("singleflash.jl")
include("scratch.jl")
include("inputs.jl")

function _main(;
               # Place-holders filled by wrap_input in inputs.jl
               _input="",
               _date=nothing,
               _git_commit=nothing,
               _git_dirty=nothing,
               
               T = 200.0,
               Te = T,
               
               h=0 * co.kilo,
               l=8 * co.kilo,
               
               rise=2 * co.micro,
               decay=10 * co.micro,
               plateau=0.0,
               pulse_speed=2e8,
               
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
               Ipeak_cutoff = Inf,

               # boundaries for integration
               hmin = 70 * co.kilo,
               hmax = 95 * co.kilo,
               
               # The storm is simulated as 2d gaussian in space. This is the standard deviation
               storm_extension = 100 * co.kilo,

               # The time distribution is either a gaussian ("gauss") or uniform ("uniform").
               storm_time_dist = "uniform",
               storm_duration = 1.5 * co.hour,
               storm_center_time = 4 * co.hour,
                              
               # Expected number of flashes in the storm (actual number may differ).
               storm_expected_flashes = 1000.0,

               # Distance of center of the storm to the observer.
               storm_distance = 100 * co.kilo,
               
               # Final time of the simulation
               final_time = storm_center_time + 5 * storm_duration,

               # Starting time of the simulation
               start_time = 0.0,
               
               # Air composition
               comp=Dict("N2" => 0.8, "O2" => 0.2),
               
               # Simulation name.
               name = splitext(splitdir(_input)[2])[1],
               
               outfolder = joinpath(splitdir(_input)[1], name,
                                    String(rand('A':'Z', 3)) * "-" * Dates.format(now(), "yyyymmdd-HHMMss")),

               # Save as CSV files
               csv_output= false,

               # Save as hdf5
               hdf5_output = true,
               
               random_seed = rand(UInt),
               
               # Time between outputs
               output_dt = 300,
               
               # resolution
               points_per_km = 2,
               
               # If not nothing, use a different gas density file
               gas_density_fname = nothing,
               
               # If not nothing use an iri profile for electron density
               iri_electron_density_fname = nothing,
               
               # Time for the fast reactions after the end of the source
               extra_time = 50e-3,
               
               # Whether to save flash data to inspect / plot / debug
               save_flash = false,
               
               #Perform a baseline run that ignores the fast reactions
               baseline_run = false,

               # If false, returns the configuration but does not run the simulation
               run = true,               
               )


    Polyester.reset_threads!()

    r1 = @SVector([0, 0, h])
    r2 = r1 + @SVector([0, 0, l])
    
    # Normalized source with an Ipeak=1A
    source = SoftStep(r2, r1, 1.0, rise, decay, plateau)

    # A time duration that includes the full source
    source_duration = rise + decay + plateau
    
    tl = trans_line(source, pulse_speed, ndipoles)

    # No absorption below this line
    hcut = 50 * co.kilo
    z = LinRange(hcut, hmax, 1 + ceil(Int, (hmax - hcut) * points_per_km / co.kilo))

    ##
    ## READ DENSITY PROFILES
    ##
    gas_profiles = load_gas_profiles(joinpath(DATA_DIR, "US_Std_1976.dat"))
    cr_profile = load_cr_profile(joinpath(DATA_DIR, "Thomas1974_Production.dat"))    
    waccm_profiles = load_waccm(joinpath(DATA_DIR, "waccm_fg_l38.dat"))
    nrlmsis = load_nrlmsis(joinpath(DATA_DIR, "nrlmsis2.dat"))

    if isnothing(iri_electron_density_fname)
        electron_density = LogInterpolatedElectronDensity(joinpath(DATA_DIR, "earth", "electrons.dat"))
    else
        electron_density = load_iri(iri_electron_density_fname)
    end
    
    ne = electron_density.(z)
    q = cr_profile.(z)
    if isnothing(gas_density_fname)
        ngas = interp_gas_profile(waccm_profiles, :nair, z)
    else
        f_ngas = load_gas_density(gas_density_fname)
        ngas = f_ngas.(z)
    end

    n_o3 = interp_gas_profile(gas_profiles, :O3, z)
    n_no = interp_gas_profile(gas_profiles, :NO, z)

    # Following Brasseur and Salomon, CO2 and H2O seem to be important for negative ion chemistry
    n_co2 = interp_gas_profile(gas_profiles, :CO2, z)
    n_h2o = interp_gas_profile(gas_profiles, :H2O, z)

    n_o = nrlmsis.o.(z)
    
    ##
    ## REACTION SETS (see chemistry.jl)
    ##    
    # O2 absorption of ionizing Lyman α and β. Flux and cross sections from Kotovsky 2016b
    R = 1e6 * co.centi^-2
    # At 1025.7 A
    lyman_β = @. 30 * R * exp(-nrlmsis.cum_o2(z) * 1.55e-18 * co.centi^2)

    # Kotovsky and Moore give 1.04e-22 cm^2 but that seems to be a typo, as Watanabe (their ref) gives
    # 1.04e-20 cm^2 citing Preston.
    # At 1215.6 A
    lyman_α = @. 5e3 * R * exp(-nrlmsis.cum_o2(z) * 1.04e-20 * co.centi^2)

    fixed_dens = [:N2 => comp["N2"] .* ngas,
                  :O2 => comp["O2"] .* ngas,
                  :M => ngas,
                  :O => n_o,
                  :O3 => n_o3,
                  :NO => n_no,
                  :CO2 => n_co2,
                  :H2O => n_h2o,
                  :Q => q,
                  :β => lyman_β,
                  :α => lyman_α]

    rs = slow_reactions(T, fixed_dens)
    @info "Number of SLOW reactions: $(length(rs.reactions))"
    
    frs = fast_reactions(T, comp; ngas)
    @info "Number of FAST reactions: $(length(frs.reactions))"

    ##
    ## CONFIG AND WORKSPACE
    ##
    kmin = searchsortedfirst(z, hmin)
    kmax = searchsortedlast(z, hmax)
    krange = kmin:kmax
    
    n1 = zeros(nspecies(frs), length(krange))
    
    conf = Config(;z, ngas, krange, r1, r2, source_duration,
                  Ipeak_median, Ipeak_log_std, Ipeak_cutoff,
                  tl, storm_distance, storm_extension, n1, rs, frs, save_flash, baseline_run,
                  extra_time)
    
    ws = [Workspace(Float64, length(tl)) for _ in krange]
    

    ##
    ## Sample flash times
    ##
    Random.seed!(random_seed)
    nevents = convert(Int, rand(Poisson(storm_expected_flashes)))
    @info "Number of flashes to simulate:" nevents

    if storm_time_dist == "gauss"
        event_times = storm_center_time .+ storm_duration .* randn(nevents)
    elseif storm_time_dist == "uniform"
        event_times = storm_center_time .+ storm_duration .* (rand(nevents) .- 0.5)
    else
        @error "Unknown storm_time_dist" storm_time_dist
    end
    
        
    sort!(event_times)
    PROGRESS_REF[] = Progress(nevents; showspeed=true)
    if length(event_times) > 0 && final_time < last(event_times)
        final_time = last(event_times) + 60.0
        @warn "final_time was adapted to include all flashes" final_time
    end

    start_time = 0.0
    if length(event_times) > 0 && start_time > first(event_times)
        start_time = first(event_times)
        @warn "start_time was adapted to include all flashes" start_time
    end
    
    ##
    ## Set up the flash sub-solver
    ##
    flash_integrator, flash_prob = singleflash_setup(conf, ws)


    ##
    ## Set up the long-time solver
    ##
    n0 = zeros(nspecies(rs), length(krange))
    for k in axes(n0, 2)
        k1 = krange[k]
        n0[:, k] = Chemise.init(rs, Dict(:e => ne[k1], :O3 => 1e9), 0.0)
    end

    cb = DiffEqCallbacks.PresetTimeCallback(event_times, flash!, save_positions=(false, false))
    
    tspan = (start_time, final_time)    
    prob = ODEProblem(derivs!, n0, tspan, (;conf, ws, flash_integrator))
    
    if !run
        return NamedTuple(Base.@locals)
    end

    ## 
    ## Make sure that output folder exists
    ##
    if !isdir(outfolder)
        mkpath(outfolder)
        @info "$(outfolder) created"
    else
        @warn "$(outfolder) already exists and output may overwrite exisiting files."
    end

    writejl(joinpath(outfolder, name * ".jl"), _main, Base.@locals)
    
    ##
    ## Solve the ODEs
    ##
    sol = solve(prob, Rodas4P(), saveat=output_dt, callback=cb)


    ##
    ## Write output
    ##
    if csv_output
        for (i, n) in enumerate(sol.u)
            fname = joinpath(outfolder, @sprintf("n-%04d.csv", i))
            
            CSV.write(fname, DataFrame(Dict(:z => z[krange],
                                            map(s -> s => n[idx(rs, s), :], species(rs))...)))
        end
        CSV.write(joinpath(outfolder, "times.csv"), DataFrame(t=sol.t))
        CSV.write(joinpath(outfolder, "flashes.csv"),
                  DataFrame(t=event_times, ipeak=IPEAK_SAMPLES, rho=RHO_SAMPLES))
        
    end

    if hdf5_output
        specs = species(rs)
        fname = joinpath(outfolder, "output.hdf5")
        h5open(fname, "w") do fp
            fp["times"] = sol.t
            fp["species"] = collect(map(string, specs))
            fp["z"] = collect(z[krange])
            fp["n"] = cat(sol.u..., dims=3)

            g = create_group(fp, "flashes")
            g["times"] = event_times
            g["ipeak"] = IPEAK_SAMPLES
            g["rho"] = RHO_SAMPLES
        end        
    end
    
    @info "Done"
    return NamedTuple(Base.@locals)
end

@kwdef mutable struct Config{Z<:AbstractVector, T, Dip, RS, FRS}
    "Altitudes (m)"
    z::Z
    
    "Gas density"
    ngas::Vector{T}

    "Range kmin:kmax where electron density is allowed to change"
    krange::UnitRange{Int}
    
    "Distance observer-center of storm"
    storm_distance::T
    
    "Standard deviation of the storm area"
    storm_extension::T

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

    "Cutoff of peak currents."
    Ipeak_cutoff::T

    "Transmission line dipoles"
    tl::Vector{Dip}
        
    "Space for densities in the fast reaction set"
    n1::Array{T, 2}
    
    "Reaction set"
    rs::RS

    "Fast reaction set"
    frs::FRS
    
    "Time in the fast system after the end of the source"
    extra_time = 50e-3
    
    "Number of gridpoints when integrating in time"
    time_gridpoints::Int = 256

    "Save ODE solution (useful only for debugging)"
    save_flash::Bool = false

    "Perform a baseline run that ignores the fast reactions"
    baseline_run::Bool = true
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

function derivs!(dn, n, p, t)
    (;conf, ws) = p
    (;rs, krange) = conf
    
    @batch for i in axes(n, 2)
        j = krange[i]
        dn[:, i] = derivs(@view(n[:, i]), rs; x=j)
    end
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
Compute the electric field at `z[i]`, `ρ` produced by a peak current `Ipeak` when the log-attenuation
is `latt`. `props` contains the propagators and `c` the attenuation factor for each dipole. 
"""
function attenuated_electric_field(r, t, Ipeak, tl, latt, props, c)    
    attenuation!(c, latt, tl, r)
    field = zero(SVector{3, typeof(t)})
    
    for j in eachindex(tl)
        field += total(remotefield(tl[j], props[j], t) * c[j])
    end

    return field * Ipeak
end


"""
Compute the free-space electric field at `z[i]`, `ρ` produced by a peak current
`Ipeak`. `props` contains the propagators. `f` usually is the function to extract the
field from the output ff `DipoleRadiators`: `total` for the full field, `induction`, `radiation` or
`static` (any other function would work though).
"""
function free_electric_field(f::Function, r, t, Ipeak, tl, props)
    field = sum(j -> f(remotefield(tl[j], props[j], t)), eachindex(tl))

    return field * Ipeak
end

nonradiation(rf) = induction(rf) + static(rf)
fieldcomponents(rf::DipoleRadiators.FieldComponents) = SA[rf.static, rf.induction, rf.radiation]

"""
Total electric field
"""
free_electric_field(r, t, Ipeak, tl, props) = free_electric_field(total, r, t, Ipeak, tl, props)


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
function trans_line(source, pulse_speed, n=100)
    # Min time considered
    tmin = 0.0
    
    # Max time considered
    tmax = 1e-3
    
    # Time interval for the discretization
    dt = 1e-9
    
    pulse = CurrentPulse(t -> current(source, t), tmin, tmax, dt)

    v = pulse_speed

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

