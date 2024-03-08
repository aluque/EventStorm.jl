#=
This file is for code useful in debugging or old code that is no longer in use for the main
algorithm.  Note that parts of this may not work any longer.
=#

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

    
    prob = ODEProblem{false}(fast_derivs, u1start, (t1, t2),
                             (conf, k1, r, Ipeak, ngas[k1], latt[k1], ws[k]))
    
    # This is type-unstable; don't know if it can be solved.
    integrator = init(prob, Tsit5(), reltol=1e-6)
    solve!(integrator)
    
    return integrator.u::typeof(u1start)
    #dlogne1(r, t1, t2, Ipeak, ngas[k1], latt[k1], conf, ws[k])
end


function fast_derivs(u, p, t)
    (conf, k1, r, Ipeak, ngas1, latt1, ws1) = p
    (;tl, frs) = conf
    (;props, c) = ws1

    en = norm(electric_field(r, t, Ipeak, tl, latt1, props, c)) / ngas1 / co.Td

    local dudt
    try
        dudt = derivs(u, frs, en, x=k1)
    catch e
        @error "error computing derivatives (field too high?)"
        @show en u
    end
    return dudt
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


