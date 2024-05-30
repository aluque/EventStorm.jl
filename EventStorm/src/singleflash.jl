#=
Simulate a single flash.
=#

"""
Call-back for the slow-system solver, whenever a flash is encountered.
"""
function flash!(integrator)
    (;conf, ws, flash_integrator) = integrator.p
    (;z, ngas, frs, rs, Ipeak_median, Ipeak_log_std, Ipeak_cutoff, storm_distance, storm_extension) = conf
    n = integrator.u
    
    next!(PROGRESS_REF[]; showvalues = [(:t, string(integrator.t))])
        
    # Sample from distributions
    rho = sqrt((storm_extension * randn() + storm_distance)^2 + (storm_extension * randn())^2)
    
    Ipeak = exp(log(Ipeak_median) + randn() * Ipeak_log_std)
    if isfinite(Ipeak_cutoff)
        while Ipeak > Ipeak_cutoff
            Ipeak = exp(log(Ipeak_median) + randn() * Ipeak_log_std)
        end
    end
    
    push!(IPEAK_SAMPLES, Ipeak)
    push!(RHO_SAMPLES, rho)
    
    singleflash_run!(flash_integrator, n, rho, Ipeak, conf, ws)
    n1 = flash_integrator.u.x[1]

    @batch for k in axes(n, 2)
        mapspecies!(@view(n[:, k]), @view(n1[:, k]), rs, frs)
    end

    u_modified!(integrator, true)
end


function singleflash_setup(conf, ws)
    # rho and Ipeak set on a per-event basis
    p = (;rho=0.0, Ipeak=1.0, conf, ws)
    (;z, tl, n1, save_flash) = conf
    
    (;frs, krange) = conf

    # These are just place-holder arrays to set the problem with the appropriate types
    n0 = zeros(nspecies(frs), length(krange))
    M0 = zeros(length(krange) + 2)
    u0 = ArrayPartition(n0, M0)
    
    # DiffEqs.jl call the derivs functions on init, so we have to set propagators to something
    # reasonable to prevent exceptions.  The values are later ignored anyway.
    for k in eachindex(krange)
        k1 = krange[k]
        r = SA[0.0, 0.0, z[k1]]
        (;props) = ws[k]
        propagators!(props, tl, r)
    end

    prob = ODEProblem{true}(self_attenuation_derivs!, u0, (0, 1e-3), p)
    integrator = init(prob, Tsit5(), reltol=1e-4,
                      dense=false, save_everystep=save_flash, save_start=save_flash,
                      save_end=save_flash,
                      initialize_save=save_flash)

    return (;integrator, prob)
end


function singleflash_run!(integrator, n, rho, Ipeak, conf, ws; extra_time=1e-3)
    (;z, r1, r2, source_duration, ngas, frs, rs, n1, krange, tl) = conf
    p = (;rho, Ipeak, conf, ws)

    (t1, _) = tminmax(SA[rho, 0.0, z[krange[begin]]], r1, r2)
    (_, t2) = tminmax(SA[rho, 0.0, z[krange[end]]], r1, r2)

    # We allow for an extra time even if the field vanishes we still may have some fast
    # reactions going on for a while.  Perhaps 1 ms is too short; we have to experiment with this.
    t2 += source_duration + extra_time

    n1 = integrator.u.x[1]

    electrons = idx(rs, :e)
    for k in axes(n, 2)
        k1 = krange[k]

        # Copy species from the slow to the fast system. Species not present in the slow system
        # are initialized with zero density.
        n1[:, k] .= 0
        mapspecies!(@view(n1[:, k]), @view(n[:, k]), frs, rs)

        r = SA[rho, 0.0, z[k1]]
        (;props) = ws[k]
        propagators!(props, tl, r)
    end

    M = integrator.u.x[2]
    M .= 0
    u0 = ArrayPartition(n1, M)
    
    reinit!(integrator, u0, t0=t1, tf=t2)
    integrator.p = p

    solve!(integrator)
end


"""
Compute the derivatives of the variables in the self-attenuation scheme: M and n.
"""
function self_attenuation_derivs!(du, u, p, t)
    (;rho, Ipeak, conf, ws) = p
    (;z, ngas, frs, rs, n1, krange, tl) = conf

    mun = 8.985943e23
    
    # Expand ArrayPartition
    (n, M) = u.x
    (dn, dM) = du.x
    
    dz = z[2] - z[1]
    dM .= 0
    
    electrons = idx(frs, :e)
    @batch for k in eachindex(krange)
        k1 = krange[k]

        # Not compatible with batch:
        # (;props, c) = ws[k]
        props = ws[k].props
        
        r = SA[rho, 0.0, z[k1]]

        # Think about this: we have to assume a planar wave so we have to select the 'best' incidence
        # angle. Here we average the incidence cosine.        
        costheta = sum(tli -> cosinc(pos(tli), r), tl) ./ length(tl)
        Eq, Ei, Er = free_electric_field(fieldcomponents, r, t, Ipeak, tl, props)

        Er *= exp(-M[k + 1])

        en = norm(Er) / ngas[k1] / co.Td
        sigma = co.elementary_charge * n[electrons, k] * mun / ngas[k1]

        dM[k + 1] = sigma / co.epsilon_0 / 2 - co.c * costheta * (M[k + 1] - M[k]) / dz
        dn[:, k] .= Chemise.derivs(@view(n[:, k]), frs, en; x=k1)
    end    
end


"""
Compute the electric field at time `t` given the status in `y` with parameters `p`.  Store the result
in `ef`. If `reduced==true` computes the reduced electric field in Td.

This is only used for debugging/plotting; normally the field is computed internally and discarded
each time step.
"""
function compute_field!(ef, u, p, t; reduced=false, free=false)
    (;rho, Ipeak, conf, ws) = p
    (;z, ngas, krange, tl) = conf

    (n, M) = u.x
    
    @batch for k in eachindex(krange)
        k1 = krange[k]
        
        # Not compatible with batch:
        # (;props, c) = ws[k]
        props = ws[k].props
        r = SA[rho, 0.0, z[k1]]
        
        # Think about this: we have to assume a planar wave so we have to select the 'best' incidence
        # angle. Here we use the angle with the first dipole
        costheta = sum(tli -> cosinc(pos(tli), r), tl) / length(tl)

        Eq, Ei, Er = free_electric_field(fieldcomponents, r, t, Ipeak, tl, props)
        if !free
            Er *= exp(-M[k + 1])
        end

        ef[k] = Er

        if reduced
            ef[k] = ef[k] / ngas[k1] / co.Td
        end
    end
end    
