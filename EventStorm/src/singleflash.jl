#=
Simulate a single flash.
=#

function singleflash_setup(conf, ws)
    # rho and Ipeak set on a per-event basis
    p = (;rho=0.0, Ipeak=1.0, conf, ws)
    (;z, tl, n1) = conf
    
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
                      dense=false, save_everystep = false, save_start=false, save_end=true,
                      initialize_save=false)

    return (;integrator, prob)
end


function singleflash_run!(n, integrator, rho, Ipeak, conf, ws)
    (;z, ngas, frs, rs, n1, krange, tl) = conf
    p = (;rho, Ipeak, conf, ws)

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

    reinit!(integrator, u0)
    integrator.p = p

    @show @allocated for i in integrator; end

    #solve!(integrator)
    
    # ts = range(0, stop = 1e-3, length=10000 + 1)
    # ef = zeros(length(krange), length(ts))
    # for (i, (u, t)) in enumerate(TimeChoiceIterator(integrator, ts))
    #     #compute_field!(@view(ef[:, i]), u, p, t)
    # end

    return NamedTuple(Base.@locals)

end


"""
Compute the derivatives of the variables in the self-attenuation scheme: M and n.
"""
function self_attenuation_derivs!(du, u, p, t)
    (;rho, Ipeak, conf, ws) = p
    (;z, ngas, frs, rs, n1, ρmin, ρmax, krange, tl) = conf

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
        c = ws[k].c
        
        r = SA[rho, 0.0, z[k1]]

        # Think about this: we have to assume a planar wave so we have to select the 'best' incidence
        # angle. Here we average the incidence cosine.        
        costheta = sum(tli -> cosinc(pos(tli), r), tl) ./ length(tl)
        en = (exp(-M[k + 1] / costheta) *
            norm(electric_field(r, t, Ipeak, tl, 0.0, props, c)) / ngas[k1] / co.Td)

        sigma = co.elementary_charge * n[electrons, k] * mun / ngas[k1]

        dM[k + 1] = sigma / co.epsilon_0 / 2 - co.c * (M[k + 1] - M[k]) / dz
        dn[:, k] .= Chemise.derivs(@view(n[:, k]), frs, en; x=k1)
    end    
end


"""
Compute the electric field at time `t` given the status in `y` with parameters `p`.  Store the result
in `ef`. If `reduced==true` computes the reduced electric field in Td.

This is only used for debugging/plotting; normally the field is computed internally and discarded
each time step.
"""
function compute_field!(ef, u, p, t, reduced=false)
    (;rho, Ipeak, conf, ws) = p
    (;z, ngas, krange, tl) = conf

    (n, M) = u.x
    
    @batch for k in eachindex(krange)
        k1 = krange[k]
        
        # Not compatible with batch:
        # (;props, c) = ws[k]
        props = ws[k].props
        c = ws[k].c
        r = SA[rho, 0.0, z[k1]]
        
        # Think about this: we have to assume a planar wave so we have to select the 'best' incidence
        # angle. Here we use the angle with the first dipole
        costheta = sum(tli -> cosinc(pos(tli), r), tl) / length(tl)
        
        ef[k] = (exp(-M[k + 1] / costheta) * norm(electric_field(r, t, Ipeak, tl, 0.0, props, c)))
        
        if reduced
            ef[k] = ef[k] / ngas[k1] / co.Td
        end
    end
end    
