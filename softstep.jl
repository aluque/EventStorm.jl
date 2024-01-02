
"""
    A representation of an axial source with soft-step polynomial functions
"""
Base.@kwdef struct SoftStep{T}
    "Lowest altitude bound."
    r1::SVector{3, T}

    "Upper altitude bound."
    r2::SVector{3, T}

    "Peak current"
    Ipeak::T

    "Rise time"
    rise::T

    "decay"
    decay::T

    "plateau"
    plateau::T=zero(decay)
end


softstep(x) = x^2 * (2 - x)^2

"""
    Current (cross-sectional, in A) at time `t`.
"""
function current(source::SoftStep, t)
    (;rise, decay, Ipeak, plateau) = source

    if t < rise
        return Ipeak * softstep(t / rise)
    elseif t < rise + plateau
        return Ipeak
    elseif t < rise + plateau + decay 
        return Ipeak * softstep(1 - (t - rise - plateau) / decay)
    else
        return zero(Ipeak)
    end
end
