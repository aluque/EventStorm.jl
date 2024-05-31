

function slow_reactions(T, fixed_dens)
    # From Kotovsky & Moore 2016
    F = 1.74e-18 + 1.93e-17 * sind(-7)^4
    ion_ion_recombination = 6e-8 * co.centi^3 * sqrt(300 / T)

    # Assuming thermal eq. for these slow times
    T_e = T
    T_i = T
    
    fix = map(((spec, v),) -> spec => j -> v[j], fixed_dens)
    
    rs = @kexprs ReactionSet([
        # From Kotovsky
        @withref(["Kotovsky2016", "Kotovsky2016[9]"],
                 "e + O3 -> O- + O2" => 1e-17,
                 "e + O3 -> O2- + O" => 1e-15),

        
        # From Gordillo-Vázquez 2016. I don't know where these come from; the reference in
        # Gordillo-Vázquez 2008 is to a measurement of cross-sections, which indeed gives
        # A larger cross section for ->O-.  
        # @withref("Gordillo-Vazquez2016",
        #          "e + O3 -> O- + O2" => 7.4e-18,
        #          "e + O3 -> O2- + O" => 1.24e-18),
        
        "e + O2 + O2 -> O2- + O2" =>
            (1.4e-41 * (300 / T_e) * exp(-600 / T) * exp(700 * (T_e - T) / (T_e * T))) .. "ref",
        
        "e + O2 + N2 -> O2- + N2" =>
            (1.07e-43 * (300 / T_e)^2 * exp(-70 / T) * exp(1500 * (T_e - T) / (T_e * T))) .. "ref",
        
        # Negative ions and electrons
        "O- + O2 -> e + O3" => 5.0e-21 .. ["Kossyi1994", "Kossyi1994[28]", "Kossyi1994[29]"],
        "O- + O -> e + O2" => 5.0e-16 .. ["Kossyi1994", "Kossyi1994[36]"],
        "O- + O3 -> O3- + O" => 5.30e-16 .. ["Kossyi1994", "Kossyi1994[20]", "Kossyi1994[25]"],
                
        "O2- + O -> e + O3" => 3.3e-16 .. ["Kotovsky2016", "Kotovsky2016[11]"],        
        
        "O2- + O -> O- + O2" => 3.30e-16 .. ["Kossyi1994", "Kossyi1994[25]", "Kossyi1994[30]"],
        "O2- + O3 -> O3- + O2" => 4e-16 .. ["Kossyi1994", "Kossyi1994[20]",
                                            "Kossyi1994[25]", "Kossyi1994[29]", "Kossyi1994[30]"],
        "O3- + O -> O2 + O2 + e" => 3e-16 .. ["Kossyi1994", "Kossyi1994[37]"],
        "O3- + O -> O2- + O2" => 2.5e-16 .. ["Kotovsky2016", "Kotovsky2016[13]"],
        
        

        # CO2- and H2O- based reactions
        @withref(["Kotovsky2016", "Kotovsky2016[13]"],
                 "O- + CO2 + M -> CO3- + M" => 2.0e-40,
                 "O2- + CO2 + M -> CO4- + M" => 4.7e-41,
                 "O3- + CO2 -> CO3- + O2" => 5.5e-16,
                 "CO3- + O -> O2- + CO2" => 1.1e-16,
                 "CO3- + O2 -> O3- + CO2" => 6.0e-21,
                 "O- + H2O -> OH- + OH" => 6.0e-19,
                 "OH- + O3 -> O3- + OH" => 9.0e-16,
                 "OH- + CO2 + M -> HCO3- + M" => 7.6e-40,
                 "OH- + O -> HO2 + e" => 2e-16),

        
        "O- + O2 + N2 -> O3- + N2" => (1.10e-42 * (300 / T)),
        "O- + O2 + O2 -> O3- + O2" => (1.10e-42 * (300 / T)),
        
        # Lyman-β
        "O2 + β -> e + O2+" => 0.90e-18 * co.centi^2,
        
        # Lyman-α
        "NO + α -> e + NO+" => 2.02e-18 * co.centi^2,
        
        # Constant rate of e/pos-ion recombination.  Must be improved
        # "e + N2+ -> " => 3e-13 * 300 / T_e + 4e-13 * (300 / T_e)^1.5 .. "Gordillo-Vazquez2016/JGR",
        # "e + O2+ -> " => 3e-13 * 300 / T_e + 4e-13 * (300 / T_e)^1.5 .. "Gordillo-Vazquez2016/JGR",
        # "e + NO+ -> " => 3e-13 * 300 / T_e + 4e-13 * (300 / T_e)^1.5 .. "Gordillo-Vazquez2016/JGR",
        @withref(["Kotovsky2016", "Sheehan2004", "Biondi1975", "Peterson1998"],
                 "e + N2+ -> N(4S) + N(2D)" => 0.37 * 2.2e-13 * (300 / T)^0.2 * (T_e / T)^-.39,
                 "e + N2+ -> N(2P) + N(4S)" => 0.11 * 2.2e-13 * (300 / T)^0.2 * (T_e / T)^-.39,
                 "e + N2+ -> N(2D) + N(2D)" => 0.52 * 2.2e-13 * (300 / T)^0.2 * (T_e / T)^-.39),
        
        @withref(["Kotovsky2016", "Sheehan2004", "Biondi1975", "Hellberg2003"],
                 "e + NO+ -> N(4S) + O(3P)" => 0.05 * 3.5e-13 * (300 / T) * (T_e / T)^-.69,
                 "e + NO+ -> N(4S) + O(3P)" => 0.95 * 3.5e-13 * (300 / T) * (T_e / T)^-.69),
        
        @withref(["Kotovsky2016", "Kotovsky2016[22]"],
                 "e + O4+ -> " => 4.2e-12 * (T_e / T)^-0.48),

        @withref(["Kotovsky2016", "Kotovsky2016[18]", "Kotovsky2016[19]"],
                 "e + O2+ -> " => 1.95e-13 * (300 / T)^0.7 * (T_e / T)^-0.7),

        @withref(["Kossyi1994"],
                 "N+ + O2 -> O2+ + N" => 2.8e-16,
                 "N+ + O2 -> NO+ + O" => 2.5e-16,
                 "N+ + O2 -> O+ + NO" => 2.8e-17,
                 "N+ + O -> O+ + N" => 1e-18,
                 "N+ + O3 -> NO+ + O2" => 5e-16,
                 "O+ + N2 -> NO+ + N" => 3e-18 * exp(-0.00311 * T),
                 "O+ + O2 -> O2+ + O" => 3.3e-17 * exp(-0.00169 * T),
                 "O+ + O3 -> O2+ + O2" => 1e-16,
                 "N2+ + O2 -> O2+ + N2" => 6e-17 * sqrt(300 / T),
                 "N2+ + O -> NO+ + N" => 1.3e-16 * sqrt(300 / T),
                 "N2+ + O3 -> O2+ + O + N2" => 1e-16,
                 "O2+ + N2 -> NO+ + NO" => 1e-23,
                 "O4+ + O2 -> O2+ + O2 + O2" => 3.3e-12 * (300 / T)^4 * exp(-5030 / T),
                 "O4+ + O -> O2+ + O3" => 3e-16),

        # Production of O4+
        "O2+ + O2 + M -> O4+ + M" => 2.6e-42 * (300 / T)^3.2 .. ["Kotovsky2016", "Brasseur1986"],

        # Hydration of positive ions and recombination with e-
        "O4+ + H2O -> Y+ + O2" => 1.5e-15 .. ["Kotovsky2016", "Brasseur1986"],
        "NO+ + M + M -> Y+ + M" => (exp(-45.17 + 0.11 * T - 4.84e-4 * T^2)
                                      .. ["Kotovsky2016", "Kotovsky2016[17]"]),
        "Y+ + e -> " => 3e-12 .. "Reid1977",
        
        # Positive ion chemistry
        # 3-body processes; slow for upper-atmosphere
        # "N2+ + N2 + N2 => N4+ + N2", 5e-41                       .. ["Kossyi1994", "Kossyi1994[20]"],
        # "N2+ + N + N2 => N3+ + N2", 0.9e-41 * exp(400 / T)       .. ["Kossyi1994", "Kossyi1994[40]"],
        # "O2+ + O2 + O2 => O4+ + O2", 2.4e-42 * (300 / T)^3.2     .. ["Kossyi1994", "Kossyi1994[24]"],
        # "O2+ + N2 + N2 => O2+.N2 + N2", 0.9e-42 * (300 / T)^2    .. ["Kossyi1994", "Kossyi1994[20]"],
        
        
        # Production rate of electrons from cosmic rays. 
        "O2 -> e + O2+" =>  F .. "ref",
        "N2 -> e + N2+" =>  F .. "ref",
        
        cartesian_product(["O-", "O2-", "O3-", "CO3-", "OH-", "HCO3-"],
                          ["N2+", "O2+", "NO+", "O4+"], "",
                          ion_ion_recombination .. ["Kotovsky2016", "Smith1977"])...,
        
        # Phony reaction to keep track of second positive emissions
        "SPS -> SPS" => 1.0
    ]; fix)
    
    return rs
end


function cartesian_product(spec1, spec2, res, rate)
    return vec(["$s1 + $s2 -> $res" => rate for s1 in spec1, s2 in spec2])
end


function fast_reactions(T, comp; ngas, swarmfile=joinpath(DATA_DIR, "swarm/LxCat_Phelps_20230914.txt"))
    lx = LxCatSwarmData.load(swarmfile)
    tbl = loadtable(eachcol(lx.data), xcol=:en)

    rec = loadtable(eachcol(CSV.read(joinpath(DATA_DIR, "swarm/rec_electron.dat"), DataFrame,
                                     header=[:en, :k])), xcol=:en)
    
    compN2 = comp["N2"]
    compO2 = comp["O2"]
    ion_ion_recombination = 1e-7 * co.centi^3
    
    frs = ReactionSet([
        "e + N2 -> 2 * e + N2+" => RateLookup(tbl, :C25),
        "e + N2 -> 2 * e + N2+" => RateLookup(tbl, :C26),
        "e + O2 -> 2 * e + O2+" => RateLookup(tbl, :C43),
        "e + O2 -> O + O-" => RateLookup(tbl, :C28),
        
        # We could have used O2 instead of O2M and then dividing the rate by 1e6
        "e + O2 + O2M -> O2- + O2M" => RateLookup(tbl, :C27),
        
        @withref(["Pancheshnyi2013/JPhD"],
                 "O2- + M -> e + O2 + M" => PancheshnyiFitEN(1.24e-11 * co.centi^3, 179.0, 8.8),
                 "O- + O2 -> O2- + O" => PancheshnyiFitEN(6.96e-11 * co.centi^3, 198.0, 5.6),
                 "O2 + O- + M -> O3- + M" => PancheshnyiFitEN2(1.1e-30 * co.centi^6, 65.0)),
        
        @withref(["Schuman2023/PhysChem"],
                 "N2 + O- -> e + N2O" => ModArrhenius(;k0=3.98e-17, T0=5097.0, d=-1.36, T),
                 "N2(v1) + O- -> e + N2O" => ModArrhenius(;k0=9.04e-18, T0=674.0, d=-0.85, T),
                 "N2(v2) + O- -> e + N2O" => ModArrhenius(;k0=2.74e-17, T0=186.0, d=-1.10, T)),
        
        @withref(["Lawton1978/JChPh", "Phelps1985/PhRvA", "Hagelaar2005/PSST"],
                 # Resonant v1
                 "e + N2 -> e + N2(v1)" => RateLookup(tbl, :C3),
                 "e + N2 -> e + N2(v1)" => RateLookup(tbl, :C4),
                 "e + N2 -> e + N2(v2)" => RateLookup(tbl, :C5),
                 "e + N2 -> e + N2(v3)" => RateLookup(tbl, :C6)),
        
        @withref(["Aleksandrov1999/PSST"],
                 "N2+ + N2 + M -> N4+ + M" => TemperaturePower(;k0=5e-29 * co.centi^6, power=2, T),
                 "N4+ + O2 -> 2 * N2 + O2+" => 2.5e-10 * co.centi^3,
                 "O2+ + O2 + M -> O4+ + M" => TemperaturePower(;k0=2.4e-30 * co.centi^6, power=3, T)),

        "e + O4+ -> O2 + O2" => (RateLookup(rec, :k) .. "Kossyi1992/PSST"),
        
        cartesian_product(["O-", "O2-"],
                          ["N2+", "N4+", "O2+", "O4+"], "",
                          ion_ion_recombination .. ["Kossyi1992/PSST"])...,
        
        # Instantaneous emission
        "e + N2 -> e + SPS" => RateLookup(tbl, :C10)
    ]; fix = [:N2 => j -> compN2 * ngas[j],
              :O2 => j -> compO2 * ngas[j],
              
              :M => j -> ngas[j],
              
              # The 1e6 comes from Bolsig+'s normalization of 3-body
              :O2M => j -> compO2 * ngas[j] / 1e6])
    
end


@kwdef struct TemperaturePower{TY,P}
    k0::TY
    power::P
    T::TY
    T0::TY = 300.0
end

Chemise.evalk(p::TemperaturePower, en) = p.k0 * (p.T0 / p.T)^p.power

@kwdef struct PancheshnyiFitEN{T}
    k0::T
    a::T
    b::T
end

Chemise.evalk(p::PancheshnyiFitEN, en) = p.k0 * exp(-(p.a / (p.b + en))^2)


@kwdef struct PancheshnyiFitEN2{T}
    k0::T
    a::T
end

Chemise.evalk(p::PancheshnyiFitEN2, en) = p.k0 * exp(-(en / p.a)^2)


@kwdef struct ModArrhenius{TY}
    k0::TY
    T0::TY
    d::TY
    T::TY
    mgas::TY = 28.02 * co.gram / co.Avogadro
    K0::TY = 4.5 * co.centi^2
end

function Chemise.evalk(p::ModArrhenius, en)
    vd = co.nair * p.K0 * en * co.Td
    Teff = p.T + p.mgas * vd^2 / (3 * co.k)
    k = p.k0 * (Teff / 300)^p.d * exp(-p.T0 / Teff)
end



# The electron transport and ionization parameters are defined for a given
# reference gas density (300 K, 1 bar).
const NGAS_REF = 2.4143235045419996e+25

const E_MOBILITY = 0.0372193 * NGAS_REF
const IONIZATION_ALPHA = 433200.0 / NGAS_REF
const IONIZATION_FIELD = 2e7 / NGAS_REF
const ATTACHMENT_ALPHA = 2000.0 / NGAS_REF
const ATTACHMENT_FIELD = 3e6 / NGAS_REF

function fast_reactions_townsend(comp; ngas)
    frs = ReactionSet(["e + M -> 2 * e + M" => en -> impact_nu(NGAS_REF, en * co.Td * NGAS_REF) /  NGAS_REF,
                       "e + M -> M" => en -> attachment_nu(NGAS_REF, en * co.Td * NGAS_REF) /  NGAS_REF,
                       ];
                      fix = [:M => j -> ngas[j]])

end



impact_nu(ngas, eabs) = (E_MOBILITY * eabs * IONIZATION_ALPHA *
                         exp(-ngas * IONIZATION_FIELD / eabs))
attachment_nu(ngas, eabs) = (E_MOBILITY * eabs * ATTACHMENT_ALPHA *
                             exp(-ngas * ATTACHMENT_FIELD / eabs))

