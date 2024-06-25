## Tue Jun 25 16:36:26 2024 : storm-baseline.jl
## This is a baseline simulation with the same parameters as storm.jl but without any flashes.
##
using ..EventStorm: DATA_DIR, co

# Resolution in points / km
points_per_km = 10

# Time between outputs
output_dt = 30

# Time in the fast model after the end of the source
extra_time = 50e-3

#
# ATMOSPHERE
#
gas_density_fname = joinpath(DATA_DIR, "earth", "stdatm.csv")
iri_electron_density_fname = joinpath(DATA_DIR, "iri", "ingrid.txt")

#
# FLASH DISTRIBUTION
#
# log-normal distribution for the peak currents from Ingrid's paper (see fit_ingrid.jl)
Ipeak_median = 20.36 * co.kilo
Ipeak_log_std = 1.14
Ipeak_cutoff = 400 * co.kilo

# The storm is simulated as gaussian in time and both spatial dimensions;
# The duration extension here are standard deviations.
storm_duration = 1 * co.hour
storm_extension = 100 * co.kilo
storm_center_time = 1.5 * co.hour

# Expected number of flashes in the storm (actual number may differ).
storm_expected_flashes = 0.0

# Distance of center of the storm to the observer.
storm_distance = 100 * co.kilo

#
# FLASH CURRENT DISTRIBUTION
#
# Discharge parameters
rise = 2 * co.micro
decay = 10 * co.micro
plateau = 0.0

