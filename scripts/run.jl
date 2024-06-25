using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using EventStorm: _main, wrap_input
wrap_input(_main, ARGS[1], pretty_print=true)
