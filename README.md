# EventStorm.jl

Supporting code for the manuscript "Cumulative effects of lightning electromagnetic pulses on the lower ionosphere" by A. Luque _et al._

## Installation

```julia
julia> using Pkg
julia> Pkg.add(url="https://github.com/aluque/Constants.jl")
julia> Pkg.add(url="https://github.com/aluque/Chemise.jl")
julia> Pkg.add(url="https://github.com/aluque/DipoleRadiators.jl")
julia> Pkg.add(url="https://github.com/aluque/EventStorm.jl")
```

## Running
The folder `utils/` contains a script called `run.jl`. Use this script to run a simulation using

```bash
julia run.jl input-file.jl
```

Alternatively, from the julia prompt:
```julia
julia> using EventStorm
julia> EventStorm.run_from_input("input-file.jl")
```

## Input files
The folder `input/` contains input the files used in the manuscript.
