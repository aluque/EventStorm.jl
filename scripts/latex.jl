# latex.jl : Thu May 23 12:36:45 2024
"""
# Latex
Produce LaTeX output from the chemical models

## Running the code

```julia
julia> includet("latex.jl")     # Needs Revise.jl
julia> Latex.main()
```
"""
module Latex
using EventStorm
using Chemise

function main(template="empstorm_si")
    out = EventStorm._main(;run=false)
    # kw = (preamble = "\\rowcolors{2}{orange!10}{white}",
    #       preheader = "\\rowcolor{orange!20}")
    
    open(joinpath(@__DIR__, "fast_chem.inc"), "w") do f
        writelatex(f, out.frs)
    end

    open(joinpath(@__DIR__, "slow_chem.inc"), "w") do f
        writelatex(f, out.rs)
    end    

    cd(@__DIR__) do 
        run(`pdflatex $(template).tex`)
        run(`bibtex $(template)`)
        run(`pdflatex $(template).tex`)
        run(`pdflatex $(template).tex`)
    end

    return NamedTuple(Base.@locals)
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    Latex.main()
end

