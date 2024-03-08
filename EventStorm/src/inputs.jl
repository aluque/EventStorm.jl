#=
This is generic code to handle inputs.  The idea is that you define the input parameters
and their default values simply by adding kwargs to your `main` function (or whatever you want to call
it). The you can use `wrap_main(main, inputfile)` to read `inputfile` and update whatever parameters
are contained there. Now inputfile can be either julia (recommended but unsafe if you do not trust
the source) or toml (safe but inconvenient). Note also that TOML cannot handle julia-specific argument
types (anything other than numbers, strings or booleans), which are forced to take default values.

Also functions writejl and writetoml are provided to write the actually used parameters, for
future reference.
=#
module WrapInput; end

using TOML
using Dates
using LibGit2

function wrap_input(f, input::String; pretty_print=false)
    _, ext = splitext(input)
    if ext == ".jl"
        ntpl = readjl(f, input)
    elseif ext == ".toml"
        ntpl = readtoml(f, input)
    end

    if pretty_print
        io = IOBuffer()
        pprint(IOContext(io, :color => true), ntpl)
        @info "Input read from $input" * "\n" * String(take!(io))
    end
    
    f(;ntpl...)
end

"""
Read a julia input file and produces a `NamedTuple` with the values defined there that are also
kwargs of the function `f`.
"""
function readjl(f, input::String)
    @eval WrapInput include(abspath($input))

    #@assert (m isa Module) "The input parameters must be contained in a `module`"
    @assert length(methods(f)) == 1 "Ambiguity: more than one method defined for _main"

    method = first(methods(f))
    
    p = Pair{Symbol, Any}[]
    
    input_vars = Base.kwarg_decl(method)
    for sym in names(WrapInput, all=true)
        if sym in input_vars
            push!(p, sym => getfield(WrapInput, sym))
        end
    end

    if :_input in input_vars
        push!(p, :_input => abspath(input))
    end

    if :_date in input_vars
        push!(p, :_date => now())
    end

    if :_git_commit in input_vars
        hash, dirty = git_commit()
        push!(p, :_git_commit => hash)
        push!(p, :_git_dirty => dirty)
    end
        
    return NamedTuple(p)
end


"""
Read a julia input file and produces a `NamedTuple` with the values defined there that are also
kwargs of the function `f`.
"""
function readtoml(f, input::String)
    params = TOML.parsefile(input)

    p = Pair{Symbol, Any}[]
    input_vars = Base.kwarg_decl(method)
    for (k, v) in params
        if Symbol(k) in input_vars
            push!(p, Symbol(k) => v)
        end
    end

    if :_input in input_vars
        push!(p, :_input => abspath(input))
    end

    if :_date in input_vars
        push!(p, :_date => now())
    end    

    if :_git_commit in input_vars
        hash, dirty = git_commit()
        push!(p, :_git_commit => hash)
        push!(p, :_git_dirty => dirty)
    end

    return NamedTuple(p)
end


"""
Writes into a file a julia expression with the kwargs of the function `f` contained in params.
Typically params will be obtained from `Base.@locals`.
"""
function writejl(io::IO, f, params)
    expr = quote end
    @assert length(methods(_main)) == 1 "Ambiguity: more than one method defined for _main"

    method = first(methods(_main))
    input_vars = Base.kwarg_decl(method)

    for (sym, val) in pairs(params)
        if sym in input_vars
            push!(expr.args, :($sym = $val))
        end
    end
    print(io, expr)
end


"""
Calls `writejl(io,...)` opening the file with name `fname`.
"""
function writejl(fname::String, args...)
    open(fname, "w") do fout
        writejl(fout, args...)
    end
end


"""
Finds the root directory of a Git repository given a starting directory.
"""
function find_git_root(start_dir::String)
    current_dir = start_dir
    while current_dir != "/"
        if isdir(joinpath(current_dir, ".git"))
            return current_dir
        end
        current_dir = dirname(current_dir)
    end
    return nothing # No Git repository found
end

"""
Returns the latest git commit and a boolean telling whether the repo is dirty.
"""
function git_commit(repo::GitRepo)
    commit = string(LibGit2.GitHash(LibGit2.peel(LibGit2.GitCommit, LibGit2.head(repo))))
    dirty = LibGit2.isdirty(repo)

    return commit, dirty
end

git_commit() = git_commit(GitRepo(find_git_root(@__DIR__)))


function pprint(io::IO, d, level=0)
    for (k,v) in pairs(d)
        if typeof(v) <: Dict
            print(io, join(fill(" ", level * 8)))
            printstyled(io, k, color=:light_yellow, bold=true)            
            printstyled(io, " => \n", color=:light_black)            
            pretty_print(io, v, level + 1)
        else
            print(io, join(fill(" ", level * 8)))
            printstyled(io, k, color=:light_yellow, bold=true)
            printstyled(io, " => ", color=:light_black)
            printstyled(io, repr(v) * "\n", color=:blue) 
        end
    end
    nothing
end

pretty_print(d::Dict, kw...) = pretty_print(stdout, d, kw...)
