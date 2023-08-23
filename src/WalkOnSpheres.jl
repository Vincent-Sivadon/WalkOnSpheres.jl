module WalkOnSpheres

using Revise
using LinearAlgebra
using CUDA
using CairoMakie
using GLMakie

# Interface to define Partial Differential Equations (Laplace or Poisson)
abstract type PDE end

"""
Data structure to store reverse walk informations
# Arguments
- `x::Vector{Float32}`: x coordinate of the sphere's center
- `y::Vector{Float32}` : y coordinate of the sphere's center
- `h::Vector{Float32}` : value (depends on the walk type)
- `r::Vector{Float32}` : sphere radius
"""
struct WalksData
    x::Vector{Float32}
    y::Vector{Float32}
    h::Vector{Float32}
    r::Vector{Float32}
end

include("PDEs.jl")

# Function that dispatch on the right type of pde
PDE(g,Ω,∂Ω)   = LaplacePDE(g,Ω,∂Ω)
PDE(f,g,Ω,∂Ω) = PoissonPDE(f,g,Ω,∂Ω)

include("Forward/Utils.jl")
include("Forward/Laplace.jl")
include("Forward/Poisson.jl")
include("Forward/Solve.jl")
include("Forward/GPU.jl")

include("Reverse/Utils.jl")
include("Reverse/Source.jl")
include("Reverse/Boundary.jl")
include("Reverse/Solve.jl")
include("Reverse/GPU.jl")

include("PlotUtils/Plot2D.jl")
include("PlotUtils/Plot3D.jl")

end # module WalkOnSpheres
