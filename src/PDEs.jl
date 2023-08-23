"""
Contains field to define a laplace problem
# Arguments
- `g::Function`: value imposed on the domain's borders
- `∂Ω::Function`: domain definition (signed distance function)
"""
struct LaplacePDE <: PDE
    g::Function
    Ω::Tuple
    ∂Ω::Function
end

"""
Contains field to define a poisson problem
# Arguments
- `f::Function`: source term
- `g::Function`: value imposed on the domain's borders
- `∂Ω::Function`: domain definition (signed distance function)
"""
struct PoissonPDE <: PDE
    f::Function
    g::Function
    Ω::Tuple{Tuple{Float64,Float64}, Tuple{Float64,Float64}}
    ∂Ω::Function
end