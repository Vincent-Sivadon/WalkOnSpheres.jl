# Green's function 2D
function greens(x1,y1,x2,y2,d)
    r = √((x2-x1)^2 + (y2-y1)^2)
    log(d/r) / 2π
end
# Green's function 3D
function greens(x1,y1,z1,x2,y2,z2,d)
    r = √((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2)
    (d-r) / (4π * r*d)
end


# Pick random cartesian coords on ball of dim 2 (disk)
function rand_on_ball(cx,cy,r)
    θ = 2π * rand()
    cx + r*cos(θ), cy + r*sin(θ)
end
# Pick random cartesian coords on ball of dim 3 (sphere)
function rand_on_ball(cx,cy,cz,r)
    θ = acos(2*rand() - 1)
    ϕ = 2π * rand()
    cx+r*sin(θ)*cos(ϕ), cy+r*sin(θ)*sin(ϕ), cz+r*cos(θ)
end

# Pick random cartesian coords in ball of dim 2 (disk)
function rand_in_ball(cx,cy,r)
    ρ = r * sqrt(rand())
    θ = 2π * rand()
    cx+ρ*cos(θ), cy+ρ*sin(θ)
end
# Pick random cartesian coords in ball of dim 3 (sphere)
function rand_in_ball(cx,cy,cz,r)
    θ = acos(2*rand() - 1)
    ϕ = 2π * rand()
    ρ = r*cbρt(ρand())
    cx+ρ*sin(θ)*cos(ϕ), cy+ρ*sin(θ)*sin(ϕ), cz+ρ*cos(θ)
end


"Makes the problem adaptable to all dimensional cases (Input:{1D,2D...}, Output:{Scalar,Vector2...})"
function preprocess(p::PDE,ngrid)
    # Input dimension (1D, 2D, 3D... ?)
    # Ex : Ω = ((-1,1),(2,3)) means x goes from -1 to 1 and y from 2 to 3
    # So it's a 2D input problem and we have 2 axes :
    #   - range(-1,1;length=ngrid)
    #   - range(2 ,3;length=ngrid)
    idim = length(p.Ω)
    axs = [range(p.Ω[i]...;length=ngrid) for i ∈ 1:idim]

    # Output dimension (The solution is a scalar field, vector field 2D...)
    idims = ntuple(x->ngrid,idim)   # ex : (ngrid,ngrid) for 2D
    idx  = ntuple(x->0,idim)        # ex : (0,0) in 2D (will be given to g as input)
    odim = length(p.g(idx...))      # ex : the size of g's output tells us the output dimension
    if odim==1 # Scalar field
        u = Array{Float64,idim}(undef, idims...)
    else    # Vector field of length odim
        u = Array{Vector{Float64},idim}(undef, idims...)
        fill!(u,Vector{Float64}(undef,odim))
    end

    idim,odim,axs,u
end
function preprocess_withgrad(p::PDE,ngrid)
    idim,odim,axs,u = preprocess(p,ngrid)
    idims = ntuple(x->ngrid,idim)

    if odim==1
        ∇u = Array{Vector{Float64},idim}(undef,idims...)
        fill!(∇u,Vector{Float64}(undef,idim))
    else
        ∇u = Array{Matrix{Float64},idim}(undef,idims...)
        fill!(∇u,Matrix{Float64}(undef,odim,idim))
    end

    idim,odim,axs,u,∇u
end

"""
    Return axes values at index

Ex : input 2D
`axs` are x = range(...) and y = range(...)
`idx` is for example (3,4)
we return x[3], y[4]
"""
function getcoords(axs::Vector,idx,idim)
    ntuple(i->axs[i][idx[i]],idim)
end