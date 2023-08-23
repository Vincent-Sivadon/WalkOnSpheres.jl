"""
    Solve a forward Walk On Spheres problem (Laplace, Poisson) with generic dimensions
    Returns a vector of the axis ranges `axs` (axs[1] -> xrange...) and the solution `u`
"""
function ForwardSolve(p::PDE; ngrid=50, nwalks=200, ϵ=0.001)
    idim,odim,axs,u = preprocess(p,ngrid)

    # Indices ex 2D with ngrid = 2
    # indices[:] = (1,1), (2,1), (1,2), (2,2)
    ranges = ntuple(x->1:ngrid, idim)
    indices = Iterators.product(ranges...)

    # Solve for each element of u
    for idx in indices
        # Coords associated with idx
        coords = getcoords(axs,idx,idim)

        # Solve at coords
        u[idx...] = walks(p,coords,nwalks,odim,ϵ)
    end
    axs,u
end

"Same as the Solve function but on GPU"
function ForwardSolveGPU(p::PDE; ngrid=50, nwalks=200, ϵ=0.001)
    idim,odim,axs,u = preprocess(p,ngrid)
    
    if idim==2 && odim==1 # 2D scalar
        return GPU_2Dscalar(p,axs,u,nwalks,ϵ)
    else
        println("WalkOnSphereSim Error : can't handle those dimensions")
        return
    end
end

"Return the solution of the laplace problem `p` for position (`xᵢ`,`yᵢ`) with `N` walks"
function walks(p::PDE,coords,nwalks,odim,ϵ)
    uᵢ = odim==1 ? 0.0 : zeros(odim)
    for _=1:nwalks
        uᵢ += walk(p,coords,odim,ϵ)
    end
    uᵢ /= nwalks
end