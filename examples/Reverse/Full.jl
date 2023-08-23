using WalkOnSpheres

# Problem definition
f(x,y) = 20.
g(x,y) = 0.
h(x,y) = 0.4
Ω = ((-1,1),(-1,1))
∂Ω(x,y) = 1.0 - sqrt(x^2+y^2)
is∂Ωn(x,y) = x>0.8
function sample_s()
    ρ = 0.1 * sqrt(rand())
    θ = 2π * rand()
    (ρ*cos(θ),ρ*sin(θ)),1/(π*0.1^2)
end
function sample_b()
    θ = 2π*rand()
    (cos(θ),sin(θ)),1/(2π)
end

ngrid = 200
nwalks = 1e5
axs,u = WalkOnSpheres.ReverseSolveGPU(f,g,h,Ω,∂Ω,is∂Ωn,sample_s,sample_b,ngrid=ngrid,nwalks=nwalks)

WalkOnSpheres.plot2D(axs,u,figure=(;resolution=(700,700)))