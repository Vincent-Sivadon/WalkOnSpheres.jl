using WalkOnSpheres

# Problem definition
g(x,y,z) = [x;y;z]
∂Ω(x,y,z) = 1.0 - sqrt(x*x + y*y + z*z)
Ω = ((-1,1),(-1,1),(-1,1))
p = WalkOnSpheres.PDE(g,Ω,∂Ω)

# Solving
ngrid=15
axs,u = WalkOnSpheres.ForwardSolve(p,ngrid=ngrid,nwalks=100)

# Display
WalkOnSpheres.plot3D(axs,u)