using WalkOnSpheres

# Problem definition
g(x,y,z) = z + 2*sin(atan(y,x))
∂Ω(x,y,z) = 1.0 - sqrt(x*x + y*y + z*z)
Ω = ((-1,1),(-1,1),(-1,1))
p = WalkOnSpheres.PDE(g,Ω,∂Ω)

# Solving
axs,u = WalkOnSpheres.ForwardSolve(p,ngrid=20,nwalks=200)

# Display
WalkOnSpheres.plot3D(axs,u)