using WalkOnSpheres

# PDE definition
g(x,y) = 1. + 2*sin(atan(y,x))
∂Ω(x,y) = 1.0 - sqrt(x*x + y*y)
Ω = ((-1,1),(-1,1))
p = WalkOnSpheres.PDE(g,Ω,∂Ω)

# Solving
axs,u = WalkOnSpheres.ForwardSolveGPU(p,ngrid=20,nwalks=200)

# Display
WalkOnSpheres.plot2D(axs,u)