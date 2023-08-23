using WalkOnSpheres

# Problem definition
g(x,y) = [x;y]
∂Ω(x,y) = 1.0 - sqrt(x*x + y*y)
Ω = ((-1,1),(-1,1))
p = WalkOnSpheres.PDE(g,Ω,∂Ω)

# Solving
axs,u = WalkOnSpheres.ForwardSolve(p,ngrid=25,nwalks=200)

# Display
WalkOnSpheres.plot2D(axs,u)