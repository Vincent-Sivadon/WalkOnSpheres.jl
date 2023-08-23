using WalkOnSpheres

# Problem Defintion
f(x,y) = 1.  # Uniform heating
g(x,y) = 0.  # Cold borders
function ∂Ω(x,y)
    vx = abs(x)-1 ; vy = abs(y)-1 ;
    dx = max(vx,0)  ; dy = max(vy,0)
    - (sqrt(dx*dx + dy*dy) + min(max(vx,vy),0))
end
Ω = ((-1,1),(-1,1))
p = WalkOnSpheres.PDE(f,g,Ω,∂Ω)

# Solving
axs,u = WalkOnSpheres.ForwardSolveGPU(p,ngrid=20,nwalks=200)

# Display
WalkOnSpheres.plot2D(axs,u)

# Analytical Solution
# function ua(x,y,n_modes=1)
#     sol = (1-x^2)/2
#     for k=1:n_modes
#         if (k%2 == 0) continue; end
#         sol -= 16/π^3 *
#             sin(k*π*(1+x)/2)/(k^3*sinh(k*π)) *
#             (sinh(k*π*(1+y)/2) + sinh(k*π*(1-y)/2))
#     end
#     sol
# end
# using CairoMakie
# heatmap(axs...,ua,aspect_ratio=1)