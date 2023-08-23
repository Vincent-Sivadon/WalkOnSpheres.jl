using WalkOnSpheres

# Problem definition
@inline ρ(x,y) = sqrt(x*x+y*y)
@inline ρ1(x,y) = sqrt((x+2.)^2+y^2)
@inline ρ2(x,y) = sqrt((x-2.)^2+y^2)
const r∞ = 10.0         # Infinite domain radius
const rw = 1.0          # Wire radius
const l  = 100.0        # Wire length
const μ₀ = 1.0          # Magnetic permeability
const I  = 1.0          # Current intensity
const j = I / (π*rw^2)  # Current density    
Ω = ((-2*r∞/3,2*r∞/3),(-2*r∞/3,2*r∞/3)) # Domain
∂Ω(x,y) = r∞ - sqrt(x*x+y*y)
# function ∂Ω(x,y)
#     vx = abs(x)-r∞ ; vy = abs(y)-r∞ ;
#     dx = max(vx,0)  ; dy = max(vy,0)
#     return - (sqrt(dx*dx + dy*dy) + min(max(vx,vy),0))
# end
@inline g(x,y) = 0.
function s(x,y)
    if ρ1(x,y) <= rw
        return μ₀*j
    elseif ρ2(x,y) <= rw
        return -μ₀*j
    else
        return 0
    end
end
p = WalkOnSpheres.PDE(s,g,Ω,∂Ω)

# Solve
ngrid = 50
axs,A = WalkOnSpheres.ForwardSolveGPU(p,ngrid=ngrid,nwalks=1000)

# # Display
fig = WalkOnSpheres.plot2D(axs,A)
fig

# Analytical Solution and tests
# σ = StandardDeviation_2Dscalar(A,ua,axs)
# fig, = WalkOnSpheres.plot2D(axs,σ)
# fig

# # Analytical solution
# using CairoMakie
# @inline ua(x,y) = - μ₀*I/(2π) * (log(ρ1(x,y)) - log(ρ2(x,y)))
# x = range(-r∞,r∞;length=50)
# y = range(-r∞,r∞;length=50)
# figa,ax,hm = heatmap(
#     x,y,ua;
#     figure = (; resolution=(700,700)),
#     axis = (; title="Infinite Wire Solution", xlabel="x", ylabel="y",aspect=1)
# )
# Colorbar(figa[1,2],hm)
# figa