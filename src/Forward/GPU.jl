# Preprocess 2D scalar problem
function PreprocessGPU_2Dscalar(p::PDE,axs,u,nwalks,ϵ)
    ngrid = size(u)[1]
    x_d = CuArray(collect(Float32,axs[1]))
    y_d = CuArray(collect(Float32,axs[2]))
    u_d = CuArray(Float32.(u))

    threads = (16,16)
    blocks = ceil.(Int64, (ngrid,ngrid) ./ threads)
    x_d,y_d,u_d,threads,blocks
end

# Run Laplace 2D scalar GPU solver
function GPU_2Dscalar(p::LaplacePDE,axs,u,nwalks,ϵ)
    x_d,y_d,u_d,threads,blocks = PreprocessGPU_2Dscalar(p,axs,u,nwalks,ϵ)
    @cuda threads=threads blocks=blocks Laplace_2Dscalar_kernel(u_d,x_d,y_d,p.g,p.∂Ω,nwalks,ϵ)
    axs,Array(u_d)
end
# Run Poisson 2D scalar GPU solver
function GPU_2Dscalar(p::PoissonPDE,axs,u,nwalks,ϵ)
    x_d,y_d,u_d,threads,blocks = PreprocessGPU_2Dscalar(p,axs,u,nwalks,ϵ)
    @cuda threads=threads blocks=blocks Poisson_2Dscalar_kernel(u_d,x_d,y_d,p.f,p.g,p.∂Ω,nwalks,ϵ)
    axs,Array(u_d)
end


# Laplace 2D scalar GPU kernel
function Laplace_2Dscalar_kernel(u,x,y,g,∂Ω,nwalks,ϵ)
    # Thread id
    ix = (blockIdx().x - 0x1) * blockDim().x + threadIdx().x
    iy = (blockIdx().y - 0x1) * blockDim().y + threadIdx().y
    
    if ix <= size(u)[1] && iy <= size(u)[2]
        uᵢⱼ = 0.0f0
        for _ ∈ 1:nwalks
            # Reset xᵢ,yᵢ
            @inbounds xᵢ = x[ix]
            @inbounds yᵢ = y[iy]
            # Walk till edge
            while true
                d = ∂Ω(xᵢ,yᵢ)
                d<0 && break 
                d<ϵ && (uᵢⱼ += g(xᵢ,yᵢ) ; break )
                theta = 2π * rand(Float32)
                xᵢ += d*cos(theta)
                yᵢ += d*sin(theta)
            end
        end
        @inbounds u[ix,iy] = uᵢⱼ / nwalks
    end

    return
end
# Laplace 2D scalar GPU kernel
function Poisson_2Dscalar_kernel(u,x,y,f,g,∂Ω,nwalks,ϵ)
    # Thread id
    ix = (blockIdx().x - 0x1) * blockDim().x + threadIdx().x
    iy = (blockIdx().y - 0x1) * blockDim().y + threadIdx().y
    
    if ix <= size(u)[1] && iy <= size(u)[2]
        uᵢⱼ = 0.0f0
        for _ ∈ 1:nwalks
            # Reset xᵢ,yᵢ
            @inbounds xᵢ = x[ix]
            @inbounds yᵢ = y[iy]
            # Walk till edge
            uₛ = 0.0f0
            while true
                d = ∂Ω(xᵢ,yᵢ)
                d<0 && break 
                d<ϵ && (uᵢⱼ += g(xᵢ,yᵢ) + uₛ ; break )

                # New inx,iny
                ρ = d * sqrt(rand(Float32))
                θ = 2π * rand(Float32)
                inx = xᵢ + ρ*cos(θ)
                iny = yᵢ + ρ*sin(θ)
                greens2D = log(abs(d) / sqrt((inx-xᵢ)^2 + (iny-yᵢ)^2)) / 2π
                uₛ += π*d^2 * f(inx,iny) * greens2D

                θ = 2π * rand(Float32)
                xᵢ += d*cos(θ)
                yᵢ += d*sin(θ)
            end
        end
        @inbounds u[ix,iy] = uᵢⱼ / nwalks
    end

    return
end