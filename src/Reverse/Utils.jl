function normal(x,y,∂Ω)
    nx = ∂Ω(x-1e-2,y) - ∂Ω(x+1e-2,y)
    ny = ∂Ω(x,y-1e-2) - ∂Ω(x,y+1e-2)
    norm = √(nx^2 + ny^2)
    nx/norm,ny/norm
end

function concatenatewalks(walks)
    walks_x = Float32[]
    walks_y = Float32[]
    walks_h = Float32[]
    walks_r = Float32[]
    for nt ∈ 1:Threads.nthreads()
        walks_x = vcat(walks_x,walks[nt].x)
        walks_y = vcat(walks_y,walks[nt].y)
        walks_h = vcat(walks_h,walks[nt].h)
        walks_r = vcat(walks_r,walks[nt].r)
    end
    WalksData(walks_x,walks_y,walks_h,walks_r)
end