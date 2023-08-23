

"Return the one-walk estimator for the laplace problem `p` for position (`xᵢ`,`yᵢ`)"
function walk(p::LaplacePDE,coords,odim,ϵ)
    while true
        # Closest point to borders
        d = p.∂Ω(coords...)

        # If d < 0 then out of bounds
        d < 0 && return (odim==1) ? 0. : zeros(odim)

        # If Distance small enough : return value
        d<ϵ && return p.g(coords...)

        # Else, take a random pos on the circle of
        # radius d and center xᵢ,yᵢ
        coords = rand_on_ball(coords...,d)
    end
end