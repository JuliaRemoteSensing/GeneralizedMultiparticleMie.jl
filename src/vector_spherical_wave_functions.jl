@inline rad_hat(θ, ϕ) =
    StaticArrays.@SVector [sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ)]

@inline theta_hat(θ, ϕ) =
    StaticArrays.@SVector [cos(θ) * cos(ϕ), cos(θ) * sin(ϕ), -sin(θ)]

@inline phi_hat(_, ϕ) =
    StaticArrays.@SVector [-sin(ϕ), cos(ϕ), 0]

@inline sph_to_cart(v, θ, ϕ) =
    StaticArrays.@SVector [
        v ⋅ rad_hat(θ, ϕ),
        v ⋅ theta_hat(θ, ϕ),
        v ⋅ phi_hat(θ, ϕ),
    ]

@inline Emn(m, n) = 
    (1im)^(n%4)*√((2n+1)*factorial(n-m)/(n*(n + 1)*factorial(n+m)))