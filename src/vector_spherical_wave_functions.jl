@enum VSWFMode Incident = 1 Interior Outgoing Ingoing

const VSWF_KERNEL = (spherical_jn, spherical_yn, spherical_hn1, spherical_hn2)
const VSWF_KERNEL_DERIV = (spherical_jn_deriv, spherical_yn_deriv, spherical_hn1_deriv, spherical_hn2_deriv)
const VSWF_KERNEL_RICATTI = (ricatti_jn, ricatti_yn, ricatti_hn1, ricatti_hn2)
const VSWF_KERNEL_RICATTI_DERIV = (ricatti_jn_deriv, ricatti_yn_deriv, ricatti_hn1_deriv, ricatti_hn2_deriv)

struct VSWFCache{T<:Real,V<:AbstractVector}
    factor::T
    A::V
    B::V
end

@inline rad_hat(θ, ϕ) = StaticArrays.@SVector [sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ)]

@inline theta_hat(θ, ϕ) = StaticArrays.@SVector [cos(θ) * cos(ϕ), cos(θ) * sin(ϕ), -sin(θ)]

@inline phi_hat(_, ϕ) = StaticArrays.@SVector [-sin(ϕ), cos(ϕ), 0]

@inline sph_to_cart(v, θ, ϕ) = StaticArrays.@SVector [v ⋅ rad_hat(θ, ϕ),
                                                      v ⋅ theta_hat(θ, ϕ),
                                                      v ⋅ phi_hat(θ, ϕ)]

@inline Emn(m, n) = (1im)^(n & 3) * √((2n + 1) * factorial(n - m) / (n * (n + 1) * factorial(n + m)))

@doc raw"""
Vector spherical wave function, electric (TM) modes.

```math
\begin{aligned}
\mathbf{N}_{m n}^{(i)}(\rho, \theta, \phi)=&\left\{\hat{\mathbf{e}}_{r} n(n+1) P_{n}^{m}(\cos \theta) \frac{z^{(i)}_{n}(\rho)}{\rho}\right.\\
&\left.+\left[\hat{\mathbf{e}}_{\theta} \tau_{m n}(\cos \theta)+\hat{\mathbf{e}}_{\phi} \mathrm{i} \pi_{m n}(\cos \theta)\right] \frac{(\rho\cdot z_n^{(i)}(\rho))^{\prime}}{\rho}\right\} \exp (\mathrm{i} m \phi)
\end{aligned}
```

where ``i=1,2,3,4`` relates to incident, interior, outgoing and ingoing VSWF.
"""
function vswf_electric(T::DataType, n::Integer, m::Integer, mode::VSWFMode, r::Number, θ::Number, ϕ::Number, k::Number)
    r = T(r)
    θ = T(θ)
    ϕ = T(ϕ)
    k = T(k)
    z = VSWF_KERNEL[Int(mode)]
    z_deriv = VSWF_KERNEL_DERIV[Int(mode)]
    kr = k * r
    expϕ = exp(1im * m * ϕ)
    H = z(T, n, kr)
    H_deriv = z_deriv(T, n, kr)
    Pnm = associated_legendre(T, n, m, cos(θ))

    factor = (H + kr * H_deriv) * expϕ / kr
    r_comp = n * (n + 1) * Pnm * H / kr * expϕ
    θ_comp = factor * tau_func(T, n, m, θ)
    ϕ_comp = 1im * factor * pi_func(T, n, m, θ)

    return StaticArrays.@SVector [r_comp, θ_comp, ϕ_comp]
end

@inline function vswf_electric(n::Integer, m::Integer, mode::VSWFMode, r::Number, θ::Number, ϕ::Number, k::Number)
    return vswf_electric(Float64, n, m, mode, r, θ, ϕ, k)
end

@doc raw"""
Vector spherical wave function, magnetic (TE) modes.

```math
\mathbf{M}_{m n}^{(i)}(\rho, \theta, \phi)=\left[\hat{\mathbf{e}}_{\theta} \mathrm{i} \pi_{m n}(\cos \theta)-\hat{\mathbf{e}}_{\phi} \tau_{m n}(\cos \theta)\right] z_{n}^{(i)}(\rho) \exp (\mathrm{i} m \phi)
```

where ``i=1,2,3,4`` relates to incident, interior, outgoing and ingoing VSWF.
"""
function vswf_magnetic(T::DataType, n::Integer, m::Integer, mode::VSWFMode, r::Number, θ::Number, ϕ::Number, k::Number)
    r = T(r)
    θ = T(θ)
    ϕ = T(ϕ)
    k = T(k)
    z = VSWF_KERNEL[Int(mode)]
    kr = k * r
    expϕ = exp(1im * m * ϕ)
    H = z(T, n, kr)

    factor = H * expϕ
    θ_comp = 1im * factor * pi_func(T, n, m, θ)
    ϕ_comp = -factor * tau_func(T, n, m, θ)

    return StaticArrays.@SVector [0, θ_comp, ϕ_comp]
end

@inline function vswf_magnetic(n::Integer, m::Integer, mode::VSWFMode, r::Number, θ::Number, ϕ::Number, k::Number)
    return vswf_magnetic(Float64, n, m, mode, r, θ, ϕ, k)
end

# TODO: Implement it
"""
Expand the electric field.
"""
function expand_E_cluster(T::DataType, mode::VSWFMode, k::Number)
    k = T(k)

    return
end

Caching.@cache function vswf_cache(T::DataType, n::Integer, m::Integer, ν::Integer, μ::Integer)
    CT = complex(T)
    factor = (-1)^(m & 1) * √((2ν + 1) * (2n + 1) * factorial(T, ν - μ) * factorial(T, n - m))
    qA = minimum((n, ν, (n + ν - abs(m + μ)) ÷ 2))
    qB = minimum((n, ν, (n + ν + 1 - abs(m + μ)) ÷ 2))
    A = CT <: Arblib.AcbLike ? Arblib.AcbRefVector(qA + 1) : zeros(CT, qA + 1)
    B = CT <: Arblib.AcbLike ? Arblib.AcbRefVector(qB + 1) : zeros(CT, qB + 1)

    for q in 0:(qA - 1)
        p = n + ν - 2q
        aq = gaunt_a(T, m, n, μ, ν, p)
        A[q + 1] = aq * (1im)^(p & 3) * T(n * (n + 1) + ν * (ν + 1) - p * (p + 1))
    end

    B[1] = 0
    for q in 1:(qB - 1)
        p = n + ν - 2q
        bq = gaunt_b(m, n, μ, ν, p)
        B[q + 1] = bq * (1im)^((p + 1) & 3) * √T(((p + 1)^2 - (n - ν)^2) * ((n + ν + 1)^2 - (p + 1)^2))
    end

    return VSWFCache(factor, A, B)
end

"""
Precompute required VSWF coefficients.
"""
function init_vswf_cache(T::DataType, lmax::Integer)
    indices = NTuple{4,Int}[(n, m, ν, μ) for n in 1:lmax for m in (-n):n for ν in 1:n for μ in (-ν):ν]

    Threads.@threads for (n, m, ν, μ) in indices
        vswf_cache(T, n, m, ν, μ)
    end
end

@inline init_vswf_cache(lmax::Integer) = init_vswf_cache(Float64, lmax)

function vsh_translation_insert_pair()
end
