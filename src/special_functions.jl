"""
Calculate factorials using the Gamma function.
"""
function factorial(T::DataType, n)
    if T <: Arblib.ArbLike
        return Arblib.gamma!(Arblib.Arb(), Arblib.Arb(n + 1))
    else
        ret = SpecialFunctions.gamma(n + 1)
        return isinf(ret) ? SpecialFunctions.gamma(big(n + 1)) : ret
    end
end

factorial(n) = factorial(Float64, n)

"""
Spherical Bessel function of the first kind.
"""
function spherical_jn(T::DataType, n::Integer, z)
    if T <: Arblib.ArbLike
        return √(Arblib.Arb(π) / (2z)) * Arblib.hypgeom_bessel_j!(Arblib.Arb(), Arblib.Arb(n + 1 // 2), Arblib.Arb(z))
    elseif T <: Arblib.AcbLike
        return √(Arblib.Arb(π) / (2z)) * Arblib.hypgeom_bessel_j!(Arblib.Acb(), Arblib.Acb(n + 1 // 2), Arblib.Acb(z))
    else
        return SpecialFunctions.sphericalbesselj(n, z)
    end
end

spherical_jn(n::Integer, z) = spherical_jn(Float64, n, z)

"""
First-order derivative of spherical Bessel function of the first kind.
"""
spherical_jn_deriv(T::DataType, n::Integer, z) =
    spherical_jn(T, n - 1, z) - (n + 1) / z * spherical_jn(T, n, z)

spherical_jn_deriv(n::Integer, z) = spherical_jn_deriv(Float64, n, z)

"""
Spherical Bessel function of the second kind.
"""
function spherical_yn(T::DataType, n::Integer, z)
    if T <: Arblib.ArbLike
        return √(Arblib.Arb(π) / (2z)) * Arblib.hypgeom_bessel_y!(Arblib.Arb(), Arblib.Arb(n + 1 // 2), Arblib.Arb(z))
    elseif T <: Arblib.AcbLike
        return √(Arblib.Arb(π) / (2z)) * Arblib.hypgeom_bessel_y!(Arblib.Acb(), Arblib.Acb(n + 1 // 2), Arblib.Acb(z))
    else
        return SpecialFunctions.sphericalbessely(n, z)
    end
end

spherical_yn(n::Integer, z) = spherical_yn(Float64, n, z)

"""
First-order derivative of spherical Bessel function of the second kind.
"""
spherical_yn_deriv(T::DataType, n::Integer, z) =
    spherical_yn(T, n - 1, z) - (n + 1) / z * spherical_yn(T, n, z)

spherical_yn_deriv(n::Integer, z) = spherical_yn_deriv(Float64, n, z)

"""
Spherical Hankel function of the first kind.
"""
spherical_hn1(T::DataType, n::Integer, z) =
    spherical_jn(T, n, z) + 1im * spherical_yn(T, n, z)

spherical_hn1(n::Integer, z) = spherical_hn1(Float64, n, z)

"""
First-order derivative of spherical Hankel function of the first kind.
"""
spherical_hn1_deriv(T::DataType, n::Integer, z) =
    spherical_jn_deriv(T, n, z) + 1im * spherical_yn_deriv(T, n, z)

spherical_hn1_deriv(n::Integer, z) = spherical_hn1_deriv(Float64, n, z)

"""
Spherical Hankel function of the second kind.
"""
spherical_hn2(T::DataType, n::Integer, z) =
    spherical_jn(T, n, z) - 1im * spherical_yn(T, n, z)

spherical_hn2(n::Integer, z) = spherical_hn2(Float64, n, z)

"""
First-order derivative of spherical Hankel function of the second kind.
"""
spherical_hn2_deriv(T::DataType, n::Integer, z) =
    spherical_jn_deriv(T, n, z) - 1im * spherical_yn_deriv(T, n, z)

spherical_hn2_deriv(n::Integer, z) = spherical_hn2_deriv(Float64, n, z)

"""
Riccati Bessel function of the first kind.
"""
ricatti_jn(T::DataType, n::Integer, z) =
    z * spherical_jn(T, n, z)

ricatti_jn(n::Integer, z) = ricatti_jn(Float64, n, z)

"""
First-order derivative of Riccati Bessel function of the first kind.
"""
function ricatti_jn_deriv(T::DataType, n::Integer, z)
    jn = spherical_jn(T, n, z)
    jnd = spherical_jn_deriv(T, n, z)
    return z * jnd + jn
end

ricatti_jn_deriv(n::Integer, z) = ricatti_jn_deriv(Float64, n, z)

"""
Riccati Bessel function of the second kind.

> Note that in `miepy`, the author used `-z⋅y(z)` instead of `z⋅y(z)`
"""
ricatti_yn(T::DataType, n::Integer, z) =
    z * spherical_yn(T, n, z)

ricatti_yn(n::Integer, z) = ricatti_yn(Float64, n, z)

"""
First-order derivative of Riccati Bessel function of the second kind.
"""
function ricatti_yn_deriv(T::DataType, n::Integer, z)
    yn = spherical_yn(T, n, z)
    ynd = spherical_yn_deriv(T, n, z)
    return z * ynd + yn
end

ricatti_yn_deriv(n::Integer, z) = ricatti_yn_deriv(Float64, n, z)

"""
Riccati Hankel function of the first kind.
"""
ricatti_hn1(T::DataType, n::Integer, z) =
    ricatti_jn(T, n, z) + 1im * ricatti_yn(T, n, z)

ricatti_hn1(n::Integer, z) = ricatti_hn1(Float64, n, z)

"""
First-order derivative of Riccati Hankel function of the first kind.
"""
ricatti_hn1_deriv(T::DataType, n::Integer, z) =
    ricatti_jn_deriv(T, n, z) + 1im * ricatti_yn_deriv(T, n, z)

ricatti_hn1_deriv(n::Integer, z) = ricatti_hn1_deriv(Float64, n, z)

"""
Riccati Hankel function of the second kind.
"""
ricatti_hn2(T::DataType, n::Integer, z) =
    ricatti_jn(T, n, z) - 1im * ricatti_yn(T, n, z)

ricatti_hn2(n::Integer, z) = ricatti_hn2(Float64, n, z)

"""
First-order derivative of Riccati Hankel function of the second kind.
"""
ricatti_hn2_deriv(T::DataType, n::Integer, z) =
    ricatti_jn_deriv(T, n, z) - 1im * ricatti_yn_deriv(T, n, z)

ricatti_hn2_deriv(n::Integer, z) = ricatti_hn2_deriv(Float64, n, z)

@doc raw"""
Associated Legendre function without the Condon-Shotley phase.

`Arblib`'s definition includes the Condon-Shotley phase, so we need to multiply the results by ``(-1)^m``.
"""
function associated_legendre(T::DataType, n::Integer, m::Integer, z)
    if T <: Arblib.ArbLike
        return (-1)^(m & 1) * Arblib.hypgeom_legendre_p!(Arblib.Arb(), Arblib.Arb(n), Arblib.Arb(m), Arblib.Arb(z), 0)
    elseif T <: Arblib.AcbLike
        return (-1)^(m & 1) * Arblib.hypgeom_legendre_p!(Arblib.Acb(), Arblib.Acb(n), Arblib.Acb(m), Arblib.Acb(z), 0)
    else
        index = GSL.sf_legendre_array_index(n, abs(m))
        val = GSL.sf_legendre_array(GSL.GSL_SF_LEGENDRE_NONE, n, z)[index+1]
        if m < 0
            val *= (-1)^(m & 1) * factorial(T, n + m) / factorial(T, n - m)
        end
        return val
    end
end

associated_legendre(n::Integer, m::Integer, z) = associated_legendre(Float64, n, m, z)

@doc raw"""
Calculate the associated Legendre function for all ``0\leq n\leq n_{\max}`` and ``-n\leq m\leq n``, and return the results as a vector.

`Arblib`'s Legendre function definition includes the Condon-Shotley phase, so we need to multiply the results by ``(-1)^m``.
"""
function associated_legendre_array(T::DataType, nmax::Integer, z)
    terms = (nmax + 1)^2

    if T <: Arblib.ArbLike
        ret = Arblib.ArbRefVector(terms)
        z_arb = Arblib.Arb(z)
        for n in 0:nmax
            n_arb = Arblib.Arb(n)
            for m in -n:n
                index = n * (n + 1) + m + 1
                Arblib.hypgeom_legendre_p!(ret[index], n_arb, Arblib.Arb(m), z_arb, 0)
                ret[index] *= (-1)^(m & 1)
            end
        end
    elseif T <: Arblib.AcbLike
        ret = Arblib.AcbRefVector(terms)
        z_acb = Arblib.Acb(z)
        for n in 0:nmax
            n_acb = Arblib.Acb(n)
            for m in -n:n
                index = n * (n + 1) + m + 1
                Arblib.hypgeom_legendre_p!(ret[index], n_acb, Arblib.Acb(m), z_acb, 0)
                ret[index] *= (-1)^(m & 1)
            end
        end
    else
        arr = GSL.sf_legendre_array(GSL.GSL_SF_LEGENDRE_NONE, nmax, z)
        ret = zeros(T, terms)
        for n in 0:nmax
            for m in 0:n
                index = GSL.sf_legendre_array_index(n, m)
                val = arr[index+1]
                ret[n*(n+1)+m+1] = val
                if m > 0
                    ret[n*(n+1)-m+1] =
                        val * (-1)^(m & 1) * factorial(T, n - m) / factorial(T, n + m)
                end
            end
        end
    end

    return ret
end

associated_legendre_array(nmax::Integer, z) = associated_legendre_array(Float64, nmax, z)

@doc raw"""
Calculate ``\pi_n^m(\theta)`` defined as

```math
\pi_n^m(\theta)=\frac{m}{\sin\theta}\cdot P_n^m(\cos\theta)
```
"""
function pi_func(n::Integer, m::Integer, θ::Number)
    T = typeof(θ)
    if iszero(θ)
        if m == 1
            return T(n * (n + 1)) / 2
        elseif m == -1
            return T(1 // 2)
        else
            return zero(T)
        end
    elseif isapprox(θ, π)
        if m == 1
            return T(n * (n + 1)) * (-1)^((n + 1) & 1) / 2
        elseif m == -1
            return T((-1)^((n + 1) & 1) // 2)
        else
            return zero(T)
        end
    else
        z = cos(θ)
        return m / sin(θ) * associated_legendre(T, n, m, z)
    end
end
