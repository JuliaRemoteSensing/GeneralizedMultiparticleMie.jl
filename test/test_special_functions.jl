@testset "Special Functions" begin
    @testset "Wigner d-function" begin
        smax = 5
        angles = rand(10) .* π

        @testset "m = $m, n = $n" for m in -2:2, n in -2:2
            for θ in angles
                d, d_deriv = GeneralizedMultiparticleMie.wigner_d_recursion(Float64, m, n, smax, θ; deriv=true)

                smin = max(abs(m), abs(n))
                for j in smin:smax
                    dd = WignerD.wignerdjmn(j, m, n, θ)
                    d1 = GeneralizedMultiparticleMie.wigner_d_naive(Float64, m, n, j, θ)
                    @test d[j] ≈ dd ≈ d1
                end

                if m == 0 && n == 0
                    _, dd_deriv = GSL.sf_legendre_deriv_array(0, smax, cos(θ))
                    for j in smin:smax
                        idx = GSL.sf_legendre_array_index(j, 0) + 1
                        @test d_deriv[j] ≈ dd_deriv[idx] * -sin(θ)
                    end
                end
            end
        end
    end

    @testset "Wigner D-function" begin
        α, β, γ = rand(3) .* π

        @testset "m = $m, m′ = $m′, n = $n" for m in -2:2, m′ in -2:2, n in max(abs(m), abs(m′)):2
            D = GeneralizedMultiparticleMie.wigner_D(Float64, m, m′, n, α, β, γ)
            DD = WignerD.wignerDjmn(n, m, m′, α, β, γ)
            @test D ≈ DD
        end
    end
end
