@testset "Special Functions" begin
    @testset "Wigner d-function" begin
        smax = 5
        angles = rand(10) .* π
        @show angles

        @testset "m = $m, n = $n" for m in -2:2, n in -2:2
            for θ in angles
                d, d_deriv = GeneralizedMultiparticleMie.wigner_d_recursion(Float64, m, n, smax, θ; deriv=true)

                smin = max(abs(m), abs(n))
                for j in smin:smax
                    dd = WignerD.wignerdjmn(j, m, n, θ)
                    d1 = GeneralizedMultiparticleMie.wigner_d(Float64, m, n, j, θ)
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
end
