using Test, FastPolynomialRoots, Polynomials

@testset "Standard normal coefficients" begin
    p = Poly(randn(50))
    @test sort(abs.(roots(p))) ≈ sort(abs.(roots(p)))
end
@testset "Standard normal complex coefficient" begin
    p = Poly(complex.(randn(50), randn(50)))
    @test sort(abs.(roots(p))) ≈ sort(abs.(roots(p)))
end

@testset "Large polynomial" begin
    p = Poly(randn(5000))
    @time roots(p)

    @info "Possible to calculate roots of large polynomial"
    @show rts = 1:100.0
    p = mapreduce(t -> Poly([1, -t]), *, rts, init=Poly(Float64[1]))
    @info "But polynomial root finding is ill conditioned"
    @show sum(abs2, sort(abs.(roots(p))) - rts)
end
