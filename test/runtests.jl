using AMVW, Base.Test, Polynomials

@testset "Standard normal coefficients" begin
    p = Poly(randn(50))
    @test sort(abs.(roots(p))) ≈ sort(abs.(AMVW.rootsAMVW(p)))
end
@testset "Standard normal complex coefficient" begin
    p = Poly(complex.(randn(50), randn(50)))
    @test sort(abs.(roots(p))) ≈ sort(abs.(AMVW.rootsAMVW(p)))
end

@testset "Possible to calculate roots of large polynomial" begin
    p = Poly(randn(5000))
    @time AMVW.rootsAMVW(p)
end

@testset "But polynomial root finding is ill conditioned" begin
    rts = 1:100.0
    p = mapreduce(t -> Poly([1, -t]), *, Poly(Float64[1]), rts)
    @show norm(sort(abs.(AMVW.rootsAMVW(p))) - rts)
end
