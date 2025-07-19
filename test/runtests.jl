using Test, FastPolynomialRoots, Polynomials, LinearAlgebra

@testset "Standard normal coefficients" begin
    p = Polynomial(randn(50))
    @test sort(abs.(FastPolynomialRoots.roots(p))) ≈ sort(abs.(roots(p)))
end

@testset "Standard normal complex coefficients" begin
    p = Polynomial(complex.(randn(50), randn(50)))
    @test sort(abs.(FastPolynomialRoots.roots(p))) ≈ sort(abs.(roots(p)))
end

@testset "Integer coefficients (Issue 19)" begin
    p = Polynomial([1, 10, 100, 1000])
    @test sort(abs.(FastPolynomialRoots.roots(p))) ≈ sort(abs.(roots(p)))
end

@testset "Large polynomial" begin
    p = Polynomial(randn(5000))
    @time FastPolynomialRoots.roots(p)

    @info "Possible to calculate roots of large polynomial"
    # @show λs = 1:100.0
    λs = sort(randn(100), rev=true)
    p = fromroots(λs)
    @info "But polynomial root finding is ill conditioned"
    @test sum(abs2, FastPolynomialRoots.roots(p) - λs) < 1000
end
