# FastPolynomialRoots.jl - Fast and backward stable computation of roots of polynomials
[![CI](https://github.com/andreasnoack/FastPolynomialRoots.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/andreasnoack/FastPolynomialRoots.jl/actions/workflows/CI.yml)
[![Coverage Status](https://coveralls.io/repos/github/andreasnoack/FastPolynomialRoots.jl/badge.svg?branch=master)](https://coveralls.io/github/andreasnoack/FastPolynomialRoots.jl?branch=master)

This package is a Julia wrapper of the Fortran programs accompanying [Fast and Backward Stable Computation of Roots of Polynomials](http://epubs.siam.org/doi/abs/10.1137/140983434) by Jared L. Aurentz, Thomas Mach, Raf Vandebril and David S. Watkins.

## Usage

The package provides the unexported function `FastPolynomialRoots.rootsFastPolynomialRoots(p::Vector{<:Union{Float64,Complex{Float64}}})`
which computes the roots of the polynomial `p[1] + p[2]*x + p[3]*x^2 + ... + p[k]*x^(k-1)`. The package also overwrites the `roots(::Polynomial)` methods in the `Polynomials` package for `Float64` and `Complex{Float64}` elements with the fast versions provided by this package. See the examples below.

## Example 1: Speed up `roots`
```julia
julia> using Polynomials, BenchmarkTools

julia> @btime roots(p) setup=(p = Polynomial(randn(500)));
  223.135 ms (23 allocations: 3.97 MiB)

julia> using FastPolynomialRoots

julia> @btime roots(p) setup=(p = Polynomial(randn(500)));
  30.786 ms (7 allocations: 26.41 KiB)
```

## Example 2: Roots of a polynomial of degree 10,000
A computation of this size would not be feasible on a desktop with the traditional method
but can be handled by FastPolynomialRoots.
```julia
julia> using Polynomials, BenchmarkTools, FastPolynomialRoots

julia> n = 10000;

julia> r = @btime roots(p) setup=(p = Polynomial(randn(n + 1)));
  10.290 s (13 allocations: 508.38 KiB)

julia> sum(isreal, r)
7

julia> 2/π*log(n) + 0.6257358072 + 2/(n*π) # Edelman and Kostlan
6.489284260212659
```
