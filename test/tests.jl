using AMVW, Base.Test	

# Standard normal coefficients
p = Poly(randn(50))
@test_approx_eq sort(abs(roots(p))) sort(abs(AMVW.rootsAMVW(p)))

# Standard normal complex coefficient
p = Poly(complex(randn(50), randn(50)))
@test_approx_eq sort(abs(roots(p))) sort(abs(AMVW.rootsAMVW(p)))

# Possible to calculate roots of large polynomial
p = Poly(randn(5000))
@time AMVW.rootsAMVW(p)

# But polynomial root finding is ill conditioned
rts = 1:100.0
p = mapreduce(t -> Poly([1, -t]), *, Poly(Float64[1]), rts)

println(norm(sort(abs(AMVW.rootsAMVW(p))) - rts))