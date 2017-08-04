p = pwd()
cd(Pkg.dir("FastPolynomialRoots/deps/"))
run(`make`)
cd(p)
