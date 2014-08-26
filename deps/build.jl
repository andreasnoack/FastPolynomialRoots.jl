p = pwd()
cd(Pkg.dir("AMVW/deps/"))
run(`make`)
cd(p)
