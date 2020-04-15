p = pwd()
cd(joinpath(@__DIR__()))
run(`make`)
cd(p)
