using Polynomial

module AMVW

	using Polynomial: Poly
	
	const dpath = joinpath(Pkg.dir("AMVW"), "deps", "libamvwdouble")

	function rootsAMVW(p::Poly{Float64})

		pl = p.a[2:end] ./ p.a[1]
		np = length(pl)
		reigs = similar(pl)
		ieigs = similar(pl)
		its = Array(Int32, np)
		flag = Int32[0]
		
		ccall((:damvw_, dpath), Void,
			(Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}),
			&np, pl, reigs, ieigs, its, flag)

		if flag[1] != 0 error("error code: $(flag[1])") end
		return complex(reigs, ieigs)
	end

end