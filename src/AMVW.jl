module AMVW

using Requires

    const dpath = joinpath(Pkg.dir("AMVW"), "deps", "libamvwdouble")
    const spath = joinpath(Pkg.dir("AMVW"), "deps", "libamvwsingle")

    function rootsAMVW(a::Vector{Float64})

        pl    = reverse!(a[1:end - 1] ./ a[end])
        np    = length(pl)
        reigs = similar(pl)
        ieigs = similar(pl)
        its   = Vector{Int32}(np)
        flag  = Int32[0]

        ccall((:damvw_, dpath), Void,
            (Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}),
            &np, pl, reigs, ieigs, its, flag)

        if flag[1] != 0
            error("error code: $(flag[1])")
        end
        return complex.(reigs, ieigs)
    end

    function rootsAMVW(a::Vector{Complex{Float64}})

        pl    = reverse!(a[1:end - 1] ./ a[end])
        plr   = real(pl)
        pli   = imag(pl)
        np    = length(pl)
        reigs = similar(plr)
        ieigs = similar(plr)
        its   = Vector{Int32}(np)
        flag  = Int32[0]

        ccall((:zamvw_, spath), Void,
            (Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}),
            &np, plr, pli, reigs, ieigs, its, flag)

        if flag[1] != 0
            error("error code: $(flag[1])")
        end
        return complex.(reigs, ieigs)
    end


    @require Polynomials begin
        using Polynomials: Poly
        import Polynomials: roots
        roots(p::Union{Poly{Float64},Poly{Complex{Float64}}}) = rootsAMVW(p.a)
        roots{T<:Integer}(p::Poly{T}) = rootsAMVW(convert(Poly{Float64}, p))
    end

end