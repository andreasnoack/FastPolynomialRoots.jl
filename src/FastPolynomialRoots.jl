module FastPolynomialRoots

using Requires

function __init__()

    @require Polynomials="f27b6e38-b328-58d1-80ce-0feddd5e7a45" begin
        Polynomials.roots(p::Union{Polynomials.Poly{Float64},Polynomials.Poly{Complex{Float64}}}) = rootsFastPolynomialRoots(p.a)
        Polynomials.roots(p::Polynomials.Poly{<:Integer}) = rootsFastPolynomialRoots(convert(Polynomials.Poly{Float64}, p))
    end
end

const dpath = joinpath(@__DIR__(), "..", "deps", "libamvwdouble")
const spath = joinpath(@__DIR__(), "..", "deps", "libamvwsingle")

function rootsFastPolynomialRoots(a::Vector{Float64})

    pl    = reverse!(a[1:end - 1] ./ a[end])
    np    = length(pl)
    reigs = similar(pl)
    ieigs = similar(pl)
    its   = Vector{Int32}(undef, np)
    flag  = Int32[0]

    ccall((:damvw_, dpath), Cvoid,
        (Ref{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}),
        np, pl, reigs, ieigs, its, flag)

    if flag[1] != 0
        error("error code: $(flag[1])")
    end
    return complex.(reigs, ieigs)
end

function rootsFastPolynomialRoots(a::Vector{Complex{Float64}})

    pl    = reverse!(a[1:end - 1] ./ a[end])
    plr   = real(pl)
    pli   = imag(pl)
    np    = length(pl)
    reigs = similar(plr)
    ieigs = similar(plr)
    its   = Vector{Int32}(undef, np)
    flag  = Int32[0]

    ccall((:zamvw_, spath), Cvoid,
        (Ref{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}),
        np, plr, pli, reigs, ieigs, its, flag)

    if flag[1] != 0
        error("error code: $(flag[1])")
    end
    return complex.(reigs, ieigs)
end

end # module