#
# Compute (& store) landen sequence, and other fun stuff
#

# TODO
# check for purely imaginary k on negative axis-> Rz < eps(Iz)

module Landen
# import ..Elliptic: K
include("misc.jl")

export LandenSeq, NonConvergedLandenSeq

# landen can be used on complex moduli, and k can be an Int
# √(eps) because 1 - k² = 1 when k² ≤ eps(k), so m is indistinguishable from 0
# need k ≤ _kmin for K′ approximation convergence
# also, see _sncndn_seq_end for √eps decision
_default_tol(k) = min(√(_one_eps(k)), _kmin)

const _kmid = 1/sqrt(2.0)
# from Lecture Notes on Elliptic Filter Design, Sophocles J. Orfanidis, 2006
const _kmin = 1e-6
const _kmax = 0.9999999999995 # _k′(_kmin)

abstract type LandenDirection end
struct DescendingLanden <: LandenDirection end
struct AscendingLanden <: LandenDirection end

swap(::Type{AscendingLanden}, a, b) = (b, a)
swap(::Type{DescendingLanden}, a, b) = (a, b)

################################################################################
# Landen Sequence
################################################################################

struct LandenSeq{N, T<:RealOrComplexFloat, D<:LandenDirection}
    # basically an SArray{(2, N), T, 2*N}
    ks::NTuple{N, T}
    k′s::NTuple{N, T}

    # non-converged case
    LandenSeq(U::Type, D::Type{<:LandenDirection}) = new{0, U, D}((), ())

    function LandenSeq(k::T, k′::Maybe{T}, N::Integer, ktol::AbstractFloat,
            D::Type{<:LandenDirection}) where {T<:RealOrComplexFloat}
        ka = abs(k)

        # all equations/integrals are in terms of k², so pick posiitve root
        (N < 0) || (0 ≤ ka ≤ 1) || !(isfinite(k)) || return LandenSeq(T, D)
        # safe to do √
        isa(k′, Nothing) && (k′ = _k′(k))

        k, k′ = swap(D, _abs_real(k), _abs_real(k′))

        # if already within range, no iterations needed
        (ka ≤ ktol) && (N = 0)

        ks = zeros(T, N+1)
        k′s = zeros(T, N+1)

        ks[1], k′s[1] = k, k′

        n = 0
        for i in 2:(N+1)
            k, k′ = _landen_kernel_stable(k, k′)
            # k, k′ = (k, k′) ./ hypot(k, k′) # dirty hack
            @inbounds ks[i], k′s[i] = k, k′

            if abs(k) ≤ ktol
                n = i
                break
            end
        end

        # N == 0 -> already converged, N > 1 && n == 0 -> did not converge
        n += (N == 0)
        ks, k′s = swap(D, ks, k′s)
        return new{n, T, D}((ks[1:n]..., ), (k′s[1:n]..., ))
    end
end

const NonConvergedLandenSeq = LandenSeq{0}
# descending is arbitrary
"""
    LandenSeq(U::Type, D::Type{<:LandenDirection}=DescendingLanden)

Return a `LandenSeq{0, U, D}() <: NonConvergedLandenSeq`.
"""
LandenSeq(U::Type) = LandenSeq(0, U, DescendingLanden)
LandenSeq(k; args...) = LandenSeq(Base.Bottom)

"""
    LandenSeq(k, [k′]; N::Int=10, ktol=√eps(k), descending=(k ≤ 1/√2))

Return the elliptic moduli, `{(kₙ, k′ₙ)}`, in a Landen sequence starting with `k₀ = k`

Compute and store `LandenSeq` containing `{(kₙ, k′ₙ)}` of a descending (ascending) Landen
sequence starting from `k₀ = k` until kₙ converges to 0 (1) via criteria
`abs(kₙ) ≤ ktol` (`abs(k′ₙ) ≤ ktol` or `abs(kₙ) ≥ 1 - ktol`), with a maximum of
`N` iterations performed. `K′ₙ/Kₙ` doubles (halves) per iteration.

Return a `landen <: LandenSeq{0} == NonConvergedLandenSeq` if `k,k′ ∉ ℂ`,
`(abs(k) > 1)`, or the sequence did not converge in `N` iterations.
"""
function LandenSeq(k::Number, k′::Maybe{Number}=nothing; N=10,
        ktol=_default_tol(k), descending=(abs(k) ≤ _kmid))
    if isa(k′, Nothing)
        k = _promote_float(k)
    else
        k, k′ = _promote_float(k, k′)
    end

    D = (descending) ? DescendingLanden : AscendingLanden

    return LandenSeq(k, k′, N, ktol, D)
end

Base.length(::LandenSeq{N}) where N = N
Base.firstindex(::LandenSeq) = 1
Base.lastindex(landen::LandenSeq) = length(landen)
Base.getindex(landen::LandenSeq, i::Int) = (landen.ks[i], landen.k′s[i])
Base.getindex(landen::LandenSeq, I) = [landen[i] for i in I]

ktype(::LandenSeq{N, T}) where {N,T} = T
Base.eltype(landen::LandenSeq)= NTuple{2, ktype(landen)}
Base.iterate(landen::LandenSeq, state=1) = state ≤ length(landen) ?
        (landen[state], state+1) : nothing

direction(::LandenSeq{N,T,D}) where {N,T,D} = D

# descending landen iteration k → 0
# ascending kernel is this with arguments swappend, and output swapped again
# not accuarte for k ≈ 1
# assumes k ∉ (-∞, 0)
@inline function _landen_kernel(k, k′)
    if abs2(k)^6 ≤ _one_eps(k)
        m = k^2
        # taylor series expansions around k = 0
        return ( ldexp(@horner(m, 0, 256, 128, 80, 56, 42, 33), -10),
            ldexp(@horner(m, 65536, 0, -2048, -2048, -1824, -1600, -1409), -16))
    else
        return ((k/(1+k′))^2, 2*sqrt(k′)/(1+k′))
    end
end

_landen_kernel_stable(k::Real, k′::Real) = _landen_kernel(k, k′)
@inline function _landen_kernel_stable(k::Complex, k′::Complex)
    k, k′ = _landen_kernel(k, k′)
    return (_abs_real(k), _abs_real(k′))
end

################################################################################
# K: quarter period
################################################################################

"""
    K(landen::LandenSeq)

Computes `K(k) * K′(k)` for a Landen sequence.
Returns `NaN` for `NonConvergedLandenSeq`
"""
function K(landen::LandenSeq)
    ks = landen.ks
    k′s = landen.k′s

    D = direction(landen)
    ks, k′s = swap(D, ks, k′s)

    KK, KK′ = _K(ks, k′s)
    return swap(D, KK, KK′)
end

K(::NonConvergedLandenSeq) = NaN

# compute K & K′ descending landen
# not accurate for k₀ ≥ kmax
# Lecture Notes on Elliptic Filter Design, Sophocles J. Orfanidis, 2006
# https://dlmf.nist.gov/19.8.E12
@inline function _K(ks, k′s)
    k = first(ks)
    k′ = first(k′s)
    kₙ = last(ks)
    k′ₙ = last(k′s)

    K = oftype(k, pi)/2
    for n in 2:length(ks)
        @inbounds k = ks[n]
        K *= 1 + k
    end

    K′ = _largeK(k′ₙ, kₙ)
    denom = oftype(k, 1)
    for n in 1:length(k′s)-1
        @inbounds k′ = k′s[n]
        denom *= 1 + k′
    end

    return (K, K′/denom)
end

# only applicable if k ≥ _kmax ≈ √(1 - 1e-12)
@inline function _largeK(::T, k′::T) where T<:RealOrComplexFloat
    K = -log(ldexp(k′, -2))
    isinf(K) && return K

    # pretty small numbers, so may not contribute much
    K += (K - 1)/2 * k′^2

    return K
end

end # module
