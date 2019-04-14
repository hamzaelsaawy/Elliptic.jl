#
# Compute (& store) landen sequence, and other fun stuff
#
module Landen

import Base.Math.@horner
import ..Elliptic: K

export LandenSeq, NonConvergedLandenSeq

################################################################################
# miscellaneous
################################################################################

const Maybe{T} = Union{Nothing, T}
const RealOrComplexFloat{T<:AbstractFloat} = Union{T, Complex{T}}

# do everything as a float
_one_eps(x) = eps(float(real(typeof(x))))

_float_typeof(x) = float(typeof(x))
_promote_float_typeof(x) = _float_typeof(x)
@inline _promote_float_typeof(x, xs...) =
        promote_type(_promote_float_typeof(x), _promote_float_typeof(xs...))
_promote_float(ts...) = convert.(_promote_float_typeof(ts...), ts)

# landen can be used on complex moduli, and k can be an Int
# √(eps) because 1 - k² = 1 when k² ≤ eps(k), so m is indistinguishable from 0
# need k ≤ _kmin for K′ approximation convergence
# also, see _sncndn_seq_end for √eps decision
_default_tol(k) = min(√(_one_eps(k)), _kmin)

# rotate numbers by pi to have phase in (-pi/2, pi/2)
_positive_root(x::Real) = abs(x)
_positive_root(x::Complex) = ifelse(real(x) < 0, -x, x)

const _kmid = 1/sqrt(2.0)
# from Lecture Notes on Elliptic Filter Design, Sophocles J. Orfanidis, 2006
const _kmin = 1e-6
const _kmax = 0.9999999999995 # _k′(_kmin)

 _k′(k) = sqrt((1-k)*(1+k))

################################################################################
# Landen Sequence
################################################################################

# TODO
# complex k -> check which root is used, push k to have |ph| < pi/2, Rz > 0
#  pq(z, k) = pq(z, -k) -> http://dlmf.nist.gov/22.17.E1
# check for purely imaginary k -> Rz < eps(Iz)

abstract type LandenDirection end
struct DescendingLanden <: LandenDirection end
struct AscendingLanden <: LandenDirection end

swap(::Type{AscendingLanden}, a, b) = (b, a)
swap(::Type{DescendingLanden}, a, b) = (a, b)

struct LandenSeq{N, T<:RealOrComplexFloat, D<:LandenDirection}
    # basically an SArray{(2, N), T, 2*N}
    ks::NTuple{N, T}
    k′s::NTuple{N, T}

    # non-converged case
    LandenSeq(U::Type, D::Type{<:LandenDirection}) = new{0, U, D}((), ())

    function LandenSeq(k::T, k′::T, N::Integer, ktol::AbstractFloat,
            D::Type{<:LandenDirection}) where {T<:RealOrComplexFloat}
        ka = abs(k)
        # all equations/integrals are in terms of k², so can pick posiitve root
        # (real(k) < 0) && (k *= -1)
        (N < 0) || (0 ≤ ka ≤ 1) || !(isfinite(k)) || return LandenSeq(T, D)

        k, k′ = swap(D, _positive_root(k), _positive_root(k′))

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
    isa(k′, Nothing) && (k′ = _k′(k))
    k, k′ = _promote_float(k, k′)

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
        return ( @horner(m, 0, 256, 128, 80, 56, 42, 33) / 1024,
            @horner(m, 65536, 0, -2048, -2048, -1824, -1600, -1409) / 65536)
    else
        return ((k/(1+k′))^2, 2*sqrt(k′)/(1+k′))
    end
end

_landen_kernel_stable(k::Real, k′::Real) = _landen_kernel(k, k′)
@inline function _landen_kernel_stable(k::Complex, k′::Complex)
    k, k′ = _landen_kernel(k, k′)
    return (_positive_root(k), _positive_root(k′))
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
    K = log(T(4)) - log(k′)
    isinf(K) && return K

    # pretty small numbers, so may not contribute much
    K += (K - 1)/2 * k′^2

    return K
end

################################################################################
# amplitude
################################################################################
"""
    gd(x)

Gudermannian funciton, `gd(x)`, equal to the integral of `sech` from `0` to `x`
"""
gd(x::Number) = asin(tanh(x))

"""
    agd(x)

Inverse Gudermannian funciton, `gd(x)`, equal to the integral of `sec` from `0` to `x`
"""
agd(x) = asinh(tan(x))

# TODO, the rest of this
# am -> small u k, k′ approximations

################################################################################
# jacobi functions
################################################################################

# TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# pq when k² ≥ 0 (or ph(k^2) = π)
# check for pole around `n = jK′` for sncndn
# condition number for _sncndn_small_z

sncndn(_, ::NonConvergedLandenSeq) = NaN

function sncndn(z, landen::LandenSeq)
    T = _promote_float_typeof(z, ktype(landen))
    z = T(z)

    (abs2(z)^4 ≤ _one_eps(z)) &&
            return _sncndn_small_z(z, last(landen)...)

    return _sncndn_seq(z, landen)
end

_sncndn_seq(z, landen::LandenSeq{1}) = _sncndn_seq_end(z, landen)

function _sncndn_seq(z, landen::LandenSeq)
    N = length(landen)
    D = direction(landen)

    ks = landen.ks
    k′s = landen.k′s

    K, K′ = _K(landen)
    u = z / K
    sn, cn, dn = _sncndn_seq_end(u * π/2, landen)

    for i in N:-1:2
        @inbounds k = ks[i]
        @inbounds k′ = k′s[i]

        sn, cn, dn = _sncndn_seq_iter(sn, cn, dn, k, k′, D)
    end

    return (sn, cn, dn)
end

# https://dlmf.nist.gov/22.7

# given a descending sequence of moduli, go back up
@inline function _sncndn_seq_iter(sn, cn, dn, k, k′, ::Type{DescendingLanden})
    n = 1 + k * sn^2
    s = (1 + k) * sn / n
    c = cn * dn / n
    d = hypot(k′* s, c)
    # avoid cancelation
    # (1 - ksn2) / (1 + ksn2))
    # (dn^2  - 1 + k) / (1 + k - dn^2)

    return (s, c, d)
end

# given an ascending sequence of moduli, go back down
@inline function _sncndn_seq_iter(sn, cn, dn, k, k′, ::Type{AscendingLanden})
    m = k^2
    dn2 = dn^2
    mdn = m * dn
    s = (1 + k′) * sn * cn / dn
    # ascending should have k′ ≈ 1
    c = m * (1 - (1+k′) * sn^2) / (m * dn)
    # (1 + k′) * (dn2 - k′) / (m * dn)
    d = hypot(k′ * s, c)
    # (1 - k′) * (dn2 + k′) / (m * dn)

    return (s, c, d)
end

#
# jacobi approximations
#
# https://dlmf.nist.gov/22.10

# |z| ≈ 0, |z|⁸ ≤ eps(T)
# converge when |z| < min(K(k), K′(k)), but π/2 ≤ min(K, K′)
function _sncndn_small_z(z, k)
    m = k^2
    zz = z^2

    # sn = @horner(z, 0, 1, 0, -(1 + m)/6, 0, (1 + m*(14 + m))/120, 0,
    #     -(1 + m*(135 + m*(135 + m)))/5040)
    # cn = @horner(z, 1, 0, -0.5, 0, (1+4*m)/24, 0, -(1 + m*(44 + 16*m))/720)
    # dn = @horner(z, 1, 0, -m/2, 0, m*(4 + m)/24, 0, -m*(16 + m*(44 + m))/720)
    sn = z * @horner(zz, 1, -(1 + m)/6, (1 + m*(14 + m))/120,
        -(1 + m*(135 + m*(135 + m)))/5040)
    cn = @horner(zz, 1, -0.5, (1+4*m)/24, -(1 + m*(44 + 16*m))/720)
    dn = @horner(zz, 1, -m/2, m*(4 + m)/24, -m*(16 + m*(44 + m))/720)

    return (sn, cn, dn)
end

# k ≈ 0, |k|³ ≤ eps(T), i guess?, its O(k⁴), but sn/cn(z) can be ≤ k
# lets not use unless |k|² ≤ eps(T)
_sncndn_seq_end(z, landen::LandenSeq{N,T,AscendingLanden}) where {N,T} =
        _sncndn_hyper(z, last(landen.k′s))

_sncndn_seq_end(z, landen::LandenSeq{N,T,DescendingLanden}) where {N,T} =
        _sncndn_circular(z, last(landen.ks))

# keep these function names around for potential use elsewhere
function _sncndn_circular(z, k)
    m = k^2
    sz = sin(z)
    cz = cos(z)
    inner = m * (z - sz * cz) / 4

    sn = sz - inner * cz
    cn = cz + inner * sz
    dn = 1 - (k * sin(z))^2 / 2

    return (sn, cn, dn)
end

function _sncndn_hyper(z, k′)
    m = k′^2 / 4
    sz = sinh(z)
    cz = cos(z)
    scz = sech(z)
    tz = tanh(z)
    sc = sz * cz
    inner = m * (z - sz * cz) / 4

    sn = tz - m * (z - sc)*scz^2
    cn = scz + m * (z - sc) * tz * scz
    dn = scz + m * (z + sc) * tz * scz

    return (sn, cn, dn)
end

end # module
