#
# Compute (& store) landen sequence, and other fun stuff
#
module Landen

import ..Elliptic: K

export LandenSeq, landenseq, NonConvergedLandenSeq

const Maybe{T} = Union{Nothing, T}
const RealOrComplex{T<:AbstractFloat} = Union{T, Complex{T}}

#
# Landen Sequence
#

abstract type LandenDirection end
struct DescendingLanded <:LandenDirection end
struct AscendingLanded <:LandenDirection end


#Returns a tuple of kₙs in descending (ascending) Landen sequence where K′ₙ/Kₙ doubles (halves)
"""
    LandenSeq(k, [k′]; N::Int=10, ktol=√eps(k), descending=(k ≤ 1/√2))

Return the elliptic moduli, `{(kₙ, k′ₙ)}`, in a Landen sequence starting with `k₀ = k`

Compute and store `LandenSeq` containing `{(kₙ, k′ₙ)}` of a descending (ascending) Landen
sequence starting from `k₀ = k` until kₙ converges to 0 (1) via criteria `kₙ ≤ ktol`
(`k′ₙ ≤ ktol` or `kₙ ≥ 1 - ktol`).
Returns a `landen <: LandenSeq{0} == NonConvergedLandenSeq` if `(k ∉ ℝ & k ∉ [0, 1])`,
`(k ∉ ℂ & abs(k) ∉ [0, 1])`, or the sequence did not converge in `N` iterations.
"""
struct LandenSeq{N, T<:RealOrComplex}
    # basically a SArray{(2, N), T, 2*N}
    ks::NTuple{N, T}
    k′s::NTuple{N, T}

    LandenSeq(U::Type) = new{0, U}((), ())
    function LandenSeq(k::Number, k′::Maybe{Number}=nothing; N=10, ktol=_default_tol(k),
            descending=(k ≤ 1/sqrt(2)))
        k = float(k)
        T = typeof(k)
        # abs for complex moduli
        ka = isreal(k) ? k : abs(k)
        (0 ≤ ka ≤ 1) || return LandenSeq(T)

        # cant compute before checking value
        isa(k′, Nothing) && (k′ = _k′(k))

        # check if already within range
        if (descending && (abs(k) ≤ ktol)) || (!descending && (abs(k′) ≤ ktol))
            return new{1, T}((k, ), (k′, ))
        end

        descending || ((k, k′) = (k′, k))

        ks = zeros(T, N+1)
        k′s = zeros(T, N+1)

        ks[1], k′s[1] = k, k′

        # landen_kernel = descending ? _desc_landen_kernel : _asc_landen_kernel
        # k_conv = descending ? ktol : (1 - ktol)

        n = 0
        @inbounds for i in 2:(N+1)
            k, k′ = _desc_landen_kernel(k, k′)
            # dirty hack
            # k, k′ = (k, k′) ./ hypot(k, k′)
            ks[i], k′s[i] = k, k′

            if abs(k) ≤ ktol # _converged(k, k′, ktol, descending)
                n = i
                break
            end
        end

        descending || ((ks, k′s) = (k′s, ks))
        return new{n, T}((ks[1:n]..., ), (k′s[1:n]..., ))
    end
end

const NonConvergedLandenSeq = LandenSeq{0}
# when `k` is not a number
LandenSeq(k; args...) = LandenSeq(Base.Bottom)

Base.length(::LandenSeq{N}) where N = N
Base.firstindex(::LandenSeq) = 1
Base.lastindex(landen::LandenSeq) = length(landen)
Base.getindex(landen::LandenSeq, i::Int) = (landen.ks[i], landen.k′s[i])
Base.getindex(landen::LandenSeq, I) = [landen[i] for i in I]

ktype(::LandenSeq{N, T}) where {N,T} = T
Base.eltype(landen::LandenSeq)= NTuple{2, ktype(landen)}
Base.iterate(landen::LandenSeq, state=1) = state ≤ length(landen) ?
        (landen[state], state+1) : nothing

# landen can be used on complex moduli, and k can be an Int
# √(eps) because 1 - k² = 1 when k² ≤ eps(k), so m is indistinguishable from 0
@inline  _default_tol(k) = √(eps(float(real(one(k)))))

const _kmid = 1/sqrt(2)

 # TODO improve accuraccy for k ≈ 1
 _k′(k) = sqrt((1-k)*(1+k))

## TODO eval stablity and error accumulation here
# esp descending


# descending landen iteration k → 0
# not accuarte for k ≈ 1
@inline function _desc_landen_kernel(k, k′=_k′(k))
    if k ≤ eps(typeof(k))^(1/5)
        kk = abs2(k)
        # taylor series expansions around k = 0
        return (1/4*kk*(1 + 1/2*kk*(1 + 5/8*kk*(1 + 7/2*kk))),
                (1 - kk^2/32*(1 - kk*(1 - 57/64*kk))))
    else
        return (abs2(k/(1+k′)), 2*sqrt(k′)/(1+k′))
    end
end

# ascending landen iteration k → 1
#  ignored second argument to match `desc_landen_kernel` call
@inline _asc_landen_kernel(k, k′) = reverse(_desc_landen_kernel(k′, k))

"""
    landenseq(k::Number; N::Int=10, ktol=eps(k), descending=true)

Return a tuple of `{kₙ}` in a Landen sequence starting from `k₀ = k`

Return `()::Tuple{}` for invalid inputs or if the sequence did not converge.
See [`LandenSeq`](@ref) for options or more details.
"""
landenseq(args...; kargs...) = LandenSeq(args...; kargs...).ks

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

    descending = last(ks) < last(k′s)
    descending || ((ks, k′s) = (k′s, ks))

    KK, KK′ = _K(ks, k′s)

    return descending ? (KK, KK′) : (KK′, KK)
end

K(::NonConvergedLandenSeq) = NaN

# compute K & K′ descending landen
# not accurate for k₀ ≥ kmax
# see Lecture Notes on Elliptic Filter Design, Sophocles J. Orfanidis, 2006
function _K(ks, k′s)
    k = first(ks)
    k′ = first(k′s)
    kₙ = last(ks)
    k′ₙ = last(k′s)

    K = oftype(k, pi)/2
    @inbounds for n in 2:length(ks)
        K *= 1 + ks[n]
    end

    K′ = _large_K_approx(k′ₙ, kₙ)
    denom = oftype(k, 1)
    @inbounds for n in 1:length(k′s)-1
        denom *= 1 + k′s[n]
    end

    return (K, K′/denom)
end

# from Lecture Notes on Elliptic Filter Design, Sophocles J. Orfanidis, 2006
const _kmin = 1e-6
const _kmax = 0.9999999999995
# only applicable if k ≥ _kmax
@inline function _large_K_approx(::T, k′::T) where T
    K = log(T(4)) - log(k′)
    # pretty small numbers, and conv criteria is k′ ≤ √eps(), but ...
    K += (K - 1)/2 * abs2(k′)

    return K
end

################################################################################
# jacobi functions
################################################################################

# TODO
# pq when u ≈ 0
# pq when k² ≥ 0 (or ph(k^2) = π)
# k ≈ 0 or ≈ 1 ? -> Landen{1, Dir}
# check for pole around `n = jK′` for trinity


# approximations
# https://dlmf.nist.gov/22.10

# z⁸ ≈ 0, k⁷ ≤ eps(T)
function _trinity_approx_small_z(z::RealOrComplex, k::RealOrComplex, ::RealOrComplex)
    kk = abs2(k)
    sn_coefs = (0, 1, 0, -(1 + kk)/6, 0, (1 + kk*(14 + kk))/120,
        -(1 + kk*(135 + kk*(135 + kk)))/5040)
    cn_coefs = (1, 0, -0.5, 0, (1+4*kk)/24, 0, -(1 + kk*(44 + 16*kk))/720)
    dn_coefs = (1, 0, -kk/2, 0, kk*(4 + kk)/24, 0, -kk*(16 + kk*(44 + kk))/720)

    return error("do me")
end

# k⁴ ≈ 0, k³ ≤ eps(T)
function _trinity_approx_circular(z::RealOrComplex, k::RealOrComplex, ::RealOrComplex)
    kk = abs2(k)
    sz = sin(z)
    cz = cos(z)
    inner = kk * (z - sz * cz) / 4

    sn = sz - inner * cz
    cn = cz + inner * sz
    dn = 1 - abs2(k* sin(z)) / 2

    return (sn, cn, dn)
end

# k′⁴ ≈ 0, k′³ ≤ eps(T)
function _trinity_approx_hyper(z::RealOrComplex, ::RealOrComplex, k′::RealOrComplex)
    kk = abs2(k′)/4
    sz = sinh(z)
    cz = cos(z)
    scz = sech(z)
    tz = tanh(z)
    sc = sz * cz
    inner = kk * (z - sz * cz) / 4

    sn = tz - kk * (z - sc)*scz^2
    cn = scz + kk * (z - sc) * tz * scz
    dn = scz + kk * (z + sc) * tz * scz

    return (sn, cn, dn)
end

# https://dlmf.nist.gov/22.7

# given a descending sequence of moduli, go back up
function _trinity_descending((sn, cn, dn)::NTuple{3, AbstractFloat}, ks::NTuple)
    N = length(ks)

    for i in N:-1:2
        @inbounds k = ks[i]
        s = (1 + k) * (sn)
        c = cn * dn
        n = 1 + k * sn^2
        d = dn^2  - 1 + k
        nn = 1 + k - dn^2

        sn = s/n
        cn = c/n
        dn = d/nn
    end
end

# given an ascending sequence of moduli, go back down
function _trinity_ascending((sn, cn, dn)::NTuple{3, AbstractFloat}, ks::NTuple)
    N = length(ks)

    for i in N:-1:2
        @inbounds k = ks[i]
        s = (1 + k) * (sn)
        c = cn * dn
        n = 1 + k * sn^2
        d = dn^2  - 1 + k
        nn = 1 + k - dn^2

        sn = s/n
        cn = c/n
        dn = d/nn
    end
end

# Elliptic Functions for Filter Design, Orchard, Willson, IEEE
# Lecture Notes on Elliptic Filter Design, Sophocles J. Orfanidis, 2006

# for both:
# `minus=true` is for `cs` (when using k′)

# start with `k₀ = first(ks) ≈ 0` and `wₙ` → `kₙ = last(ks)` and `w₀`
# ie increasing `kᵢ`, in magnitude, so typically a reversed `ks` of a descending landen seq
function _asc_gauss(w::Number, invw::Number, ks::NTuple, minus::Bool=false)
    # TODO: compile-time versions based on minus. hopefully the compiler can optimize
    N = length(ks)
    coef = ifelse(minus, -1, 1)

    # this does one extra calculation, so if last `k` is the target, skip N
    for n in 1:(N-1)
        k = ks[n]
        num = w + coef * k * invw
        den = 1 + k
        # rather than have to do 1/(1/w) = 1/invw repeatdely
        w = num/den
        invw = den/num
    end

    return w, invw
end

# start with `k₀ = first(ks)` and `w₀` → `kₙ = last(ks) ≈ 0` and `wₙ`
# ie decreasing `kᵢ`, in magnitude, so typically a reversed `ks` of a descending landen seq
function _desc_gauss(w::Number, ::Number, ks::NTuple, minus::Bool=false)
    # TODO: compile-time versions based on minus. hopefully the compiler can optimize
    N = length(ks)

    for n in 2:N
        k = ks[n]
        kp = ks[n-1]

        t = ifelse(minus,
            sqrt(abs2(w) + abs2(kp)),
            sqrt((w-kp) * (w+kp)) )
        w = (1+k)*(w + t)/2
    end

    return w, inv(w)
end




end # module
