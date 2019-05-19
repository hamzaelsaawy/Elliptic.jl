#
# jacobi elliptic functions
#
# sn, cn, dn:
# https://dlmf.nist.gov/22.2
# am (amplitude):
# http://dlmf.nist.gov/22.16.i

module Jacobi

import ..Landen
using ..Landen: LandenSeq, NonConvergedLandenSeq,
        AscendingLanden, DescendingLanden, direction, ktype, K

export am, sncndn,
    sn, cn, dn

include("misc.jl")
include("jacobi_approx.jl")

################################################################################
# amplitude
################################################################################

@inline function am(z, k, k′=_k′(k))
    z, k = _promote_float(z, k)
    tol = _one_eps(z)

    return if abs2(z)^3 ≤ tol
        _am_small_z(z, k)
    elseif abs(k)^3 ≤ tol
        _am_circular(z, k)
    elseif abs(k′)^3 ≤ tol
        _am_hyper(z, k′)
    else
        asin(sn(z, k))
    end
end

################################################################################
# jacobi landen/gauss sequence
################################################################################

# TODO
# pq when k² ≥ 0 (or ph(k^2) = π)
# check for pole around `n = jK′` for sncndn
# special versions of sn,dn,cn iterations?

sn(z, k) = sncndn(z, k)[1]
cn(z, k) = sncndn(z, k)[2]
dn(z, k) = sncndn(z, k)[3]

function sncndn(z, k::Number, k′::Maybe{Number}=missing)
    ma = abs2(z)
    swapped = false

    if ma ≥ 1
        z /= k
        k, k′ = 1/k, im*k′/k
        swapped = true
    end

    s, c, d = sncndn(z, LandenSeq(k, k′))

    return (swapped ? (k*s, d, c) : (s, c, d))
end

sncndn(_, ::NonConvergedLandenSeq) = NaN, NaN, NaN
function sncndn(z, landen::LandenSeq)
    T = promote_type(_float_typeof(z), ktype(landen))
    z = T(z)

    (abs2(z)^4 ≤ eps(T)) && return _sncndn_small_z(z, last(landen.ks))

    return _sncndn_seq(z, landen)
end

_sncndn_seq(z, landen::LandenSeq{1}) = _sncndn_seq_end(z, landen)

@inline function _sncndn_seq(z, landen::LandenSeq)
    N = length(landen)
    D = direction(landen)

    ks = landen.ks
    k′s = landen.k′s

    KK, _ = K(landen)
    sn, cn, dn = _sncndn_seq_end(z * π / (2 * KK), landen)

    for i in N:-1:2
        @inbounds k = ks[i]
        @inbounds k′ = k′s[i]

        sn, cn, dn = _sncndn_seq_iter(sn, cn, dn, k, k′, D)
    end

    return (sn, cn, dn)
end

_sncndn_seq_end(z, landen::LandenSeq{N,T,DescendingLanden}) where {N,T} =
        _sncndn_circular(z, last(landen.ks))

_sncndn_seq_end(z, landen::LandenSeq{N,T,AscendingLanden}) where {N,T} =
        _sncndn_hyper(z, last(landen.k′s))

#
# landen/gauss transformations
#
# https://dlmf.nist.gov/22.7

# given a descending sequence of moduli, go back up
@inline function _sncndn_seq_iter(sn, cn, dn, k, k′, ::Type{DescendingLanden})
    n = 1 + k * sn^2
    s = (1 + k) * sn / n
    c = cn * dn / n
    # avoid cancelation
    sn2 = sn^2
    d = (1 - k * sn2) / (1 + k * sn2)
    # d = hypot(k′ * s, c)
    # d = (dn^2  - 1 + k) / (1 - dn^2 + k)

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

end # module
