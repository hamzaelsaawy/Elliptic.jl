#
# jacobi elliptic function approximations
#
# sn, cn, dn:
# https://dlmf.nist.gov/22.10
# https://dlmf.nist.gov/22.13
# https://dlmf.nist.gov/22.5#ii
# am (amplitude):
# http://dlmf.nist.gov/22.16.i

################################################################################
# small z
################################################################################

# |z| ≈ 0, |z|⁸ ≤ eps(T)
# converge when |z| < min(K(k), K′(k)), but π/2 ≤ min(K, K′)
@inline function _sncndn_small_z(z, k)
    m = k^2
    z2 = z^2

    sn = z * @horner(z2, 1, -(1 + m)/6, (1 + m*(14 + m))/120,
        -(1 + m*(135 + m*(135 + m)))/5040)
    cn = @horner(z2, 1, -0.5, (1+4*m)/24, -(1 + m*(44 + 16*m))/720)
    dn = @horner(z2, 1, -m/2, m*(4 + m)/24, -m*(16 + m*(44 + m))/720)

    return (sn, cn, dn)
end

# z ≈ 0, O(z⁷), |k|⁶ ≤ eps
@inline function _am_small_z(z, k)
    m = k^2
    z2 = z^2
    return z * @horner(z, 1, -m/6, m*(4+k)/120)
end

################################################################################
# small k or k′
# most equations below are unstable when 2z ≈ sin(2*z) (or sinh), ie z ≈ 0
################################################################################

# k ≈ 0, O(k⁴), |k|³ ≤ eps; sn (& cn) can be ≤ k
# use when |k|² ≤ eps
@inline function _sncndn_circular(z, k)
    s, c = sincos(z)
    inner = k^2 * (z - s * c) / 4
    sqrt2 = √2

    sn = s - inner * c
    cn = c + inner * s
    dn =  (sqrt2 - k * s) * (sqrt2 + k * s) / 2

    return (sn, cn, dn)
end

# k′ ≈ 0, O(k′⁴), so |k′|³ ≤ eps?
# use when |k′|² ≤ eps
@inline function _sncndn_hyper(z, k′)
    m′ = k′^2 / 4
    c = cosh(z)
    t = tanh(z)
    s2 = sinh(2*z)

    sn = t - m′ * (2*z - s2) / (8 * c^2)
    cn = (8 + m′ * (2*z - s2) * t) / (8 * c)
    dn = (8 + m′ * (2*z + s2) * t) / (8 * c)

    return (sn, cn, dn)
end

# k ≈ 0, O(k⁴), |k|³ ≤ eps
_am_circular(z, k) = z - k^2 * (2*z - sin(2*z)) / 8

# k′ ≈ 0, O(k′⁴), |k′|³ ≤ eps
_am_hyper(z, k′) = gd(z) - k′^2 * (2*x - sinh(2*z)) / (8 * cosh(z))
