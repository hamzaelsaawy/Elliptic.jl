#
# miscellaneous functions
#

import Base.Math.@horner

const Maybe{T} = Union{Nothing, T}
const RealOrComplexFloat{T<:AbstractFloat} = Union{T, Complex{T}}

Base.ldexp(z::Complex, n::Integer) =
        Complex(ldexp(real(z), n), ldexp(imag(z), n))

# do everything as a float
_one_eps(::Type{T}=Float64) where T<:Number = eps(float(real(T)))
_one_eps(x::Number) = _one_eps(typeof(x))

_float_typeof(x::Number) = float(typeof(x))
_float_typeof(::Nothing) = Nothing
_promote_float_typeof(x) = _float_typeof(x)
@inline _promote_float_typeof(x, xs...) =
        promote_type(_promote_float_typeof(x), _promote_float_typeof(xs...))
_promote_float(ts...) = convert.(_promote_float_typeof(ts...), ts)

# rotate numbers by pi to have phase in (-pi/2, pi/2)
_abs_real(x::Real) = abs(x)
# if real(x) == 0, can either be ±pi/2 ... wont check for that
_abs_real(x::Complex) = ifelse(real(x) < 0, -x, x)

 _k′(k) = sqrt((1-k)*(1+k))

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
