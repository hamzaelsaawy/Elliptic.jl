module Elliptic


# elliptic integrals of 1st/2nd/3rd kind
export E, F, K, Pi
# jacobi elliptic functions
export Jacobi
# matlab compatible
export ellipj, ellipke

include("slatec.jl")
include("integrals.jl")

include("landen.jl")
include("jacobi.jl")

end # module
