module Elliptic


# elliptic integrals of 1st/2nd/3rd kind
export E, F, K, Pi
# jacobi elliptic functions
export Jacobi, Landen
# matlab compatible
export ellipj, ellipke

include("slatec.jl")
include("integrals.jl")

include("landen.jl")
include("jacobi.jl")
include("jacobi.old.jl")

end # module
