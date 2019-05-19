module Elliptic
# elliptic integrals of 1st/2nd/3rd kind
export E, F, K, Pi

# jacobi elliptic functions
export Jacobi

# matlab compatible
export ellipj, ellipke

include("slatec.jl")
include("integrals.jl")
include("jacobi.jl")
include("landen.jl")

end # module
