module Elliptic
export
    # elliptic integrals of 1st/2nd/3rd kind
    E, F, K, Pi,
    # jacobi elliptic functions
    Jacobi,
    # matlab compatible function calls
    ellipj, ellipke

include("slatec.jl")
include("integrals.jl")
include("jacobi.jl")

end # module
