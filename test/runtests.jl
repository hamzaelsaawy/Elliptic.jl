using Test
using SpecialFunctions: gamma
using DelimitedFiles: readdlm

using Elliptic

include("integrals_tests.jl")
include("jacobi_tests.jl")
include("landen_tests.jl")
