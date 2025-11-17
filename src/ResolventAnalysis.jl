module ResolventAnalysis

using LinearAlgebra, TSVD

export Resolvent, svd

export DivideAndConquer, QRIteration, Lanczos, Adaptive

include("resolvent.jl")
include("cholesky.jl")
include("svd.jl")
include("generate_modes.jl")

end
