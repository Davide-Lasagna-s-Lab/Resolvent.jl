# Definitions for the custom SVD of the resolvent matrix

const DivideAndConquer = LinearAlgebra.DivideAndConquer
const QRIteration = LinearAlgebra.QRIteration
struct Lanczos <: LinearAlgebra.Algorithm end
struct Arnoldi <: LinearAlgebra.Algorithm end
struct Adaptive <: LinearAlgebra.Algorithm end


LinearAlgebra.svd(H::Matrix{ComplexF64}, ws::Vector{Float64}, nvals::Int=size(H, 2); alg::LinearAlgebra.Algorithm=DivideAndConquer(), debug::Bool=false) =
    _svd(H, cholesky(ws)..., nvals, alg, debug)

function _svd(H, L, L_inv, nvals, ::Adaptive, debug)
    try
        return _mysvd(H, L, L_inv, nvals, DivideAndConquer(), debug)
    catch end

    try
        return _mysvd(H, L, L_inv, nvals, QRIteration(), debug)
    catch end

    try
        return _mysvd(H, L, L_inv, nvals, Lanczos(), debug)
    catch end

    try
        return _mysvd(H, L, L_inv, nvals, Arnoldi(), debug)
    catch end
    
    throw(LinearAlgebra.LAPACKException)
end

function _svd(H, L, L_inv, nvals, alg::LinearAlgebra.Algorithm, debug)
    # perform SVD
    decomp = _svd(L*H*L_inv, nvals, alg, debug)

    # scale outputs
    mul!(decomp.U, L_inv, copy(decomp.U))
    mul!(decomp.Vt, copy(decomp.Vt), L_inv)

    return decomp
end

_svd(H::Matrix{ComplexF64}, nvals, ::Lanczos, debug) = ((U, S, V) = tsvd(H, nvals, debug=debug); return SVD(U, S, V'))
_svd(H::Matrix{ComplexF64}, nvals, ::Arnoldi, debug) = ((U, S, V) = svds(H, nsv=nvals, tol=1e-12, maxiter=1000); return SVD(U, S, V'))
_svd(H::Matrix{ComplexF64}, nvals, alg::LinearAlgebra.Algorithm, ::Bool) = truncateSVD(svd(H, alg=alg), nvals)

truncateSVD(svd::SVD, trunc::Int) = SVD(svd.U[:, 1:trunc], svd.S[1:trunc], svd.V[:, 1:trunc]')
