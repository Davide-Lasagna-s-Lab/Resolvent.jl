using Test
using Random
using LinearAlgebra

using ResolventAnalysis
using ChebUtils

# # set up the PyCall interface to be able to use Sean's Matlab code directly
using PyCall
o = pyimport("oct2py").octave
o.addpath("./primitive_resolvent/")

@testset "Truncation performed correctly    " begin
    # contruct SVD of random matrix and perform truncation
    trunc_length = rand(1:4)
    A = rand(rand(5:10), rand(5:10))
    U, S, V = svd(A)
    svd_trunc = ResolventAnalysis.truncateSVD(svd(A), trunc_length)

    # check SVDs are the same as their individually truncated versions
    @test svd_trunc.U == U[:, 1:trunc_length]
    @test svd_trunc.S == S[1:trunc_length]
    @test svd_trunc.Vt == V[:, 1:trunc_length]'
end

@testset "Cholesky decomposition of weights " begin
    # generate weights
    N = rand(3:100)
    ws = chebws(N)

    # compute cholesky decomposition matrices
    L, L_inv = ResolventAnalysis.cholesky(ws)

    Z = zeros(N, N)
    @test L*L' ≈   [Diagonal(ws) Z            Z           ;
                    Z            Diagonal(ws) Z           ;
                    Z            Z            Diagonal(ws);]
    
    @test L_inv*L ≈    [I Z Z Z;
                        Z I Z Z;
                        Z Z I Z;]
end

@testset "Resolvent at specific mode number " begin
    # generate random flow parameters
    Re = rand()*200
    kz = rand()*100
    kt = rand()*100
    N = rand(4:200)

    # compute Resolvents
    H0 = o.quick_example(Re, kz, kt, N, nout=1)
    H_me = Resolvent(N, chebdiff(N), chebddiff(N))

    @test H_me(kz, kt, ones(N), Re, 0.0) ≈ H0[1:4*N, 1:3*N]
end

@testset "Matrix SVD at specific mode number" begin
    # generate random flow parameters
    Re = rand()*200
    kz = rand()*100
    kt = rand()*100
    N = rand(4:200)
    ws = chebws(N)

    # compute Resolvents
    _, U_sean, S_sean, V_sean = o.quick_example(Re, kz, kt, N, nout=4); S_sean = reshape(S_sean, 4*N)
    H_me = Resolvent(N, chebdiff(N), chebddiff(N))
    SVD_me = svd(H_me(kz, kt, ones(N), Re, 0.0), ws, 5)

    @test abs.(SVD_me.U) ≈ abs.(U_sean[1:3*N, 1:5])
    @test SVD_me.S ≈ S_sean[1:5]
    @test abs.(SVD_me.V) ≈ abs.(V_sean[1:3*N, 1:5])
end

@testset "generate_modes components checked against reference" begin
    N   = rand(4:20)
    M   = 3
    ws  = chebws(N)
    Dy  = chebdiff(N)
    Dy2 = chebddiff(N)
    Re  = rand()*200
    Ro  = 0.0
    β   = rand(); ω = rand()
    S   = (N, 1, 4, 4)  # Nx=1 so kx=0 throughout, matching quick_example

    Ψ = generate_modes(S, M, Re, Ro, Dy, Dy2, ws, 0.0, β, ω; verbose=false)

    # structure: tuple of three (Ny × M × Nx' × Nz × Nt) arrays
    @test Ψ isa NTuple{3, Array{ComplexF64, 5}}
    @test all(size(Ψ[n]) == (N, M, 1, S[3], S[4]) for n in 1:3)

    # compare u/v/w components at nz=1, nt=1 against Sean's reference SVD
    nz = 1; nt = 1
    _, U_sean, _, _ = o.quick_example(Re, nz*β, nt*ω, N, nout=4)
    @test abs.(Ψ[1][:, :, 1, nz + 1, nt + 1]) ≈ abs.(U_sean[1:N,          1:M])
    @test abs.(Ψ[2][:, :, 1, nz + 1, nt + 1]) ≈ abs.(U_sean[(N + 1):2N,   1:M])
    @test abs.(Ψ[3][:, :, 1, nz + 1, nt + 1]) ≈ abs.(U_sean[(2N + 1):3N,  1:M])
end
