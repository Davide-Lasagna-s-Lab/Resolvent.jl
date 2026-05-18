# Utility method to generate modes for a set of parameters

export generate_modes!, generate_modes

function generate_modes!(Ψ::NTuple{3, Array{Complex{T}, 5}},
                        Re::T,
                        Ro::T,
                        Dy::AbstractMatrix{T},
                       Dy2::AbstractMatrix{T},
                        ws::Vector{T},
                         α::T,
                         β::T,
                         ω::T;
                      base::Vector{T}=ones(T, size(Ψ[1], 1)),
                   verbose::Bool=true) where {T}
    # unpack useful variables from inputs
    Ny, M, Nx, Nz, Nt = size(Ψ[1])

    # initialise resolvent operator
    H = Resolvent(Ny, Dy, Dy2)

    # loop over frequencies computing response modes
    verbose && println("0:$(Nx-1), $(-(Nz >> 1)):$(Nz >> 1), $(-(Nt >> 1)):$(Nt >> 1)")
    for nt in -(Nt >> 1):(Nt >> 1), nz in -(Nz >> 1):(Nz >> 1), nx in 0:Nx-1
        verbose && print("$nx, $nz, $nt                  \r"); flush(stdout)
        _nx = nx + 1
        _nz = nz >= 0 ? nz + 1 : Nz + nz + 1
        _nt = nt >= 0 ? nt + 1 : Nt + nt + 1
        U = svd(H(nx*α, nz*β, nt*ω, base, Re, Ro), ws, M).U
        Ψ[1][:, :, _nx, _nz, _nt] .= U[1:Ny, :]
        Ψ[2][:, :, _nx, _nz, _nt] .= U[(Ny + 1):(2*Ny), :]
        Ψ[3][:, :, _nx, _nz, _nt] .= U[(2*Ny + 1):(3*Ny), :]
    end

    return Ψ
end

generate_modes(S, M, Re, Ro, Dy, Dy2, ws, α, β, ω, ::Type{T}=Float64; base=ones(T, S[1]), verbose=true) where {T} =
    generate_modes!(ntuple(_ -> zeros(Complex{T}, S[1], M, (S[2] >> 1) + 1, S[3], S[4]), 3),
                    T(Re), T(Ro), T.(Dy), T.(Dy2), T.(ws), T(α), T(β), T(ω), base=T.(base), verbose=verbose)
