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
    for kt in -(Nt >> 1):(Nt >> 1), kz in -(Nz >> 1):(Nz >> 1), kx in 0:Nx-1
        verbose && print("$kx, $kz, $kt                  \r"); flush(stdout)
        ind_kx = kx + 1
        ind_kz = kz >= 0 ? kz + 1 : Nz + kz + 1
        ind_kt = kt >= 0 ? kt + 1 : Nt + kt + 1
        U = svd(H(kx*α, kz*β, kt*ω, base, Re, Ro), ws, M).U
        Ψ[1][:, :, ind_kx, ind_kz, ind_kt] .= U[1:Ny, :]
        Ψ[2][:, :, ind_kx, ind_kz, ind_kt] .= U[(Ny + 1):(2*Ny), :]
        Ψ[3][:, :, ind_kx, ind_kz, ind_kt] .= U[(2*Ny + 1):(3*Ny), :]
    end

    return Ψ
end

generate_modes(S, M, Re, Ro, Dy, Dy2, ws, α, β, ω, ::Type{T}=Float64; base=ones(T, S[1]), verbose=true) where {T} =
    generate_modes!(ntuple(_ -> zeros(Complex{T}, S[1], M, (S[2] >> 1) + 1, S[3], S[4]), 3),
                    T(Re), T(Ro), T.(Dy), T.(Dy2), T.(ws), T(α), T(β), T(ω), base=T.(base), verbose=verbose)
