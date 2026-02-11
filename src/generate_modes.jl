# Utility method to generate modes for a set of parameters

export generate_modes!, generate_modes

function generate_modes!(Ψ::Array{Complex{T}, 5},
                        Re::T,
                        Ro::T,
                        Dy::AbstractMatrix{T},
                       Dy2::AbstractMatrix{T},
                        ws::Vector{T},
                         α::T,
                         β::T,
                         ω::T;
                      base::Vector{T}=ones(T, size(Ψ, 1)÷3),
                   verbose::Bool=true) where {T}
    # unpack useful variables from inputs
    N, M, Nx, Nz, Nt = size(Ψ)

    # initialise resolvent operator
    H = Resolvent(N÷3, Dy, Dy2)

    # loop over frequencies computing response modes
    verbose && println("0:$(Nx-1), $(-(Nz >> 1)):$(Nz >> 1), $(-(Nt >> 1)):$(Nt >> 1)")
    for nt in -(Nt >> 1):(Nt >> 1), nz in -(Nz >> 1):(Nz >> 1), nx in 0:Nx-1
        verbose && print("$nx, $nz, $nt                  \r"); flush(stdout)
        _nx = nx + 1
        _nz = nz >= 0 ? nz + 1 : Nz + nz + 1
        _nt = nt >= 0 ? nt + 1 : Nt + nt + 1
        Ψ[:, :, _nx, _nz, _nt] .= svd(H(nx*α, nz*β, nt*ω, base, Re, Ro), ws, M).U
    end

    return Ψ
end

generate_modes(S, M, Re, Ro, Dy, Dy2, ws, α, β, ω, ::Type{T}=Float64; base=ones(T, S[1]), verbose=true) where {T} = 
    generate_modes!(zeros(Complex{T}, 3*S[1], M, (S[2] >> 1) + 1, S[3], S[4]), T(Re), T(Ro), T.(Dy), T.(Dy2), T.(ws), T(α), T(β), T(ω), base=T.(base), verbose=verbose)
