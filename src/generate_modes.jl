# Utility method to generate modes for a set of parameters

export generate_modes!, generate_modes

function generate_modes!(Ψ::Array{Complex{T}},
                        Re::T,
                        Ro::T,
                        Dy::AbstractMatrix{T},
                       Dy2::AbstractMatrix{T},
                        ws::Vector{T},
                         β::T,
                         ω::T;
                      base::Vector{T}=ones(T, size(Ψ, 1)÷3),
                   verbose::Bool=true) where {T}
    # unpack useful variables from inputs
    N, M, Nz, Nt = size(Ψ)

    # initialise resolvent operator
    H = Resolvent(N÷3, Dy, Dy2)

    # loop over frequencies computing response modes
    for nt in 1:(Nt >> 1) + 1, nz in 1:Nz
        verbose && print("$nz/$Nz, $nt/$((Nt >> 1) + 1)                    \r")
        if nt == 1
            Ψ[:, :, nz, 1]        .= svd(H((nz - 1)*β, 0,          base, Re, Ro), ws, M).U
        elseif nz == 1
            Ψ[:, :, 1, nt]        .= svd(H(0,          (nt - 1)*ω, base, Re, Ro), ws, M).U
            Ψ[:, :, 1, end-nt+2]  .= conj.(Ψ[:, :, 1, nt])
        else
            Ψ[:, :, nz, nt]       .= svd(H((nz - 1)*β, (nt - 1)*ω, base, Re, Ro), ws, M).U
            Ψ[:, :, nz, end-nt+2] .= svd(H((nz - 1)*β, (1 - nt)*ω, base, Re, Ro), ws, M).U
        end
        flush(stdout)
    end

    return Ψ
end

generate_modes(S, M, Re, Ro, Dy, Dy2, ws, β, ω, ::Type{T}=Float64; base=ones(T, S[1]), verbose=true) where {T} = 
    generate_modes!(zeros(Complex{T}, 3*S[1], M, (S[2] >> 1) + 1, S[3]), T(Re), T(Ro), T.(Dy), T.(Dy2), T.(ws), T(β), T(ω), base=T.(base), verbose=verbose)

function generate_modes(path)
    throw(error("Not implemented"))
end

function load_mdoes!(Ψ, path) end # read file and assign to array
load_modes(path, S, M, ::Type{T}=Float64) where {T} = load_modes!(zeros(Complex{T}, S[1], M, (S[2] >> 1) + 1, S[3]), path)
function load_modes(path) end # create a mmap of the array

function write_modes(path, Ψ) end
