# Utility method to generate modes for a set of parameters

export generate_modes!, generate_modes

function generate_modes!(Ψ::NTuple{3, Array{ComplexF64, 5}},
                        Re::Real,
                        Ro::Real,
                        Dy::AbstractMatrix{<:Real},
                       Dy2::AbstractMatrix{<:Real},
                        ws::Vector{<:Real},
                         α::Real,
                         β::Real,
                         ω::Real;
                      base::Vector{<:Real}=ones(size(Ψ[1], 1)),
                   verbose::Bool=true)
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

"""
    generate_modes(Nx, Nz, Nt, M, Re, Ro, Dy, Dy2, ws, α, β, ω; base, verbose) -> NTuple{3, Array{ComplexF64, 5}}

Compute the `M` leading resolvent response modes at every wavenumber on the
grid `(kx, kz, kt) = (0:Nx-1, -(Nz÷2):Nz÷2, -(Nt÷2):Nt÷2)` scaled by `(α, β, ω)`.

Returns a tuple `(Ψu, Ψv, Ψw)` where each array has shape `(Ny, M, Nx÷2+1, Nz, Nt)`:
- axis 1: wall-normal points
- axis 2: mode index (ordered by singular value, descending)
- axes 3–5: rfft storage for kx, and full FFT storage for kz and kt

# Arguments
- `Nx, Nz, Nt`: number of streamwise, spanwise, and temporal wavenumbers
- `M`: number of modes to retain
- `Re, Ro`: Reynolds and Rotation numbers
- `Dy, Dy2`: first and second wall-normal derivative matrices (size `Ny × Ny`)
- `ws`: wall-normal quadrature weights (length `Ny`)
- `α, β, ω`: wavenumber scales for x, z, t directions
- `base`: mean streamwise velocity gradient profile (default: uniform flow `ones(Ny)`)
- `verbose`: print progress to stdout (default: `true`)
"""
function generate_modes(Nx, Nz, Nt, M, Re, Ro, Dy, Dy2, ws, α, β, ω;
                        base=ones(size(Dy, 1)), verbose=true)
    Ny = size(Dy, 1)
    return generate_modes!(ntuple(_ -> zeros(ComplexF64, Ny, M, (Nx >> 1) + 1, Nz, Nt), 3),
                           Re, Ro, Dy, Dy2, ws, α, β, ω, base=base, verbose=verbose)
end
