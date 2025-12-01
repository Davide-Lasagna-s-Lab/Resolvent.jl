# Definition of the function to compute the resolvent matrix at a specific
# frequency.

struct Resolvent{Ny}
    H::Matrix{ComplexF64}
    H_inv::Matrix{ComplexF64}
    M::Matrix{Float64}
    Δ::Matrix{Float64}
    Dy::Matrix{Float64}
    Dy2::Matrix{Float64}
    id::Diagonal{Bool, Vector{Bool}}

    function Resolvent(Ny::Int, Dy, Dy2)
        Z = zeros(Ny, Ny)
        M = [I Z Z;
             Z I Z;
             Z Z I;
             Z Z Z]
        M[1, :] .= 0.0
        M[Ny:(Ny + 1), :] .= 0.0
        M[(2*Ny):(2*Ny + 1), :] .= 0.0
        M[3*Ny, :] .= 0.0
        new{Ny}(zeros(ComplexF64, 4*Ny, 3*Ny), zeros(ComplexF64, 4*Ny, 4*Ny), M, zeros(ComplexF64, Ny, Ny), Dy, Dy2, I(Ny))
    end
end
Base.size(::Resolvent{Ny}) where {Ny} = Ny

(f::Resolvent)(kz, kt, dūdy, Re, Ro) = f(0, kz, kt, dūdy, Re, Ro)

function (f::Resolvent{Ny})(kx, kz, kt, dūdy, Re, Ro) where {Ny}
    # define laplacian operator
    f.Δ .= .-f.id.*(kx^2) .+ f.Dy2 .- f.id.*(kz^2)

    # fill resolvent matrix
    f.H_inv[1:Ny, 1:Ny]                           .= 1im.*kt.*f.id   .- f.Δ./Re
    f.H_inv[1:Ny, (Ny + 1):(2*Ny)]                .= Diagonal(dūdy) .- f.id.*Ro
    f.H_inv[1:Ny, (3*Ny + 1):end]                 .= 1im.*kx.*f.id
    f.H_inv[(Ny + 1):(2*Ny), 1:Ny]                .= f.id.*Ro
    f.H_inv[(Ny + 1):(2*Ny), (Ny + 1):(2*Ny)]     .= 1im.*kt.*f.id   .- f.Δ./Re
    f.H_inv[(Ny + 1):(2*Ny), (3*Ny + 1):end]      .= f.Dy
    f.H_inv[(2*Ny + 1):(3*Ny), (2*Ny + 1):(3*Ny)] .= 1im.*kt.*f.id  .- f.Δ./Re
    f.H_inv[(2*Ny + 1):(3*Ny), (3*Ny + 1):end]    .= 1im.*kz.*f.id
    f.H_inv[(3*Ny + 1):end, 1:Ny]                 .= -1im.*kx.*f.id
    f.H_inv[(3*Ny + 1):end, (Ny + 1):(2*Ny)]      .= .-f.Dy
    f.H_inv[(3*Ny + 1):end, (2*Ny + 1):(3*Ny)]    .= -1im.*kz.*f.id

    # apply boundary conditions
    f.H_inv[1, :]                 .= 0.0
    f.H_inv[1, 1]                  = 1.0
    f.H_inv[Ny:(Ny + 1), :]       .= 0.0
    f.H_inv[Ny, Ny]                = 1.0
    f.H_inv[Ny + 1, Ny + 1]        = 1.0
    f.H_inv[(2*Ny):(2*Ny + 1), :] .= 0.0
    f.H_inv[2*Ny, 2*Ny]            = 1.0
    f.H_inv[2*Ny + 1, 2*Ny + 1]    = 1.0
    f.H_inv[3*Ny, :]              .= 0.0
    f.H_inv[3*Ny, 3*Ny]            = 1.0

    # invert resolvent and multiply by influence matrix
    mul!(f.H, inv(f.H_inv), f.M)

    return f.H
end
