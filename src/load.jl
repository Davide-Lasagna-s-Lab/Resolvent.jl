# Functionality to read and write modes from disk.

function slab_modes(Ψ, np, rank)
    # initialise slab of modes
    Ny    = div(size(Ψ, 1), 3)
    Ny_sb = div(Ny, np)
    Ψ_sb  = zeros(eltype(Ψ), Ny_sb*3, size(Ψ)[2:end]...)

    # copy data from complete modes
    @views Ψ_sb[(1        ):1*Ny_sb, :, :, :, :] .= Ψ[(1+     rank*Ny_sb):(     (rank+1)*Ny_sb), :, :, :, :]
    @views Ψ_sb[(1+  Ny_sb):2*Ny_sb, :, :, :, :] .= Ψ[(1+  Ny+rank*Ny_sb):(  Ny+(rank+1)*Ny_sb), :, :, :, :]
    @views Ψ_sb[(1+2*Ny_sb):3*Ny_sb, :, :, :, :] .= Ψ[(1+2*Ny+rank*Ny_sb):(2*Ny+(rank+1)*Ny_sb), :, :, :, :]

    return Ψ_sb
end

function load_modes(np, rank, path)
    open(path, "r") do f
        # read type and size data
        T     = eval(Symbol(readuntil(f, "\n")))
        shape = ntuple(_->read(f, Int), 5)

        # initialise slab of modes
        Ny    = div(shape[1], 3)
        Ny_sb = div(Ny, np)
        Ψ_sb  = zeros(T, Ny_sb*3, shape[2:end]...)

        # copy data from file
        skip(f, rank*Ny_sb*sizeof(T)) # initial offset according the MPI rank
        for I in CartesianIndices(shape[2:end])
            for n in 0:2
                for ny in 1:Ny_sb
                    Ψ_sb[ny + n*Ny_sb, I] = read(f, T)
                end
                skip(f, (Ny - Ny_sb)*sizeof(T))
            end
        end

        return Ψ_sb
    end
end

function load_modes!(Ψ, path)
    open(path, "r") do f
        # go past size and type data
        readuntil(f, "\n")
        [read(f, Int) for _ in 1:5]

        # read data
        read!(f, Ψ)
    end
    return Ψ
end
function load_modes(path)
    open(path, "r") do f
        # read type and size data
        T     = eval(Symbol(readuntil(f, "\n")))
        shape = ntuple(_->read(f, Int), 5)

        # read data
        Ψ = zeros(T, shape...)
        read!(f, Ψ)

        return Ψ
    end
end

function save_modes(Ψ, path)
    open(path, "w") do f
        # write data type
        write(f, string(eltype(Ψ))*"\n")

        # write size of data
        [write(f, size(Ψ, n)) for n in 1:5]

        # write data
        write(f, Ψ)
    end
    return nothing
end
