using SparseArrays
using LinearAlgebra
using ITensors
using HDF5

include("src/Lattice.jl")
include("src/Parameter.jl")


para = GetParameter("input.inp")
Nsite, mesh, nn, _, _ = Honeycomb(para)
Nstates = para.Nstates  #number of EXCITED states to probe

# @time Ham = Hamiltonian(param)
# @time eigvals, eigvecs = eigs(Ham, nev = Nstates, which=:SR, tol = 0.001) # SR = Smallest Real 
# println("\nEigen values:"); show(stdout, "text/plain", real.(eigvals)); println()

let
    sites = siteinds("S=1/2", Nsite)
    ampo = OpSum()

    # coupling
    for j=1:Nsite
        nx::Int64 = nn[j, 1]
        ny::Int64 = nn[j, 2]
        nz::Int64 = nn[j, 3]
        
        if j < nn[j, 1]
            ampo += para.Kxx, "Sx", j, "Sx", nx
        end
        if j < nn[j, 2]
            ampo += para.Kyy, "Sy", j, "Sy", ny
        end
        if j < nn[j, 3]
            ampo += para.Kzz, "Sz", j, "Sz", nz
        end
    end

    # magnetic field
    for j=1:Nsite
        ampo += para.Hx, "Sx", j
        ampo += para.Hy, "Sy", j
        ampo += para.Hz, "Sz", j
    end
    
    H = MPO(ampo, sites)

    sweeps = Sweeps(12)
    setmaxdim!(sweeps, 100, 200, 200, 300, 400, 500, 500, 500, 500, 500, 600, 600)
    setcutoff!(sweeps, 1E-8)

    f = h5open("data.h5","r")
    @time ψ = read(f,"psi",MPS)
    close(f)
    ψ = replace_siteinds(ψ,sites)

    @show energy0 = inner(ψ', H, ψ)

    # excited states
    weight::Float64 = 1.0
    es = [real(energy0)]
    ψs = [ψ]
    for i=1:Nstates
        ψ_init = randomMPS(sites, 100)
        ψ_lower = [ψs[j] for j=1:i]
        energy, ψ = dmrg(H, ψ_lower, ψ_init, sweeps; weight)
        push!(es, real(energy))
        push!(ψs, ψ)
    
        f = h5open("Excited/data"*string(i)*".h5","w")
        write(f, "psi", ψ)
        close(f)
    end

    @show es


end





