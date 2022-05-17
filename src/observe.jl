using LinearAlgebra
using ITensors
using HDF5

let

    f = h5open("../data/data56.h5","r")
    @time ψ = read(f,"psi",MPS)
    close(f)
    sites = siteinds(ψ) 
    Nsite = 56    
#     # Magnetization
#     @time magz = expect(ψ,"Sz")
#     for (j,mz) in enumerate(magz)
#         println("$j $mz")
#     end
#     # Correlation x-x
#     @time Cx = correlation_matrix(ψ,"Sx","Sx") 
#     show(stdout, "text/plain", real(Cx)); println()






  #https://itensor.github.io/ITensors.jl/stable/examples/MPSandMPO.html
    Mag = Matrix{Float64}(undef, Nsite, 3)
    for i=1:Nsite
        orthogonalize!(ψ, i)
        ket::ITensor = ψ[i]
        bra::ITensor = dag(prime(ket, "Site"))
        opSx::ITensor = op(sites, "Sx", i)
        opSy::ITensor = op(sites, "Sy", i)
        opSz::ITensor = op(sites, "Sz", i)

        sx = scalar(bra * opSx * ket)
        sy = scalar(bra * opSy * ket)
        sz = scalar(bra * opSz * ket)

        @show Mag[i, 1] = real(sx)
        @show Mag[i, 2] = real(sy)
        @show Mag[i, 3] = real(sz)
    end
	show(stdout, "text/plain", Mag); println()
# 
#     orthogonalize!(ψ, 1)
#     op1 = op(sites, "Sx", 1)
#     op2 = op(sites, "Sx", 2)
#     C::ITensor = ψ[1]
#     C *= op1
#     ir = commonind(ψ[1],ψ[1+1],"Link")
#     C *= dag(prime(prime(ψ[1],"Site"),ir))
# 
#     C *= ψ[2]
#     C *= op2
#     jl = commonind(ψ[2], ψ[2-1], "Link")
#     C *= dag(prime(prime(ψ[2], "Site"), jl))
#     res = scalar(C)
# 
#     # correlation x-x
#     Cx = Matrix{Float64}(undef, Nsite, Nsite)
#     for i=1:Nsite
#         orthogonalize!(ψ, i)
#         op_i::ITensor = op(sites, "Sx", i)
#         for j=i+1:Nsite
#             println(i, " ", j)
#             op_j::ITensor = op(sites, "Sx", j)
#             C::ITensor = ψ[i]
#             C *= op_i
#             ir = commonind(ψ[i],ψ[i+1],"Link")
#             C *= dag(prime(prime(ψ[i],"Site"),ir))
# 
#             for k=i+1:j-1
#                 println("k=", k)
#                 C *= ψ[k]
#                 C *= dag(prime(ψ[k], "Link"))
#             end
# 
#             C *= ψ[j]
#             C *= op_j
#             jl = commonind(ψ[j], ψ[j-1], "Link")
#             C *= dag(prime(prime(ψ[j], "Site"), jl))
# 
#             res = scalar(C)
# #                 Cx[i, j] = res
# #                 println(res)
#         end
#     end
end
