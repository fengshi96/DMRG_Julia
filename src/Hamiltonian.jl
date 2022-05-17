using LinearAlgebra
using SparseArrays
using Arpack
include("Lattice.jl")
include("Dofs.jl")
include("Parameter.jl")

function Hamiltonian(param::Parameter)  #LLX, LLY, BC::String = "OBC"

	Kxx::Float64 = param.Kxx 
	Kyy::Float64 = param.Kyy 
	Kzz::Float64 = param.Kzz
	Hx::Float64 = param.Hx 
	Hy::Float64 = param.Hy
	Hz::Float64 = param.Hz
	
	nsite, mesh_, nn_, indx_, indy_= Honeycomb(param)
	KxxGraph_ = fill!(Matrix{Float64}(undef, nsite, nsite),0)
	KyyGraph_ = fill!(Matrix{Float64}(undef, nsite, nsite),0)
	KzzGraph_ = fill!(Matrix{Float64}(undef, nsite, nsite),0)
	
	for i::Int8 in 1:nsite
		# Kxx_Conn
		j::Int8 = nn_[i,1]
		if i < j && j>=1
			KxxGraph_[i,j] = Kxx
			KxxGraph_[j,i] = Kxx
		end
	
		# Kyy_Conn
		j = nn_[i,2]
		if i < j && j>=1
			KyyGraph_[i,j] = Kyy
			KyyGraph_[j,i] = Kyy
		end	

		# Kzz_Conn
		j = nn_[i,3]
		if i < j && j>=1
			KzzGraph_[i,j] = Kzz
			KzzGraph_[j,i] = Kzz
		end
	
	end
	println("KxxGraph_:"); show(stdout, "text/plain", KxxGraph_); println()
	println("KyyGraph_:"); show(stdout, "text/plain", KyyGraph_); println()
	println("KzzGraph_:"); show(stdout, "text/plain", KzzGraph_); println()
	
	
	xbonds::Int8 = length(KxxGraph_[KxxGraph_ .> 0])/2  # only need the upper half of KnnGraph_
	ybonds::Int8 = length(KyyGraph_[KyyGraph_ .> 0])/2
	zbonds::Int8 = length(KzzGraph_[KzzGraph_ .> 0])/2
	
	KxxPair_ = Matrix{Int8}(undef, xbonds, 2)  # pairwise non-zero coupling \\
	KyyPair_ = Matrix{Int8}(undef, ybonds, 2)  # 1st and 2nd cols are site indices of pairs 
	KzzPair_ = Matrix{Int8}(undef, zbonds, 2)  
	
	Kxxcoef_ = Vector{Float64}(undef,xbonds)  # pairwise non-zero coupling strength
	Kyycoef_ = Vector{Float64}(undef,ybonds)
	Kzzcoef_ = Vector{Float64}(undef,zbonds)
	
	# extract non-zero x-coupling pairs 
	counter::Int8 = 1
	for i::Int8 in 1:nsite
		for j::Int8 in i:nsite
			if KxxGraph_[i,j]!=0
				KxxPair_[counter,1] = i
				KxxPair_[counter,2] = j
				Kxxcoef_[counter] = KxxGraph_[i,j] 
				counter += 1
			end		
		end
	end
	#println("KxxPair_:"); show(stdout, "text/plain", KxxPair_); println()
	
	# extract non-zero y-coupling pairs 
	counter = 1
	for i::Int8 in 1:nsite
		for j::Int8 in i:nsite
			if KyyGraph_[i,j]!=0
				KyyPair_[counter,1] = i
				KyyPair_[counter,2] = j
				Kyycoef_[counter] = KyyGraph_[i,j] 
				counter += 1
			end	
		end
	end

	#println("KyyPair_:"); show(stdout, "text/plain", KyyPair_); println()
	# extract non-zero z-coupling pairs 
	counter = 1
	for i::Int8 in 1:nsite
		for j::Int8 in i:nsite
			if KzzGraph_[i,j]!=0
				KzzPair_[counter,1] = i
				KzzPair_[counter,2] = j
				Kzzcoef_[counter] = KzzGraph_[i,j] 
				counter += 1
			end	
		end
	end
	#println("KzzPair_:"); show(stdout, "text/plain", KzzPair_); println()

	# ---------------------Build Hamiltonian as Sparse Matrix-------------------
	# ---------------------Build Hamiltonian as Sparse Matrix-------------------
	# ---------------------Build Hamiltonian as Sparse Matrix-------------------
	
	println("[Hamiltonian.jl] Building Hamiltonian as Sparse Matrix...")
	sx::SparseMatrixCSC{Complex{Float64},Int8} = Sx()
	sy::SparseMatrixCSC{Complex{Float64},Int8} = Sy()
	sz::SparseMatrixCSC{Complex{Float64},Int8} = Sz()
	id::SparseMatrixCSC{Complex{Float64},Int8} = sparse(I,2,2)
	
	Ham::SparseMatrixCSC{Complex{Float64},Int64} = spzeros(2^nsite,2^nsite)
	Hamx::SparseMatrixCSC{Complex{Float64},Int64} = spzeros(2^nsite,2^nsite)
	Hamy::SparseMatrixCSC{Complex{Float64},Int64} = spzeros(2^nsite,2^nsite)
	Hamz::SparseMatrixCSC{Complex{Float64},Int64} = spzeros(2^nsite,2^nsite)
	
	
	for i::Int8 in 1:xbonds
		ia::Int8 = KxxPair_[i,1]
		ib::Int8 = KxxPair_[i,2]
		coef::Float64 = Kxxcoef_[i]
		# println("x-x: ia=",ia,", ib=",ib,", coef=",coef)
		ida::SparseMatrixCSC{Complex{Float64},Int64} = sparse(I,2^(ia-1),2^(ia-1))
		idm::SparseMatrixCSC{Complex{Float64},Int64} = sparse(I,2^(ib-ia-1),2^(ib-ia-1))
		idb::SparseMatrixCSC{Complex{Float64},Int64} = sparse(I,2^(nsite-ib),2^(nsite-ib))		
		Hamx += kron(ida, sx, idm, sx, idb)*coef
	end

	for i::Int8 in 1:ybonds
		ia::Int8 = KyyPair_[i,1]
		ib::Int8 = KyyPair_[i,2]
		coef::Float64 = Kyycoef_[i]
		# println("y-y: ia=",ia,", ib=",ib,", coef=",coef)
		ida::SparseMatrixCSC{Complex{Float64},Int64} = sparse(I,2^(ia-1),2^(ia-1))
		idm::SparseMatrixCSC{Complex{Float64},Int64} = sparse(I,2^(ib-ia-1),2^(ib-ia-1))
		idb::SparseMatrixCSC{Complex{Float64},Int64} = sparse(I,2^(nsite-ib),2^(nsite-ib))
		Hamy += kron(ida, sy, idm, sy, idb)*coef	
	end

	for i::Int8 in 1:zbonds
		ia::Int8 = KzzPair_[i,1]
		ib::Int8 = KzzPair_[i,2]
		coef::Float64 = Kzzcoef_[i]
		# println("z-z: ia=",ia,", ib=",ib,", coef=",coef)
		ida::SparseMatrixCSC{Complex{Float64},Int64} = sparse(I,2^(ia-1),2^(ia-1))
		idm::SparseMatrixCSC{Complex{Float64},Int64} = sparse(I,2^(ib-ia-1),2^(ib-ia-1))
		idb::SparseMatrixCSC{Complex{Float64},Int64} = sparse(I,2^(nsite-ib),2^(nsite-ib))
		Hamz += kron(ida, sz, idm, sz, idb)*coef	
	end

	Ham = Hamx + Hamy + Hamz
	
	
	# --------------------------- Add external field -------------------------
	
	for i::Int8 in 1:nsite
		ida::SparseMatrixCSC{Complex{Float64},Int64} = sparse(I,2^(i-1),2^(i-1))
		idb::SparseMatrixCSC{Complex{Float64},Int64} = sparse(I,2^(nsite-i),2^(nsite-i))
		Ham += kron(ida, sx, idb)*Hx
		Ham += kron(ida, sy, idb)*Hy
		Ham += kron(ida, sz, idb)*Hz
	end
	
	
	return Ham


end

#Ham = Hamiltonian(2,1,"PBC")
#eigvals, eigvecs = eigs(Ham, nev = 4, which=:SR)
#println(eigvals)
#show(stdout, "text/plain", nn_); println()

