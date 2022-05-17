using LinearAlgebra
include("Parameter.jl")

function Honeycomb(param::Parameter)
	LLX = param.LLX; LLY = param.LLY
	IsPeriodicX = param.IsPeriodicX
	IsPeriodicY = param.IsPeriodicY

	#LLX = 2; LLY = 2
	nsite::Int8 = LLX * LLY * 2
	number1neigh::Int8 = 3
	println("[Lattice.jl] creating Honeycomb lattice... ")

	scalex::Float64 = 2; scaley::Float64 = 4.0/sqrt(3)
	t1 = [1.0 * scalex; 0]
	t2 = [0.5 * scalex; sqrt(3)/2.0 * scaley]

	indx_ = Vector{Int8}(undef, nsite)  # x coordinate of sites
	indy_ = Vector{Int8}(undef, nsite)  # y coordinate of sites
	mesh_ = fill!(Matrix{Int8}(undef, LLX*2+LLY,LLY*2),0)  # lattice on meshgrid
	nn_ = fill!(Matrix{Int8}(undef, nsite, number1neigh),0)  # nearest neighbors
	
	# ---------------- Construct Lattice mesh ------------------------ 
	# ---------------- Construct Lattice mesh ------------------------ 
	# ---------------- Construct Lattice mesh ------------------------ 
	xv::Int8 = 1  # xv is unit cell length in x (vertival) direction
	counter::Int8 = 1	
	for i in 1: LLX
		if i!=1
			xv += t1[1]
		end
				
		xa::Int8 = xv; xb::Int8 = xv + 1
		ya::Int8 = 1; yb::Int8 = 2
		
		for j::Int8 in 1:LLY
			indx_[counter] = xa
			indy_[counter] = ya
			mesh_[xa,ya] = counter
			#println("xa=", xa, " ya=",ya, " counter=",counter)
			counter += 1
			
			indx_[counter] = xb
			indy_[counter] = yb
			mesh_[xb,yb] = counter
			#println("xb=", xb, " yb=",yb, " counter=",counter)
			counter += 1
			
			xa = xv + j*t2[1]; xb = (xv + 1) + j*t2[1]
			ya = 1 + j*t2[2]; yb = 2 + j*t2[2]
										
		end		
	end
	show(stdout, "text/plain", mesh_); println()

	# ------------ Construct matrix of 1st n.n.s --------------------- 
	# ------------ Construct matrix of 1st n.n.s --------------------- 
	# ------------ Construct matrix of 1st n.n.s --------------------- 
	println("[Lattice.jl] Looking for nearest neighbors... ")
	xmax = findmax(indx_)[1]   # findmax returns both max value and its index
	ymax = findmax(indy_)[1]
	
	for i::Int8 in 1:nsite
		ix::Int8 = indx_[i]  # coordinate of n-th site in matrix
		iy::Int8 = indy_[i]
				
		#----------------------------OBC-----------------------------------
		# n.n in x-bond
		jx::Int8 = ix + 1; jy::Int8 = iy + 1  # move 1 step forward in x = (1,1) direction
		if jx <= xmax && jy <= ymax && mesh_[jx,jy] != 0
			j = mesh_[jx,jy]  # site index of n.n. in x direction
			nn_[i,1] = j
			nn_[j,1] = i 	
		end
	
		# n.n in y-bond
		jx = ix + 1; jy = iy - 1  # move 1 step in y = (1,-1) direction		
		if jx <= xmax && jy <= ymax && jy >= 1 && mesh_[jx,jy] != 0
			j = mesh_[jx,jy]  # site index of n.n. in x direction
			nn_[i,2] = j
			nn_[j,2] = i 	
		end	
				
		# n.n in z-bond
		jx = ix; jy = iy + 1  # move 1 step in z = (0,1) direction		
		if jx <= xmax && jy <= ymax && mesh_[jx,jy] != 0
			j = mesh_[jx,jy]  # site index of n.n. in x direction
			nn_[i,3] = j
			nn_[j,3] = i 	
		end
		#----------------------------OBC-----------------------------------
		
		
		#--------------------------Apply PBC-------------------------------	
		if IsPeriodicY == true
			# z-bond
			jx = ix - LLY 
			jy = 1
			if jx >= 1 && iy == ymax && mesh_[jx,jy]!=0
				j = mesh_[jx,jy]
				nn_[i,3] = j
				nn_[j,3] = i
			end
		end
			
		if IsPeriodicX == true
			# y-bond
			jx = ix + 2 * LLX - 1 
			jy = iy + 1
			if jx <= xmax && iy <= ymax && iy%2 !=0 && mesh_[jx,jy]!=0
				j = mesh_[jx,jy]
				nn_[i,2] = j
				nn_[j,2] = i
			end
		end			
	
	end
	
	show(stdout, "text/plain", nn_); println()
	return nsite, mesh_, nn_, indx_, indy_

end

#para = GetParameter("../input.inp")
#nsite, mesh_, nn_, indx_, indy_ = Honeycomb(para)	
#show(stdout, "text/plain", nn_); println()
#show(stdout, "text/plain", mesh_); println()

