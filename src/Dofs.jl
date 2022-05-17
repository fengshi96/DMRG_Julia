using SparseArrays

	
function Sx()
	I = Int8[1,2]; J = Int8[2,1]; V = Complex{Float64}[0.5,0.5]
	Sx = sparse(I,J,V)
	return Sx

end


function Sy()
	I = Int8[1,2]; J = Int8[2,1]; V = Complex{Float64}[-0.5im,0.5im]
	Sy = sparse(I,J,V)
	return Sy

end

function Sz()
	I = Int8[1,2]; J = Int8[1,2]; V = Complex{Float64}[0.5,-0.5]
	Sz = sparse(I,J,V)
	return Sz

end

#show(stdout, "text/plain", Sx()); println()
#stemp = kron(Sx(),Sx())
#kron(Sx(),Sx())
#println(kron(  Sx(),  Sx()  ) )












