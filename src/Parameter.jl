struct Parameter
	LLX::Int8; LLY::Int8
	IsPeriodicX::Bool
	IsPeriodicY::Bool
	Model::String
	Kxx::Float64 
	Kyy::Float64
	Kzz::Float64
	Hx::Float64
	Hy::Float64
	Hz::Float64
	Nstates::Int64	
end


function GetParameter(path)
	file = open(path,"r")
	input = read(file, String)
	input = split(input)

	LLX::Int8 = 1; LLY::Int8 = 1
	IsPeriodicX::Bool = true
	IsPeriodicY::Bool = true
	Model::String = "Kitaev"
	Kxx::Float64 = 0 
	Kyy::Float64 = 0 
	Kzz::Float64 = 0
	Hx::Float64 = 0 
	Hy::Float64 = 0 
	Hz::Float64 = 0
	Nstates::Int64 = 1
	
	for i in 1:length(input)
		line = split(input[i],"=")

		if line[1]=="LLX"
			LLX = parse(Int8, line[2])
		elseif line[1]=="LLY"
			LLY = parse(Int8, line[2])
		elseif line[1]=="IsPeriodicX"
			IsPeriodicX = parse(Bool, line[2])
		elseif line[1]=="IsPeriodicY"
			IsPeriodicY = parse(Bool, line[2])
		elseif line[1]=="Model"
			Model = line[2]
		elseif line[1]=="Kxx"
			Kxx = parse(Float64, line[2])
		elseif line[1]=="Kyy"
			Kyy = parse(Float64, line[2])
		elseif line[1]=="Kzz"
			Kzz = parse(Float64, line[2])	
		elseif line[1]=="Bxx"
			Hx = parse(Float64, line[2])	
		elseif line[1]=="Byy"
			Hy = parse(Float64, line[2])
		elseif line[1]=="Bzz"
			Hz = parse(Float64, line[2])	
		elseif line[1]=="Nstates"
			Nstates = parse(Int64, line[2])		
		end		
	end
	
	


	para = Parameter(LLX, LLY, IsPeriodicX, IsPeriodicY, Model, Kxx, Kyy, Kzz, Hx, Hy, Hz, Nstates)
	return para













end


#para = GetParameter("../input.inp")
#show(stdout, "text/plain", para.Nstates); println()



