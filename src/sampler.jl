function projector(i::Int, axis::String, ind::Index)
    # i = 1 :: SpinUp;  i = 0 :: SpinDown
    if axis == "z"
        if i == 1
#             return (op("Sz", ind) * 2 + op("Id", ind)) / 2
            data = zeros(2)
            data[1] = 1
            return ITensor(data, ind)
        elseif i == 0
#             return (-op("Sz", ind) * 2 + op("Id", ind)) / 2
            data = zeros(2)
            data[2] = 1
            return ITensor(data, ind)
        else
            error("invalid spin dir")
        end
    elseif axis == "x"
        if i == 1
            data = zeros(2)
            data[1] = 1/sqrt(2)
            data[2] = 1/sqrt(2)
            return ITensor(data, ind)
#             return (op("Sx", ind) * 2 + op("Id", ind)) / 2
        elseif i == 0
            data = zeros(2)
            data[1] = 1/sqrt(2)
            data[2] = -1/sqrt(2)
            return ITensor(data, ind)
#             return (-op("Sx", ind) * 2 + op("Id", ind)) / 2
        else
            error("invalid spin dir")
        end
    elseif axis == "y"
        if i == 1
            data = zeros(Complex{Float64}, 2)
            data[1] = 1/sqrt(2)
            data[2] = Complex(0, 1)/sqrt(2)
            return ITensor(data, ind)
#             return (op("Sy", ind) * 2 + op("Id", ind)) / 2
        elseif i == 0
            data = zeros(Complex{Float64}, 2)
            data[1] = 1/sqrt(2)
            data[2] = Complex(0, -1)/sqrt(2)
            return ITensor(data, ind)
#             return (-op("Sy", ind) * 2 + op("Id", ind)) / 2
        else
            error("invalid spin dir")
        end
    else
        error("invalid axis")
    end
end


function snapShot(sites::Vector{Int64}, ψ::MPS, axes::Vector{String}, rngIndx=1234)
    rng = MersenneTwister(rngIndx)
    sample = Vector{Int8}(undef, length(sites))  # 1 = spin up; 0 = spin down
    probs = Vector{Float64}(undef, length(sites))  # stores probability of a local spin up (in any axis)
    clamped = nothing  # to be attached to MPS after each measurement
    clamps = Vector{String}(undef, length(sites))
    
    orthogonalize!(ψ, sites[1])
    for i = 1:length(sites)
        ψi::ITensor = ψ[sites[i]]
        isnothing(clamped) ? nothing : ψi *= clamped

        ψi1::ITensor = projector(1, axes[i], findindex(ψi, "Site")) * ψi
        prob1::Float64 = real(scalar(ψi1 * dag(ψi1)))  # prob of spin up
        ψi0::ITensor = projector(0, axes[i], findindex(ψi, "Site")) * ψi
        prob0::Float64 = real(scalar(ψi0 * dag(ψi0)))  # prob of spin down
        println(prob1, " ", prob0)
        @assert prob0 + prob1 ≈ 1
        
        if prob1 > rand(rng)
            sample[i] = 1
            clamped = ψi1 / sqrt(prob1)
            clamps[i] =  "P" * axes[i] * "↑ @" * string(i)
        else
            sample[i] = 0
            clamped = ψi0 / sqrt(prob0)
            clamps[i] = "P" * axes[i] * "↓ @" * string(i)
        end
        
        probs[i] = prob1
    end

    return sample, probs, clamps
end



