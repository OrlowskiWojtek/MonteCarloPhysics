
function basic_mc(k::Vector, f::Function, a::Number, b::Number)
    s = 0
    s_sq = 0

    counter = 1
    N = k[end]    

    s_vec = Vector{Float64}(undef,length(k))
    s_sq_vec = Vector{Float64}(undef, length(k))

    M = 10
    table = Array{Float64}(undef,length(k),M)
    hist_data = zeros(M)
    width = (b-a) / 10.

    for i in 1:N
        x = a + (b-a)*rand()

        interval = Int(floor((x-a)/width)) + 1
        hist_data[interval] += 1
        
        val = f(x)
        s += val
        s_sq += val^2
        if(i == k[counter])
            table[counter,:] = hist_data
            s_vec[counter] = s
            s_sq_vec[counter] = s_sq
            counter += 1
        end
    end

    ḡ = 1 ./ k .* (b-a) .* s_vec
    ḡ_2 = 1 ./ k .* (b-a)^2 .* s_sq_vec 

    σ = sqrt.((ḡ_2 .- ḡ .^2) ./ k)

    return (ḡ, σ, table)
end

function basic_mc(N::Int, f::Function, a::Number, b::Number)
    M = 10
    s = 0
    s_sq = 0

    hist_data = zeros(M)
    width = (b-a) / M

    for i in 1:N
        x = a + (b-a)*rand()

        interval = Int(floor((x-a)/width)) + 1
        hist_data[interval] += 1
        
        val = f(x)
        s += val
        s_sq += val^2
    end

    ḡ = 1 / N * (b-a) * s
    ḡ_2 = 1 / N * (b-a)^2 * s_sq 

    σ = sqrt((ḡ_2 - ḡ ^2) / N)

    return (ḡ, σ, hist_data)
end

function sys_mc(N::Int, f::Function, a::Number, b::Number)
    M = 10
    width = (b-a) / 10.
    
    s = zeros(M)
    s_sq = zeros(M)    
    
    Nₘ = Int(N/M)
    x_vec = LinRange(a,b,M+1)
    hist_data = zeros(M)
    
    for j in 1:M
        for i in 1:Nₘ
            x = x_vec[j] + (x_vec[j+1]-x_vec[j])*rand()        
            val = f(x)
            s[j] += (b-a)*val
            s_sq[j] += ((b-a)*val)^2
            interval = Int(floor((x-a)/width)) + 1
            hist_data[interval] += 1
        end
    end
    
    ḡ = s ./ Nₘ
    ḡ_2 = s_sq ./ Nₘ
    
    σ_m = ḡ_2 .- ḡ .^2

    C = sum(ḡ)/M
    pₘ = 1/M
    σ = sqrt(pₘ^2 / Nₘ * sum(σ_m) )
    

    return (C, σ, hist_data)
end

function layer_mc(N::Int, f::Function, a::Number, b::Number, estimator::Function, n_estimate)
    M = 10
    width = (b-a) / M
    
    s = zeros(M)
    s_sq = zeros(M)    
    
    x_vec = LinRange(a,b,M+1)
    hist_data = zeros(M)
    
    σ_vec = Vector{Float64}(undef,M)

    for i in 1:M
        (μ,σ_vec[i],table) = estimator(n_estimate,f,x_vec[i],x_vec[i+1])
    end  

    σ_sum = sum(σ_vec)
    N_vec = round.(σ_vec ./ σ_sum * N)
    
    for j in 1:M 
        Nₘ = N_vec[j]
        for i in 1:Nₘ
            x = x_vec[j] + (x_vec[j+1]-x_vec[j])*rand()        
            val = f(x)
            s[j] += (b-a)*val
            s_sq[j] += ((b-a)*val)^2
            interval = Int(floor((x-a)/width)) + 1
            hist_data[interval] += 1
        end
    end
    
    ḡ = Vector{Float64}(undef, length(s))
    ḡ_2 = Vector{Float64}(undef, length(s))
    for i in eachindex(s)
        if(N_vec[i] == 0)
            ḡ[i] = 0
            ḡ_2[i] = 0
            continue
        end
        ḡ[i] = s[i] / N_vec[i]
        ḡ_2[i] = s_sq[i] / N_vec[i]
    end

    σ²_m = ḡ_2 .- ḡ .^2

    C = sum(ḡ)/M
    pₘ = 1/M
    
    s_val = 0
    for i in eachindex(N_vec)
        if(N_vec[i] == 0)
            continue
        end
        s_val += σ²_m[i] / N_vec[i]
    end

    σ = sqrt(pₘ^2 * s_val)
    
    return (C, σ, hist_data)
end
