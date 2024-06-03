function calc_Γ(k::Vector{Float64},x::Vector{Float64})
    return [k[1], k[2], k[3]*x[1]*x[2], k[4]*x[3]]
end

function change_state(m,x)
    if(m == 1) x[1] += 1 end    
    if(m == 2) x[2] += 1 end
    if(m == 3)
        x[1] = x[1] - 1
        x[2] = x[2] - 1
        x[3] = x[3] + 1
    end
    if(m == 4) x[3] = x[3] - 1 end    
end


function perform(k,x_0,tmax,N,pmax)
    x1 = Vector{Vector{Float64}}(undef,pmax)
    x2 = Vector{Vector{Float64}}(undef,pmax)
    x3 = Vector{Vector{Float64}}(undef,pmax)
    time = Vector{Vector{Float64}}(undef,pmax)

    h0 = Array{Float64}(undef, N)
    h1 = zeros(N)
    h2 = zeros(N)
    ncount = Array{Int}(undef, N)
    δt = tmax/N

    ##

    for p = 1:Pmax
        x1[p] = zeros(0)
        x2[p] = zeros(0)
        x3[p] = zeros(0)
        time[p] = zeros(0)

        t = 0
        x = copy(x_0)

        h0 = zeros(N)
        ncount = zeros(N)
        while(t < tmax)
            Γ = calc_Γ(k,x)
            Γ_max = sum(Γ)
            Δt = - 1 / Γ_max * log(rand())
            U2 = rand()
            Γ_norm = Γ ./ Γ_max
            if(U2 <= Γ_norm[1]) m = 1
                elseif(U2 <= Γ_norm[1] + Γ_norm[2]) m = 2
                elseif(U2 <= Γ_norm[1] + Γ_norm[2] + Γ_norm[3]) m = 3   
                else m = 4
            end
            change_state(m,x)
             
            push!(time[p],t)
            push!(x1[p],x[1])
            push!(x2[p],x[2])
            push!(x3[p],x[3])
            
            l = Int(floor(t/δt)+1)
            h0[l] += x[3]
            ncount[l] += 1
            
            t += Δt 
        end
        
        for l = 1:N
            x3_temp = h0[l] / ncount[l]
            h1[l] += x3_temp
            h2[l] += (x3_temp)^2
        end
    end

    x3_avg = Vector{Float64}(undef, N)
    time_avg = Vector{Float64}(undef, N)
    std_err_avg = Vector{Float64}(undef, N)
    for l in 1:N
        μ1 = h1[l]/Pmax
        μ2 = h2[l]/Pmax
        x3_avg[l] = μ1
        std_err_avg[l] = sqrt((μ2 - μ1^2)/Pmax)
        time_avg[l] = (l -1/2)*δt 
    end

    return (time,x1,x2,x3,time_avg,x3_avg,std_err_avg)
end


