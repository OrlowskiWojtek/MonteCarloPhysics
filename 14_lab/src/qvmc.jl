struct QVMC{F1,F2}
    var_params::Vector{Float64}    
    Ψ_trial::F1
    εloc::F2 
    N::Int64
    Δr::Float64
end

function calc_energy(p::QVMC)
    r = 0.1
    loc_en = 0.
    loc_en2 = 0.

    for i in 1:p.N
        r_new  = r + p.Δr * (2*rand() - 1)
        p_acc = min((r_new/r)^2 * abs(p.Ψ_trial(p.var_params, r_new))^2/abs(p.Ψ_trial(p.var_params, r))^2, 1.)
        
        U2 = rand()
        
        if(r_new <= 0)
            # do nothing
        elseif(U2 <= p_acc)
            r = r_new
        elseif(U2 > p_acc)
            # do nothing
        end
        
        local_energy = p.εloc(p.var_params, r)
        loc_en += local_energy
        loc_en2 += local_energy^2
    end

    return (loc_en / p.N, loc_en2 / p.N - (loc_en / p.N)^2)
end

function calc_hist(p::QVMC)
    M = 200
    r_max = 8.
    δ_r = r_max / M
    r_hist = zeros(Int64,M)
    r = 0.1

    for i in 1:p.N
        r_new  = r + p.Δr * (2*rand() - 1)
        p_acc = min((r_new/r)^2 * abs(p.Ψ_trial(p.var_params, r_new))^2/abs(p.Ψ_trial(p.var_params, r))^2, 1.)
        
        U2 = rand()
        
        if(r_new <= 0)
            # do nothing
        elseif(U2 <= p_acc)
            r = r_new
        elseif(U2 > p_acc)
            # do nothing
        end
        k = floor(Int64,r/δ_r) + 1
        if(k > 200) continue end
        r_hist[k] += 1
    end

    return r_hist ./ (p.N *δ_r)
end


