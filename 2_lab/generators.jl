## Plik z funkcjami generatorów o fgp f(x) = 4/5 (1 + x - x³)

function compound(g₁::Float64)
    U₁ = rand()
    U₂ = rand()

    if(U₁ <= g₁)
        return U₂
    else
        return √(1 - √(1 - U₂))
    end
end

## łańcuch Markowa

function p_acc(f::Function, X₁, X₂)
    return min(f(X₂)/f(X₁),1)
end

function markow(N::Int, f::Function, Δ::Float64)
    
    chain = Vector{Float64}(undef, N)
    chain[1] = rand()

    for i in 2:N
        U₁ = rand()
        U₂ = rand()  
        x_new = chain[i-1] + (2*U₁ - 1)*Δ
        if(x_new >= 0 && x_new <= 1 && U₂ < p_acc(f, chain[i-1], x_new))
            chain[i] = x_new
        else
            chain[i] = chain[i-1]
        end
    end
    return chain
end

## metoda eliminacji

function eliminate(f::Function)
    U₁ = rand()
    G₂ = 1.15*rand()

    if(G₂ <= f(U₁))
        return U₁
    else
        eliminate(f)
    end
end

function χ²(F::Function,hist, crit_value)
    
    N = sum(hist.weights)
    Fᵢ = [F(i) for i = hist.edges[begin]]
    
    χ = 0
    
    for i ∈ eachindex(hist.weights)
        p = Fᵢ[i+1] - Fᵢ[i]
        χ += (hist.weights[i] - p * N)^2 / (p * N) 
    end

    hypothesis_failed = χ > crit_value

    return (χ, !hypothesis_failed)
end
