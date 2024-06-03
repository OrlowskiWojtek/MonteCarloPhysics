using Plots

include("funcs.jl")

N = 10_000_000

p = [0.1, 0.5, 0.9]

exp_val = Array{Float64}(undef, 3, 6)
var = Array{Float64}(undef, 3, 6)
sum_val = zeros(3)
sum_squares = zeros(3)

## zadanie 1.3

@time for pᵢ ∈ eachindex(p)
    war = 100 
    for i ∈ 1:N
        rand_val = bernoulli(p[pᵢ])
        sum_val[pᵢ] += rand_val
        sum_squares[pᵢ] += rand_val^2
        if(i == var)
            exp_val[pᵢ,i] = sum_val[pᵢ]/i
            var[pᵢ,i] = sum_squares[pᵢ]/i - exp_val[pᵢ,i]^2 
            war = war * 10   
        end
    end
end

## 

@time for pᵢ ∈ eachindex(p)
    for i ∈ 1:N
        rand_val = bernoulli(p[pᵢ])
        sum_val[pᵢ] += rand_val
        sum_squares[pᵢ] += rand_val^2
        exp_val[pᵢ,i] = sum_val[pᵢ]/i
        var[pᵢ,i] = sum_squares[pᵢ]/i - exp_val[pᵢ,i]^2 
    end
end

## wykonywanie wykresów

make_exp_val_plot(p,exp_val)
make_var_plot(p,var)
