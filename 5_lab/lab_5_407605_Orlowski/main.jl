using Plots, Printf

include("mc.jl")

f1(x) = 1. + tanh(x) 
f2(x) = 1. / (1. + x^2) 
f3(x) = cos(π*x)^10

f = [f1,f2,f3]
bound_a = [-3, 0, 0]
bound_b = [3, 10, 1]
exact_values = [6.0, atan(10) - atan(0), 0.24609375]
funcs_string = ["1 + tanh(x)", "1 / (1 + x^2)", "cos(π x)¹⁰"]
## 

k = [10^k for k in 2:5]
(μ,σ, table) = basic_mc(k,f1, -3,  3)

plot(k,μ, yerror = σ, xaxis = :log, xlabel = "n", ylabel = "Wartość całki",title = "Metoda podstawowa",lw = 2, color = :red, label = false)
plot!(k, μ .- σ, fillrange = μ .+ σ, fc=:red, fillalpha = 0.35, label = false)
hline!([6.0], label = "wartość dokładna", fc = :blue, lw = 1)

savefig("gallery/metoda_podstawowa_f1.png")

## 

k⃗ = [10^i for i in 2:5]
μ = Vector{Float64}(undef,length(k⃗))
σ = Vector{Float64}(undef,length(k⃗))
table = Array{Float64}(undef, length(k⃗), 10)

for i in eachindex(k⃗)
    (μ[i], σ[i] ,table[i,:]) = sys_mc(k⃗[i],f1, -3,  3)
end 

plot(k⃗,μ, yerror = σ, xaxis = :log, xlabel = "n", ylabel = "Wartość całki",title = "Metoda warstwowa",lw = 2, color = :red, label = false)
plot!(k⃗, μ .- σ, fillrange = μ .+ σ, fc=:red, fillalpha = 0.35, label = false)
hline!([6.0], label = "wartość dokładna", fc = :blue, lw = 1)

savefig("gallery/metoda_warstwowa_f1.png")

## 

k⃗ = [10^i for i in 2:5]
μ = Vector{Float64}(undef,length(k⃗))
σ = Vector{Float64}(undef,length(k⃗))
table = Array{Float64}(undef, length(k⃗), 10)

for i in eachindex(k⃗)
    (μ[i], σ[i] ,table[i,:]) = layer_mc(k⃗[i],f1, -3,  3, basic_mc,100)
end 

plot(k⃗,μ, yerror = σ, xaxis = :log, xlabel = "n", ylabel = "Wartość całki",title = "Metoda warstwowa, optymalna",lw = 2, color = :red, label = false)
plot!(k⃗, μ .- σ, fillrange = μ .+ σ, fc=:red, fillalpha = 0.35, label = false)
hline!([6.0], label = "wartość dokładna", fc = :blue, lw = 1)

savefig("gallery/metoda_warstwowa_opt_f1.png")

## do all (wlasciwie zad 3)
k⃗ = [10^i for i in 2:5]

for i in eachindex(f)
    for k in eachindex(k⃗)
        println("N = $(k⃗[k]); f(x) = $(funcs_string[i]); Dokładna wartość = $(exact_values[i])")        
        (μ₁,σ₁,table₁) =  basic_mc(k⃗[k],f[i],bound_a[i], bound_b[i])
        (μ₂,σ₂,table₂) =  sys_mc(k⃗[k],f[i],bound_a[i], bound_b[i])
        if(k < 10^3)
            (μ₃,σ₃,table₃) =  layer_mc(k⃗[k],f[i],bound_a[i], bound_b[i], basic_mc, 100)
        else
            (μ₃,σ₃,table₃) =  layer_mc(k⃗[k],f[i],bound_a[i], bound_b[i], basic_mc, 1000)
        end
        println(@sprintf "Metoda podstawowa   :  C = %.5f ; sigma = %.5f; \\\\" μ₁ σ₁)
        println(@sprintf "Metoda systematycza :  C = %.5f ; sigma = %.5f; \\\\" μ₂ σ₂)
        println(@sprintf "Metoda warstwowa    :  C = %.5f ; sigma = %.5f; \\\\" μ₃ σ₃)
        println()
    end
end

## plot few histograms

k⃗ = [10^i for i in 2:2]

for i in eachindex(f)
    for k in eachindex(k⃗)
        (μ₃,σ₃,table₃) =  layer_mc(k⃗[k],f[i],bound_a[i], bound_b[i], basic_mc, 100)
        (μ₂,σ₂,table₂) =  sys_mc(k⃗[k],f[i],bound_a[i], bound_b[i])
        (μ₁,σ₁,table₁) =  basic_mc(k⃗[k],f[i],bound_a[i], bound_b[i])
        range = LinRange(bound_a[i],bound_b[i],11)[begin:end-1] .+ 0.5*(LinRange(bound_a[i],bound_b[i],11)[2] - LinRange(bound_a[i],bound_b[i],11)[1])
        dense_range = LinRange(bound_a[i],bound_b[i], 200)

        width = range[2]-range[1]
        
        plt = plot(bar(range,table₃, label = "Metoda warstwowa",fillalpha = 0.3,bar_width = width ),title = "N = $(k⃗[k])" , xlabel = "Zakres", ylabel = "Ilość zliczeń")
        plot!(bar!(range,table₁, label = "Metoda podstawowa", bar_width = width, fillalpha = 0.3))
        plot!(bar!(range,table₂, label = "Metoda systematyczna", bar_width = width, fillalpha = 0.3), legend = :top)
        plot!(twinx(),dense_range, f[i],lw = 2, lc = :red, label = funcs_string[i], ylabel = "f(x) = "*funcs_string[i], legend = :bottomright )
        display(plt)
        savefig(plt,"gallery/f$(i)_low.png")
    end
end

