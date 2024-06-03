using Plots

include("gen.jl")

gr(size = (1024,768))

## zadanie 2
Ra = 2. 
Rb = 2. *sqrt(2)

circle_1 = (Ra,Ra + Rb, 0.0)
circle_2 = (sqrt(2)*Ra,0.0,0.0)
n = 10_000

gen_circles_plot(circle_1,circle_2,n)
savefig("gallery/2_test.png")

##
save_at = [10^k for k in 2:7]
@time (μ,σ) = calc_common_area(circle_1,circle_2, save_at)

plot(save_at,μ,yerror = σ, xaxis = :log, xlabel = "n", ylabel = "Wspólne pole", lw = 2, color = :red, label = gen_circle_label(circle_1)*"\n"*gen_circle_label(circle_2))
plot!(save_at, μ .- σ, fillrange = μ .+ σ, fc=:red, fillalpha = 0.35, label = false)

## zad 3
# a)
Ra = 2. 
Rb = 2. *sqrt(2)
save_at = [10^k for k in 2:6]

circle_1 = (Ra,Rb + 0.5 * Ra,0.0)
circle_2 = (Rb, 0., 0.)

@time (μ,σ) = calc_common_area(circle_1,circle_2, save_at)
plt = plot(save_at,μ,yerror = σ, xaxis = :log, xlabel = "n", ylabel = "Wspólne pole", lw = 2, color = :red, label = "a)", markerstrokecolor = :red, legendfontsize = 20)
plot!(save_at, μ .- σ, fillrange = μ .+ σ, fc=:red, fillalpha = 0.35, label = false)

# c) 

@time (μ,σ) = calc_common_area(circle_2,circle_1, save_at)
plot!(save_at, μ, yerror = σ, xaxis = :log, xlabel = "n", ylabel = "Wspólne pole", lw = 2, color = :blue, label = "c)", markerstrokecolor = :blue)
plot!(save_at, μ .- σ,  fillrange = μ .+ σ, fc=:blue, fillalpha = 0.35, label = false)
savefig("gallery/a_c.png")

gen_circles_plot(circle_1,circle_2,save_at[end-2])
savefig("gallery/a_c_circ.png")
## b)

Ra = 2. 
Rb = 2. *sqrt(2)
save_at = [10^k for k in 2:6]

circle_1 = (Ra, 0., 0.)
circle_2 = (Rb, 0., 0.)

@time (μ,σ) = calc_common_area(circle_1,circle_2, save_at)
plt = plot(save_at,μ,yerror = σ, xaxis = :log, xlabel = "n", ylabel = "Wspólne pole", lw = 2, color = :red, label = "b)", markerstrokecolor = :red, legendfontsize = 20)
plot!(save_at, μ .- σ, fillrange = μ .+ σ, fc=:red, fillalpha = 0.35, label = false)

# d) 

@time (μ,σ) = calc_common_area(circle_2,circle_1, save_at)
plot!(save_at, μ, yerror = σ, xaxis = :log, xlabel = "n", ylabel = "Wspólne pole", lw = 2, color = :blue, label = "d)", markerstrokecolor = :blue)
plot!(save_at, μ .- σ,  fillrange = μ .+ σ, fc=:blue, fillalpha = 0.35, label = false)

savefig("gallery/b_d.png")

gen_circles_plot(circle_1,circle_2,save_at[end-2])
savefig("gallery/b_d_circ.png")

