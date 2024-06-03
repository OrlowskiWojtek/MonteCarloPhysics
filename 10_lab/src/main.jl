## testing pseudoalgorithm before wrapping everything in function
using Makie, CairoMakie
using ProgressMeter

include("u_exact_wrap.jl")
include("utils.jl")
include("mc.jl")

##
L = 0.25 * 10^(-6) # μH
C = 100. * 10^(-12) # pF
R = 12.5 # Ω
G = 0.5  * 10^(-3) # mS
l = 2.   # m
Rl= 12.5 # Ω
Rg= 75. # Ω

ν = 1. * 10^(9) # GHz
t₀= 7.5  * 10^(-9) # ns
σ = 0.75  * 10^(-9) # ns

params = (L = L, C = C, R = R, G = G, l = l, Rl = Rl, Rg = Rg, ν = ν, t₀ = t₀, σ = σ)

## 
x = Vector(LinRange{Float64}(0,2,200))
nchains = [10^5]#[10^3, 10^4, 10^5]

for time in [10.,15,25,35,50]
    time *= 10^(-9)
    uexact = exact(x,time,params)
    sol = Vector{Vector{Float64}}(undef,3)
    for j in eachindex(nchains)
        sol[j] = Vector{Float64}(undef,length(x))

        println("times for nchains = $(nchains[j])")
        @time for i in eachindex(x)
            (b1,b2,sol[j][i]) = mc_solve(nchains[j],x[i],time,params)
        end
    end
    #fig = do_whole_plot(x,sol,uexact,"czas = $(Int(round(time*10^(9)))) ns")
    #save("plots/whole_plots/time_$(Int(round(time*10^(9)))).png",fig)
end
## 
calc_times = [3.89, 4.97, 7.126, 9.4085, 12.834]
times = [10.,15,25,35,50]

f = Figure()
ax = Axis(f[1,1], xlabel = "czas propagacji fali [ns]", ylabel = "czas obliczeń [s]")
scatter!(times, calc_times, color = :red, label = "czasy obliczeń")

A = [ones(length(times)) times]
sol = A' * A \ (A' * calc_times) 
lin(x) = sol[2]*x + sol[1]
lines!(times,lin, label = "dopasowanie liniowe; a = $(round(sol[2],digits = 2))")

axislegend(position = :lt)

display(f)
save("plots/times.svg",f)

##
p = params
N = 111000
λ =  0.5 * (p.R/p.L - p.G/p.C)
μ = p.G / p.C
expval_s = 1/(λ + ν)

a = N / expval_s
