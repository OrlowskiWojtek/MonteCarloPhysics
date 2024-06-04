using Makie, CairoMakie, DelimitedFiles

include("qvmc.jl")

## 

a_begin = 0.3
a_end   = 1.2
c_begin = -0.7
c_end   = 0.3
Δa = Δc = 0.02

a = [i for i in a_begin:Δa:a_end]
c = [i for i in c_begin:Δc:c_end]

hydrogen_atom = QVMC(
    [a_begin, c_begin],
    (par::Vector{Float64},r::Float64) -> (1 + par[2]*r)*exp(-par[1]*r),
    (par::Vector{Float64},r::Float64) -> (- par[1]^2*par[2]*r^2 + (-par[1]^2 + 4*par[1]*par[2] - 2*par[2])*r + 2*par[1] - 2*par[2] -2)/(2 * par[2] * r^2 + 2 * r),
    1_000_000,
    0.1
)


##

energy = Array{Float64}(undef, length(a), length(c))
variance = Array{Float64}(undef, length(a), length(c))

for i in eachindex(a)
    for j in eachindex(c)
        hydrogen_atom.var_params[1] = a[i]
        hydrogen_atom.var_params[2] = c[j]
        (energy[i,j],variance[i,j]) = calc_energy(hydrogen_atom)
    end
end

save_energy(energy)
save_variance(variance)

ene = plot_energy_variance(energy,variance,a,c)
save("plots/maps.pdf",ene)

## 

hydrogen_atom.var_params[1] = a[36]
hydrogen_atom.var_params[2] = c[36]
r_hist = calc_hist(hydrogen_atom)

h = plot_r_hist(r_hist)
save("plots/hist.pdf", h)