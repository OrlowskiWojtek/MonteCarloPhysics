using Plots, LaTeXStrings
using StatsBase

include("generators.jl")

N = 1_000_000
f(x) = 4/5 * (1 + x - x^3)
F(x) = 4/5 * (x + x^2/2. - x^4/4.)
## rozkład złożony

g₁ = 4/5

cmp = Vector{Float64}(undef, N)
@time for i ∈ eachindex(cmp)
    cmp[i] = compound(g₁)
end

plt = histogram(cmp,bins = range(0,1,11), normalize =:pdf, color = :violet, label = L"Metoda rozkładu złożonego, $N=10^6$ próbek", ylim = (0.8, 1.15))
plot!([0:0.01:1],f, lw = 3, color = :blue, label = L"f(x) = \frac{4}{5} \left(1 + x - x^3\right)", legend=:bottom)

#savefig("gallery/compound.svg")
display(plt)
## łańcuch Markowa (a)

@time mrk_005 = markow(N,f,0.05)
histogram(mrk_005,bins = range(0,1,11), normalize =:pdf, color = :violet, label = L"Metoda łańcucha Markowa, $N=10^6$ próbek", ylim = (0.8, 1.15))
plot!([0:0.01:1],f, lw = 3, color = :blue, label = L"f(x) = \frac{4}{5} \left(1 + x - x^3\right)", legend=:bottom)

#savefig("gallery/markov_005.svg")

# (b)

@time mrk_05 = markow(N,f,0.5)
histogram(mrk_05,bins = range(0,1,11), normalize =:pdf, color = :violet, label = L"Metoda łańcucha Markowa, $N=10^6$ próbek", ylim = (0.8, 1.15))
plot!([0:0.01:1],f, lw = 3, color = :blue, label = L"f(x) = \frac{4}{5} \left(1 + x - x^3\right)", legend=:bottom)

#savefig("gallery/markov_05.svg")

# (a) vs (b)
plt = histogram(mrk_005,bins = range(0,1,11),alpha = 0.4, normalize =:pdf, color = :violet , label = L"Metoda łańcucha Markowa, $N=10^6$ próbek, $\Delta = 0.05$", ylim = (0.8, 1.15))
histogram!(mrk_05,bins = range(0,1,11),alpha = 0.4, normalize =:pdf, color = :blue, label = L"Metoda łańcucha Markowa, $N=10^6$ próbek,  $\Delta = 0.5$")
plot!([0:0.01:1],f, lw = 3, color = :blue, label = L"f(x) = \frac{4}{5} \left(1 + x - x^3\right)", legend=:bottom)

#savefig("gallery/markov_05_005.svg")
display(plt)
##

el = Vector{Float64}(undef, N)
@time for i in 1:N
    el[i] = eliminate(f)
end

plt = histogram(el,bins = range(0,1,11), normalize =:pdf, color = :violet, label = L"Metoda eliminacji, $N=10^6$ próbek", ylim = (0.8, 1.15))
plot!([0:0.01:1],f, lw = 3, color = :blue, label = L"f(x) = \frac{4}{5} \left(1 + x - x^3\right)", legend=:bottom)

#savefig("gallery/elimination.svg")
display(plt)
## Χ² analiza dla poziomu ufności - 0.95 (α = 0.05)

cmp_hist = fit(Histogram, cmp, nbins = 10)
mrk_005_hist = fit(Histogram, mrk_005, nbins = 10)
mrk_05_hist = fit(Histogram, mrk_05, nbins = 10)
el_hist = fit(Histogram, el, nbins = 10)

println(χ²(F,cmp_hist,16.919))
println(χ²(F,mrk_005_hist,16.919))
println(χ²(F,mrk_05_hist,16.919))
println(χ²(F,el_hist,16.919))


## 

for j in 1:5
    mrk_005 = markow(N,f,0.05)
    mrk_005_hist = fit(Histogram, mrk_005, nbins = 10)
    println(χ²(F,mrk_005_hist,16.919))
end