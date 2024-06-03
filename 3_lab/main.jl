using Plots

include("generators.jl")

n = 10_000

gr(size =(800,800))
## zadanie 1

(x,y) = gen_n_box_muller(n,false)
r = gen_n_box_muller(n,true)

plt = scatter(x,y, markersize = 1, color = :grey, title = "Rozkład normalny, metoda Boxa-Mullera", xticks = [-4,-2,0,2,4], label = false)
display(plt)
savefig("galeria/box_muller.png")

## zadanie 2

(x,y) = gen_n_on_circle(n)
plt = scatter(x,y,markersize = 1, color = :grey, title ="Rozkład jednorodny na kole", label = false)
display(plt)

savefig("galeria/na_kole.png")

(x,y) = gen_n_in_circle(n)
plt = scatter(x,y,markersize = 1, color = :grey, title ="Rozkład jednorodny w kole", label = false)
display(plt)

savefig("galeria/w_kole.png")

## zadanie 3

(x,y) = gen_n_in_circle(n)
(x,y) = move_to_ellipse(π/4,1,0.2,x,y)
plt = scatter(x,y,markersize = 1, color = :grey, title ="Rozkład jednorodny w elipsie", label = false)
display(plt)

savefig("galeria/elipsa_jednorodny.png")

cov_coeff(x,y)

## zadanie 4

α = π/4
b₁ = 1
b₂ = 0.2

(x,y) = gen_n_box_muller(n,false)
(x,y) = move_to_ellipse(α,b₁,b₂,x,y)
plt = scatter(x,y,markersize = 1, color = :grey, title ="Rozkład normalny po dokonanej transformacji", label = false)
display(plt)

savefig("galeria/elipsa_normalny.png")

println("Współczynnik kowariancji z wygenerowanego rozkładu: $(cov_coeff(x,y))")

A = gen_trans_matrix(α,b₁,b₂)
println("Współczynnik kowariancji z obliczeń: $(cov_coeff(A*A'))")
