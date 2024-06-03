using Makie
using CairoMakie
using OffsetArrays

include("mc_solve.jl")
include("utils.jl")
include("poisson_solve.jl")

## parameterers

nx = 30
ny = 30
Δ = 0.05
ρmax = 10.
V_boundary = (Vl = 1.,Vb =-1., Vt =-1.)
p = (ε=1, ρmax = ρmax, Δ = Δ)

##
Aρ = fill_rho(nx,ny,Δ,ρmax) ;
V_rel = init_V(V_boundary...,nx,ny,Δ);
@time solve(V_rel,Aρ,p)
f = plot_V(V_rel, false)
save("plots/relaxation.svg",f)
f = plot_V(V_rel, true)
save("plots/relaxation_smooth.svg",f)

## monte carlo calucations
N = 100
n = 100
block = false
Δ = 0.1

mc = init_mc(nx,ny,N,n,block,Δ,ρmax,V_boundary);
@time solve(mc)
f = plot_V(mc.V, false)

do_all_plots(mc,V_rel,"100_100_v2")

##
N = 100
n = 100
block = true 
mc = init_mc(nx,ny,N,n,block,Δ,ρmax,V_boundary);
@time solve(mc)

f = plot_V(mc.V, false)

do_all_plots(mc,V_rel,"100_100_block_v2")

##
N = 300
n = 300
block = false
mc = init_mc(nx,ny,N,n,block,Δ,ρmax,V_boundary);
@time solve(mc)
f = plot_V(mc.V, false)

do_all_plots(mc,V_rel,"test")

