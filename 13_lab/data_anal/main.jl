## data analysis in julia
using Makie, CairoMakie, DelimitedFiles, DataFrames
include("utils.jl")
## loading data - rozkłady 

dir = "data/data_evol/"
loc = particle_location(dir)
save("plots/pos_end.pdf", loc)

##

dir = "data/data_nptv/"
iters = [i for i in 1:100:301]
f = plot_nptv(dir, iters)
save("plots/nptv_$(iters[end]).pdf",f)

iters = [i for i in 401:100:701]
f = plot_nptv(dir, iters)
save("plots/nptv_$(iters[end]).pdf",f)

iters = [i for i in 801:200:1801]
f = plot_nptv(dir, iters)
save("plots/nptv_$(iters[end]).pdf",f)
