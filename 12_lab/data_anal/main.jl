## data analysis in julia
using Makie, CairoMakie, DelimitedFiles, DataFrames
include("utils.jl")
## loading data - rozkłady 

tasks = "task_8"
dir = joinpath("tasks",tasks)
plots_dir = joinpath("plots",tasks)
name = tasks*"pt.pdf"
pt = plot_nptv(dir)
save(joinpath(plots_dir,name),pt)

## to jest trajektoria zerowej cząstki 

file = "xy_1.dat"
xy = readdlm(joinpath(dir,file))

plot(xy[:,1],xy[:,2], lw = 2,color = :black, xlabel ="x", ylabel = "y", label = false, title = "trajektoria pierwszej cząstki", xlims = (0,1), ylims = (0,1))

## położenie i prędkość cząstek XY - x y vx vy rc (rozmiar cząstek) mc (masa cząstek)
dir = "tasks/task_3/data_evol/"
evol = particle_location(dir)
save(joinpath(dir,"evol.pdf"),evol)

##
## histogram prędkosci 

dir = "tasks/task_2"
hist = hist_plot(dir)
save(joinpath(dir,"hists.pdf"),hist)

## iteracji jest git
## histogram prędkość końcowa i początkowa

task = "task_8"
dir = joinpath("tasks",task)
hist = hist_begin_end_plot(dir)
plots_dir = joinpath("plots",task)
name = task*"hists.pdf"
save(joinpath(plots_dir,name),hist)

