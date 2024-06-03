using Plots

include("funcs.jl")

## static parameters

k = [1., 1., 0.001, 0.01]
x_0 = [120., 80., 1.]
tmax = 200.
N = 50

## task1
# variable
Pmax = 1
params = (k = k, x = x_0, tmax = tmax, N = N, Pmax = Pmax)

(time,x1,x2,x3,time_avg,x3_avg,std_err_avg) = perform(params...)

plt = plot(time, x1, label = false, color = :blue, xlabel = "Czas")
plot!([],[], color = :blue, label = "x1")
plot!(time,x2, label = false, color = :red)
plot!([],[], color = :red, label = "x2")
plot!(time,x3, label = false, color = :green)
plot!([],[], color = :green, label = "x3")

display(plt)

savefig("plots/task_1.png")
savefig("plots/svgs/task_1.svg")

## task2

Pmax = 5
params = (k = k, x = x_0, tmax = tmax, N = N, Pmax = Pmax)

(time,x1,x2,x3,time_avg,x3_avg,std_err_avg) = perform(params...)

plt = plot(time, x1, label = false, color = :blue, xlabel = "Czas")
plot!([],[], color = :blue, label = "x1")
plot!(time,x2, label = false, color = :red)
plot!([],[], color = :red, label = "x2")
plot!(time,x3, label = false, color = :green)
plot!([],[], color = :green, label = "x3")

display(plt)

savefig("plots/task_2.png")
savefig("plots/svgs/task_2.svg")

## task_3 

Pmax = 5
params = (k = k, x = x_0, tmax = tmax, N = N, Pmax = Pmax)

(time,x1,x2,x3,time_avg,x3_avg,std_err_avg) = perform(params...)

plt = plot(xlabel = "Czas")
plot!(time,x3, label = false, color = :green)
plot!([],[], color = :green, label = "x3")
plot!(time_avg, x3_avg, ribbon =  std_err_avg, yerror = std_err_avg, label = "x3 uśrednione",title = "Pmax = $Pmax", color = :black, xlabel = "Czas")

display(plt)

savefig("plots/additional_1.png")
savefig("plots/svgs/additional_1.svg")

## some additional plots 
Pmax = 10000
params = (k = k, x = x_0, tmax = tmax, N = N, Pmax = Pmax)

(time,x1,x2,x3,time_avg,x3_avg,std_err_avg) = perform(params...)


hist = Vector{Float64}(undef, length(x3))
for i in eachindex(x3)
    hist[i] = x3[i][end] 
end

histogram(hist, nbins = 20, title = "Histogram wartości x₃ dla t = 200", label = false, xlabel = "x₃")

savefig("plots/additional_2.png")
savefig("plots/svgs/additional_2.svg")

