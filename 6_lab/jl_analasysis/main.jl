using Plots, LsqFit

std_error_xx = [0.00295, 0.0015, 0.00037, 0.0000757]
N = [100., 1000., 10000., 100000.]


p0 = [1.]
m(x, p) = p[1] ./ sqrt.(x)


fit = curve_fit(m,N,std_error_xx,p0)

fit.param
f(x) = fit.param[begin] / sqrt(x)

fit_n = [i for i in N[begin]:100:N[end]]

plt = plot(fit_n,f, label = "f(x) = 0.0311/sqrt(x)", xlabel = "Ilość próbek", ylabel = "Odchyelenie standardowe Dxx")
scatter!(N,std_error_xx, label = "Dane")
savefig("std_dxx.png")

##

std_error_yy = [0.00380, 0.00062, 0.000259, 0.000094718]
N = [100., 1000., 10000., 100000.]


p0 = [1.]
m(x, p) = p[1] ./ sqrt.(x)


fit = curve_fit(m,N,std_error_yy,p0)

fit.param
f(x) = fit.param[begin] / sqrt(x)
fit.param[begin]

fit_n = [i for i in N[begin]:100:N[end]]

plt = plot(fit_n,f, label = "f(x) = 0.03623/sqrt(x)", xlabel = "Ilość próbek", ylabel = "Odchyelenie standardowe Dyy")
scatter!(N,std_error_yy, label = "Dane")
savefig("std_dyy.png")


