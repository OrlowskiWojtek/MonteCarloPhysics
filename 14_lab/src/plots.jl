function plot_energy_variance(energy::Array{Float64}, variance::Array{Float64}, a::Vector{Float64}, c::Vector{Float64})
    fig = Figure(size=(768,512))

    ax1 = Axis(fig[1,1], title = "Energia lokalna")
    h = heatmap!(ax1,a,c,energy, colormap = :solar)
    Colorbar(fig[2,1],h, vertical = false)

    ax2 = Axis(fig[1,2], title = "σ")
    h = heatmap!(ax2,a,c,variance, colormap = :dense)
    Colorbar(fig[2,2],h, vertical = false)
    hideydecorations!(ax2)

    ax3 = Axis(fig[1,3], title = "log(σ + 10⁻²⁰)")
    h = heatmap!(ax3,a,c,log.(variance .+ 10^(-20)), colormap = :dense)
    Colorbar(fig[2,3],h, vertical = false)
    hideydecorations!(ax3)
    rowsize!(fig.layout, 1, Aspect(3,1))

    return fig
end

function plot_r_hist(r_hist::Vector{Float64})
    xs = LinRange(0,8,200)
    Ψ_exact(r) = r^2* abs(2 * exp(-r))^2


    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "r [a.u.]", ylabel = "r²⋅|Ψ₁₀₀(r)|² ")
    lines!(ax,xs,Ψ_exact, label = "dokładny", color = :blue, linewidth = 3)
    lines!(ax,xs,r_hist, label = "mc: n = 10⁶", color = :red)


    axislegend()

    return fig
end
