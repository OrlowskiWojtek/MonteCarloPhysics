
## Histogram plot recipe

function hist_plot(dir::String)
    file_end  = "hist_end.dat"
    file_500  = "hist500.dat"
    file_1000 = "hist1000.dat"
    h_end = readdlm(joinpath(dir, file_end))
    h_500 = readdlm(joinpath(dir, file_500))
    h_1000= readdlm(joinpath(dir, file_1000))

    fig = Figure()
    
    d = h_500
    ax_500 = Axis(fig[1,2])
    lines!(ax_500,d[:,1],d[:,3],label = "dokładne", linewidth = 2, color = :black)
    lines!(ax_500,d[:,1],d[:,2], label = "mc")
    info_500 = Label(fig[1,3], text = "It = 500", rotation = pi/2, halign = :left)
    info_500.tellheight = false
    axislegend()


    d = h_1000
    ax_1000 = Axis(fig[2,2])
    lines!(ax_1000,d[:,1],d[:,3],label = "dokładne", linewidth = 2, color = :black)
    lines!(ax_1000,d[:,1],d[:,2], label = "mc")
    info_1000 = Label(fig[2,3], text = "It = 1000", rotation = pi/2, halign = :left)
    info_1000.tellheight = false


    d = h_end
    ax_end = Axis(fig[3,2])
    lines!(ax_end,d[:,1],d[:,3],label = "dokładne",linewidth = 2, color = :black)
    lines!(ax_end,d[:,1],d[:,2], label = "mc")
    info_end = Label(fig[3,3], text = "It = 20000", rotation = pi/2, halign = :left)
    info_end.tellheight = false

    ylabel = Label(fig[1:3,1],text ="rozkład cząstek", font = :bold, rotation = pi/2)
    xlabel = Label(fig[4,1:3],text="prędkość",font = :bold)
    
    return fig
end

function particle_location(dir::String)
    
    fig = Figure()
    ax1 = Axis(fig[1,1], title = "początek")

    skip = 5

    file = "data_evolution_it1"
    rv = readdlm(joinpath(dir,file), Float64, skipblanks = true)[begin:end-1,:]
    scatter!(ax1, rv[begin:skip:end,1],rv[begin:skip:end,2], title = "It = 1", markersize = 1.0, color = :black)
    xlims!(0,1)
    ylims!(0,1)

    ax2 = Axis(fig[1,2], title = "iteracja 20")
    file = "data_evolution_it21"
    rv = readdlm(joinpath(dir,file), Float64, skipblanks = true)[begin:end-1,:]
    scatter!(ax2, rv[begin:skip:end,1],rv[begin:skip:end,2], title = "It = 1", markersize = 1.0, color = :black)
    xlims!(0,1)
    ylims!(0,1)
    
    ax3 = Axis(fig[1,3], title = "iteracja 40" )
    file = "data_evolution_it41"
    rv = readdlm(joinpath(dir,file), Float64, skipblanks = true)[begin:end-1,:]
    scatter!(ax3, rv[begin:skip:end,1],rv[begin:skip:end,2], title = "It = 1", markersize = 1.0, color = :black)
    xlims!(0,1)
    ylims!(0,1)
    
    ax4 = Axis(fig[2,1], title = "iteracja 60")
    file = "data_evolution_it61"
    rv = readdlm(joinpath(dir,file), Float64, skipblanks = true)[begin:end-1,:]
    scatter!(ax4, rv[begin:skip:end,1],rv[begin:skip:end,2], title = "It = 1", markersize = 1.0, color = :black)
    xlims!(0,1)
    ylims!(0,1)

    ax5 = Axis(fig[2,2], title = "iteracja 80")
    file = "data_evolution_it81"
    rv = readdlm(joinpath(dir,file), Float64, skipblanks = true)[begin:end-1,:]
    scatter!(ax5, rv[begin:skip:end,1],rv[begin:skip:end,2], title = "It = 1", markersize = 1.0, color = :black)
    xlims!(0,1)
    ylims!(0,1)

    ax6 = Axis(fig[2,3], title = "koniec")
    file = "../rv.dat"
    rv = readdlm(joinpath(dir,file), Float64, skipblanks = true)[begin:end-1,:]
    scatter!(ax6, rv[begin:skip:end,1],rv[begin:skip:end,2], title = "It = 1", markersize = 1.0, color = :black)
    xlims!(0,1)
    ylims!(0,1)

    return fig
end

function plot_nptv(dir::String)
    fig = Figure()
    wyniki = "wyniki"

    its = [100, 1000, 10000, 20000]

    for i in eachindex(its)
        xpos = Int(floor(i / 3)) + 1
        ypos = Int(i%2 == 0) + 1 

        axl = Axis(fig[xpos, ypos], yticklabelcolor = :blue, title = "It = $(its[i])")
        axr = Axis(fig[xpos, ypos], yticklabelcolor = :red , yaxisposition = :right)

        nptv = readdlm(joinpath(dir,wyniki,"nptv_$(its[i]).dat"))
        lines!(axl, nptv[:,1], nptv[:,3] .* 10^16, color = :blue)
        lines!(axr, nptv[:,1], nptv[:,4], color = :red)
    end

    llabel = Label(fig[1:2,0], text = "ciśnienie [Pa*10^16]", color = :blue, rotation = pi/2)
    rlabel = Label(fig[1:2,3], text = "temperatura [K]", color = :red, rotation = pi/2)
    xlabel = Label(fig[3,1:2], text = "x [m]", color = :black, font = :bold)
    
    return fig
end

# x | dens | p | temp | v | nr 


function hist_begin_end_plot(dir::String)
    file_end  = "hist_end.dat"
    file_begin  = "hist_begin.dat"
    h_end = readdlm(joinpath(dir, file_end))
    h_begin = readdlm(joinpath(dir, file_begin))
    fig = Figure()
    
    ax = Axis(fig[1,1],ylabel = "rozkłąd cząstek", xlabel = "prędkość")
    lines!(ax,h_end[:,1],h_end[:,2], label = "rozkład końcowy",  color = :black, linewidth = 2)
    lines!(ax,h_begin[:,1],h_begin[:,2], label = "rozkład początkowy",  color = :blue, linewidth = 2)
    axislegend()

    return fig
end
