
## particle locations 
function particle_location(dir::String)
    
    fig = Figure()
    xpos = [1,1,2,2,3,3]
    ypos = [1,2,1,2,1,2]
    its  = [601,701,801,901,1001,1101]

    skip = 20
    for i in eachindex(its)
        ax1 = Axis(fig[xpos[i],ypos[i]], title = "iteracja $(its[i])")
        file = "data_evolution_it$(its[i])"
        rv = readdlm(joinpath(dir,file), Float64, skipblanks = true)[begin:end-1,:]
        scatter!(ax1, rv[begin:skip:end,1],rv[begin:skip:end,2], markersize = 1.0, color = :black)
        xlims!(0,2)
        ylims!(0,0.5)
    end

    return fig
end

function plot_nptv(dir::String, its)
    fig = Figure(size = (1024,768))

    xpos = [1,2,3,4]
    ypos = [1,1,1,1]
    col  = [2, 3 ,4 ,5]
    col_name =  ["x", "gęstość", "ciśnienie", "temperatura", "prędkość"]

    for i in eachindex(col)
        title = col_name[col[i]]
        axl = Axis(fig[xpos[i], ypos[i]], title = title)
        if(col[i] < 5) hidexdecorations!(axl,grid = false) end

        for j in eachindex(its)
            nptv = readdlm(joinpath(dir,"data_nptv_it$(its[j])"))
            lines!(axl, nptv[:,1], nptv[:,col[i]], label = "$(its[j])", linewidth = 2)
        end
        fig[1:5, 2] = Legend(fig, axl, "iteracja", framevisible = false)
    end

    xlabel = Label(fig[5,1], text = "x", color = :black, font = :bold, tellwidth = false)

    return fig
end

