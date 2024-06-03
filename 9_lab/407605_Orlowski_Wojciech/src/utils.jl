function ρ(ρmax::Float64,ρ_p::Float64,rmax::Vector{Float64},r::Vector{Float64})
    return ρmax * exp(-((rmax[1]/2 - r[1])^2 + (rmax[2]/2 - r[2])^2 )/(2*ρ_p^2))
end


function fill_rho(nx,ny,Δ,ρmax)
    xmax = Δ*nx
    ymax = Δ*ny
    ρ_p = xmax/10
    
    A = Array{Float64}(undef,nx+1,ny+1)

    ρ_array = OffsetArray(A,0:nx,0:ny)

    for i in 0:nx, j in 0:ny
        ρ_array[i,j] = ρ(ρmax,ρ_p,[xmax,ymax],[i*Δ,j*Δ])
    end

    return ρ_array
end

## plot scripts

function plot_V(V, interpolate::Bool = false)
    
    A = Matrix{Float64}(OffsetArrays.no_offset_view(V))
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Potencjał V", xticks = 0:5:size(A,1), yticks = 0:5:size(A,2))
    hm = heatmap!(ax, A, interpolate = interpolate, colormap = Reverse(:RdBu))
    Colorbar(fig[:, end+1], hm)
    fig
end


function do_all_plots(mc::MC_poisson, V_rel, name::String)
    fig = plot_V(mc.V,false)
    save("plots/"*name*"Vmc.svg",fig)

    A = Matrix{Float64}(OffsetArrays.no_offset_view(abs.(mc.V .- V_rel)))
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Różnica potencjałów", xticks = 0:5:size(A,1), yticks = 0:5:size(A,2))
    hm = heatmap!(ax, A, interpolate = false, colormap = :acton)
    Colorbar(fig[:, end+1], hm)

    save("plots/"*name*"Vdiff.svg",fig)

    A = Matrix{Float64}(OffsetArrays.no_offset_view(mc.σ))
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Odchylenie standardowe σ", xticks = 0:5:size(A,1), yticks = 0:5:size(A,2))
    hm = heatmap!(ax, A, interpolate = false)
    Colorbar(fig[:, end+1], hm)

    save("plots/"*name*"stddiff.svg",fig)

    A = Matrix{Float64}(OffsetArrays.no_offset_view(mc.S))
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Ułamek zaabsorbowanych łańcuchów", xticks = 0:5:size(A,1), yticks = 0:5:size(A,2))
    hm = heatmap!(ax, A, interpolate = false, colormap = :Oranges)
    Colorbar(fig[:, end+1], hm)
    
    save("plots/"*name*"chain_map.svg",fig)
end