@inline function Vg(t::Float64)
    # static parameters 
    ν = 1. *10^(9) # GHz
    t₀= 7.5 *10^(-9) # ns
    σ = 0.75 *10^(-9) # ns
    
    return sin(2π*ν*t)*exp(-(t-t₀)^2/(2*σ^2))
end


function u_plot(x::Vector{Float64},u::Vector{Float64})
    fig = Figure()
    ax = Axis(fig[1,1],xlabel = "x", ylabel = "u")
    lines!(x,u)
    fig
end


function u_plot(x::Vector{Float64},u::Vector{Float64},u_exact::Vector{Float64}, title::String)
    fig = Figure()
    ax = Axis(fig[1,1],xlabel = "x", ylabel = "u", title = title)
    scatter!(x,u,label = "monte carlo", color = :red)
    lines!(x,u_exact, label = "dokładny")
    axislegend()
    fig
end

function solve_mc_for_x(x,time,params,nchains)
    sol = Vector{Float64}(undef,length(x))

    for i in eachindex(x)
        (b1,b2,sol[i]) = mc_solve(nchains,x[i],time,params)
    end
    return sol
end

function animate_in_time(params)
    nchains = 10^3
    x = Vector(LinRange{Float64}(0,2,200))
    
    time = Observable(15*10^(-9))

    sol = @lift solve_mc_for_x(x,$time,params,nchains)
    uexact = @lift exact(x,$time,params)

    fig = Figure()
    ax = Axis(fig[1,1],xlabel = "x", ylabel = "u", title = @lift("t = $(round($time*10^9, digits = 1))"))
    scatter!(x,sol,label = "monte carlo", color = :red)
    lines!(x,uexact, label = "dokładny")
    axislegend()

    framerate = 100
    timestamps = range(10., 30., step=1/framerate)

    record(fig, "time_animation_v2.mp4", timestamps;
        framerate = framerate) do t
        time[] = t*10^(-9)
    end
end

function diff_plot(x::Vector{Float64},diff::Vector{Float64}, title::String)
    fig = Figure()
    ax = Axis(fig[1,1],xlabel = "x", ylabel = "diff(u)", title = title)
    lines!(x,diff, label = "Różnica")
    axislegend()
    fig
end

function do_whole_plot(x::Vector{Float64},sol::Vector{Vector{Float64}},exact::Vector{Float64}, title_text::String)
    f = Figure(size = (2048,1024))
    
    ax1 = Axis(f[1,1],xlabel = "x", ylabel = "u", title = "npaths=1000")
    scatter!(f[1,1],x,sol[1],label = "monte carlo",color = :red)
    lines!(f[1,1],x,exact,label = "dokładny")
    axislegend(ax1)

    ax2 = Axis(f[1,2],xlabel = "x", ylabel = "u", title = "npaths=10000")
    scatter!(f[1,2],x,sol[2],label = "monte carlo",color = :red)
    lines!(f[1,2],x,exact,label = "dokładny")
    axislegend(ax2)

    ax3 = Axis(f[1,3],xlabel = "x", ylabel = "u", title = "npaths=100000")
    scatter!(f[1,3],x,sol[3],label = "monte carlo",color = :red)
    lines!(f[1,3],x,exact,label = "dokładny")
    axislegend(ax3)    

    max_diff =1.2*maximum(abs.(sol[1] .- exact))

    ax4 = Axis(f[2,1],xlabel = "x", ylabel = "abs(u - u_exact)")
    lines!(f[2,1],x,abs.(sol[1] .- exact), color = :black)
    ylims!(ax4,0,max_diff)

    ax5 = Axis(f[2,2],xlabel = "x", ylabel = "abs(u - u_exact)")
    lines!(f[2,2],x,abs.(sol[2] .- exact), color = :black)
    ylims!(ax5,0,max_diff)

    ax6 = Axis(f[2,3],xlabel = "x", ylabel = "abs(u - u_exact)")
    lines!(f[2,3],x,abs.(sol[3] .- exact), color = :black)   
    ylims!(ax6,0,max_diff)

    f[0, :] = Label(f, title_text, fontsize = 20)

    display(f)

    f
end
