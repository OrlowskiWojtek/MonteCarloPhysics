using Makie, CairoMakie, DelimitedFiles

## plotowanie i obróbka danych

function do_all_plots(name)
    trans_file_name = "data/trans"*name*".dat"
    abs_file_name = "data/abs_"*name*".dat"
    refl_file_name = "data/refl"*name*".dat"
    
    abs_data = readdlm(abs_file_name, comments = true)
    y = LinRange(0,0.06,length(abs_data[:,1]))
    x = LinRange(0,0.2, length(abs_data[1,:]))
    fig = Figure(size = (1024,768))
    ax, h = heatmap(fig[2:3,1],x,y,abs_data, colorrange = (0,70))
    ax.ylabel = "y [cm]"
    Colorbar(fig[2:3,2],h)

    trans_data = readdlm(trans_file_name)
    x = LinRange(0,0.2,length(trans_data))
    ax = Axis(fig[1,1], ylabel = "Transmisja")
    lines!(ax,x,trans_data[:,1], color = :black, linewidth = 3)

    refl_data = readdlm(refl_file_name)
    x = LinRange(0,0.2,length(refl_data))
    ax = Axis(fig[4,1],xlabel = "x [cm]", ylabel = "Odbicie")
    lines!(ax,x,refl_data[:,1], color = :red, linewidth = 3)

    save("plots/"*name*".pdf",fig)
    fig    
end

do_all_plots("task2_0")
do_all_plots("task2_1")
do_all_plots("task2_2")
do_all_plots("task2_3")
do_all_plots("task3_0")
do_all_plots("task3_1")
do_all_plots("task3_2")
do_all_plots("task3_3")

do_all_plots("test")
## checking if T + R + A / N == 1 -> all correct

function data_check(name)
    trans_file_name = "data/trans"*name*".dat"
    abs_file_name = "data/abs_"*name*".dat"
    refl_file_name = "data/refl"*name*".dat"

    N = 200_000

    A_sum = sum(readdlm(abs_file_name, comments = true))
    T_sum = sum(readdlm(trans_file_name))
    R_sum = sum(readdlm(refl_file_name))

    s = (A_sum + T_sum + R_sum)/N 
    println("Suma składowych wynosi: $s")
end

data_check("task2_0")
data_check("task2_1")
data_check("task2_2")
data_check("task2_3")
data_check("task3_0")
data_check("task3_1")
data_check("task3_2")
data_check("task3_3")

##

atan(0.6/0.8)*180/pi

