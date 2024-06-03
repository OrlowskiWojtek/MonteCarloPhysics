# input all files
using Makie, DelimitedFiles

#for 3D visualization
using CairoMakie
CairoMakie.activate!(type = "svg")

#for high quality svg
using GLMakie
GLMakie.activate!()

## loading Libraries

include("atoms.jl")
include("brenner.jl")
include("pcf.jl")
include("plot.jl")
include("optimize.jl")

Threads.nthreads()

## test task

atoms = read_atoms(readdlm("data.dat",comments = true, comment_char = '#'))
@time println("Testowa energia wynosi: $(V_total(atoms,Brenner_params)) eV") 

GLMakie.activate!()
scene = plot_atoms(atoms,2.0,3.3)
save("plots/test_fulleren.png",scene, px_per_unit = 2)

CairoMakie.activate!(type = "svg")

(PCF, radius) = pcf(atoms,100)
plt=lines(radius,PCF,  axis = (; title = "PCF", xlabel = "r [Å]"))
save("plots/test_pcf.svg",plt)

## task 3

params = (n = 60, β_min = 1.0, β_max = 100., p = 2, itmax = 10^5, Wall = 10^(-4))
wr = 10^(-4)
wϕ = 0.05
wθ = 0.05
w = (wr = wr, wθ = wθ, wϕ = wϕ)

r_start = 3.5
# random atoms 
atoms = rand_atoms(r_start,params.n)

GLMakie.activate!()
scene = plot_atoms(atoms,2.0,3.3)
save("plots/task3_random_start.png",scene, px_per_unit = 2, size = (500, 500))

(v,β_v) = optimize(atoms,w, Brenner_params,params);

scene = plot_atoms(atoms,2.0,3.3)
save("plots/task3_after_optim.png",scene, px_per_unit = 2, size = (500, 500))

CairoMakie.activate!(type = "svg")

(PCF, radius, Δr) = pcf(atoms,100)
plt=lines(radius,PCF, axis = (; title = "PCF", xlabel = "r [Å]"))
save("plots/task3_pcf.svg",plt)

it_vec = [i for i in 100:100:params.itmax]
plt = plot_β_v(β_v,v,it_vec)
save("plots/task3_beta.svg",plt)

println("Średnia odległość od środka wyniosła $(calc_r_avg(atoms)) Å")

## task 4 - need to change ζ function

function ζ(atoms::Vector{Atom},i::Int,j::Int,p)
    sum = 0.
    for k in eachindex(atoms)
        if(k == i || k == j)
            continue
        end
        f_cut = fcut(distance(atoms[i],atoms[k]),p)
        if(f_cut == 0)
            continue
        end
        cosθ_val = cosθ(atoms,i,j,k)
        if(cosθ_val > 0)
            return 10.
        end
        sum += f_cut * g(cosθ_val,p)        
    end
    return sum
end

params = (n = 60, β_min = 1.0, β_max = 100., p = 2, itmax = 10^5, Wall = 10^(-4))
wr = 10^(-4)
wϕ = 0.05
wθ = 0.05
w = (wr = wr, wθ = wθ, wϕ = wϕ)

r_start = 3.5
# random atoms 
atoms = rand_atoms(r_start, params.n)

GLMakie.activate!()
scene = plot_atoms(atoms,2.0,3.3)
save("plots/task4_random_start.png",scene, px_per_unit = 2,size = (500, 500))

(v,β_vec) = optimize(atoms,w, Brenner_params,params);

println("Średnia odległość od środka wyniosła $(calc_r_avg(atoms)) Å")

#
scene = plot_atoms(atoms,2.0,3.3)
save("plots/task4_after_optim.png",scene, px_per_unit = 2,size = (500, 500))

CairoMakie.activate!()
it_vec = [i for i in 100:100:params.itmax]
plt = plot_β_v(β_vec,v,it_vec)
save("plots/task4_beta.svg",plt)

(PCF, radius, Δr) = pcf(atoms,100)
plt=lines(radius,PCF, axis = (; title = "PCF", xlabel = "r [Å]"))
save("plots/task4_pcf.svg",plt)


## task 5

r_start = 2.5
# random atoms 
atoms = rand_atoms(r_start, params.n)

GLMakie.activate!()
scene = plot_atoms(atoms,1.6,2.4)
save("plots/task5_random_start.png",scene, px_per_unit = 2, size = (500,500))

(v,β_vec,r_v) = optimize(atoms,w, Brenner_params,params);
println("Średnia odległość od środka wyniosła: $(r_v[end]) Å")

scene = plot_atoms(atoms,1.8,2.4)
save("plots/task5_after_optim.png",scene, px_per_unit = 2)

CairoMakie.activate!()
it_vec = [i for i in 100:100:params.itmax]
plt = plot_β_v(β_vec,v,it_vec)
save("plots/task5_beta.svg",plt)

plt = lines(it_vec,r_v; axis = (xlabel = "Iteracja", ylabel = "r_śr [Å]"))
save("plots/task5_r.svg",plt)

(PCF, radius, Δr) = pcf(atoms,100)
plt=lines(radius,PCF, axis = (; title = "PCF", xlabel = "r [Å]"))
save("plots/task5_pcf.svg",plt)

## task 6 

params = (n = 60, β_min = 1.0, β_max = 200. , p = 2.0, itmax = 10^4, Wall = 10^(-5))
wr = 10^(-5)
wϕ = 0.50
wθ = 0.50
w = (wr = wr, wθ = wθ, wϕ = wϕ)

r_start = 3.3

## random atoms 

atoms = rand_atoms(r_start, params.n)
(v,β_vec,r_v) = optimize(atoms,w, Brenner_params,params);    

GLMakie.activate!()
scene = plot_atoms(atoms,2.0,2.4)
save("plots/task6_try_3.png", scene)
save_atoms(atoms,"task_6_try_3",v[end],params,w)

CairoMakie.activate!()
it_vec = [i for i in 100:100:params.itmax]
plt = plot_β_v(β_vec,v,it_vec)
save("plots/task6_beta_3.svg",plt)

plt = lines(it_vec, r_v, axis = (xlabel = "iteracja", ylabel = "r [Å]"))
hlines!(r_v[end],label = "r = $(r_v[end]) Å", color = :red)
axislegend()
save("plots/task6_r_3.svg",plt)

## task 7

E_b = Vector{Float64}(undef,0)
r_b = Vector{Float64}(undef,0)

for n in 30:60
    params = (n = n, β_min = 1.0, β_max = 100. , p = 2.0, itmax = 10^5, Wall = 10^(-4))
    wr = 10^(-4)
    wϕ = 0.05
    wθ = 0.05
    w = (wr = wr, wθ = wθ, wϕ = wϕ)

    r_start = 2.5
    atoms = rand_atoms(r_start, params.n)
    
    (v,β_vec,r_v) = optimize(atoms,w, Brenner_params,params);    
    save_atoms(atoms,"task7/$n",v[end],params,w)

    scene = plot_atoms(atoms,2.0,2.4);
    save("plots/task7/$n.png",scene)

    push!(E_b,v[end]/n)
    push!(r_b,r_v[end])
end

CairoMakie.activate!()
plt = lines(30:60,E_b, axis = (xlabel = "Ilość atomów ", ylabel = "Energia wiązania na jeden atom [eV]"))
save("plots/eb.svg",plt)

plt = lines(30:60,r_b, axis = (xlabel = "Ilość atomów ", ylabel = "Odległość od środka układu współrzędnych [Å]"))
save("plots/rb.svg",plt)

## additional tasks

params = (n = 60, β_min = 1.0, β_max = 100., p = 2.0, itmax = 10^4, Wall = 10^(-4))
wr = 10^(-4)
wϕ = 0.05
wθ = 0.05
w = (wr = wr, wθ = wθ, wϕ = wϕ)

r_start = 2.5 
atoms = rand_atoms(r_start, params.n)

(v,β_vec,r_v, atom_pos) = optimize(atoms,w, Brenner_params,params);    

CairoMakie.activate!()
for i in [1,10,100,1000]
    (PCF, radius, Δr) = pcf(atom_pos[i],100)
    plt=lines(radius,PCF, axis = (; title = "PCF; it = $i", xlabel = "r [Å]"))
    display(plt)  
    save("plots/task_add_pcf_it$i.svg",plt)
end

GLMakie.activate!()
for i in [1,10,100,1000]
    scene = plot_atoms(atom_pos[i], 1.7,2.4)
    save("plots/task_add_it$i.png",scene)
end

plot_atoms(atoms, 2.0,3.5)
save_atoms_avogadro(atoms)

animate_atoms(atom_pos,"test2_anim")
