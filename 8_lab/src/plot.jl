using CairoMakie

function plot_atoms(atoms::Matrix,r, R)
    plt = scatter(atoms[:,1],atoms[:,2],atoms[:,3], markersize = 20)
    mesh!(Sphere(Point3f(0), R * 1f0), color = :white)
    for i in eachindex(atoms[:,1])
        for j in i+1:length(atoms[:,1])
            xi = atoms[i,1]
            yi = atoms[i,2]
            zi = atoms[i,3]
            xj = atoms[j,1]
            yj = atoms[j,2]
            zj = atoms[j,3]
    
            r0 = sqrt((xi-xj)^2 + (yi - yj)^2 + (zi - zj)^2)
            if(r0 < r)
                lines!([xi, xj], [yi, yj], [zi, zj], label = false, color = :black, lw = 2)
            end
        end
    end
    
    #display(plt)
    return plt
end


function plot_atoms(atoms_vec::Vector{Atom},r, R)
    atoms = Matrix{Float64}(undef, length(atoms_vec), 3)
    for i in eachindex(atoms_vec)
        atoms[i,1] = atoms_vec[i].x
        atoms[i,2] = atoms_vec[i].y
        atoms[i,3] = atoms_vec[i].z
    end
    return plot_atoms(atoms,r,R)
end

function plot_β_v(β,v,it_vec)    
    f = Figure()

    ax1 = Axis(f[1, 1], yticklabelcolor = :blue, xlabel = "Iteracja", ylabel = "Energia układu [eV]")
    ax2 = Axis(f[1, 1], yticklabelcolor = :red, yaxisposition = :right, ylabel = "Parametr β")
    hidespines!(ax2)
    hidexdecorations!(ax2)

    lines!(ax1, it_vec, v, color = :blue)
    lines!(ax2, it_vec, β, color = :red)

    return f
end

function generate_tex_input()
    com::String = "\\begin{figure}[H] \n"
    counter = 1

    for n in 30:2:60
        subfig =  "  \\begin{subfigure}{0.19\\textwidth}\n"
        subfig *= "    \\centering \n"
        subfig *= "    \\includegraphics[width=\\textwidth]{../plots/task7/cropped_$n.png}\n"
        subfig *= "    \\caption{n = $n}\n"
        subfig *= "  \\end{subfigure}\n"
        if(counter % 5 == 0)
            subfig *= "\\\\ \n"
        end
        com *= subfig
        counter += 1
    end

    com *= "\\end{figure}"
end

println(generate_tex_input())

function animate_atoms(atom_pos::Vector{Vector{Atom}}, name)    
    atoms_vec = atom_pos[begin]
    atoms = Matrix{Float64}(undef, length(atoms_vec), 3)
    
    for i in eachindex(atoms_vec)
        atoms[i,1] = atoms_vec[i].x
        atoms[i,2] = atoms_vec[i].y
        atoms[i,3] = atoms_vec[i].z
    end

    plot_atoms_x = Observable(atoms[:,1])
    plot_atoms_y = Observable(atoms[:,2])
    plot_atoms_z = Observable(atoms[:,3])
    
    frames = 1:length(atom_pos)
    fig, ax, plt = scatter(plot_atoms_x,plot_atoms_y,plot_atoms_z, markersize = 20)
    mesh!(Sphere(Point3f(0), 2.4f0), color = :white)

    record(fig, name*".mp4", frames;
        framerate = 30) do frame
            atoms_vec = atom_pos[frame]
            atoms = Matrix{Float64}(undef, length(atoms_vec), 3)
            
            for i in eachindex(atoms_vec)
                atoms[i,1] = atoms_vec[i].x
                atoms[i,2] = atoms_vec[i].y
                atoms[i,3] = atoms_vec[i].z
            end

            plot_atoms_x[] = atoms[:,1]
            plot_atoms_y[] = atoms[:,2]
            plot_atoms_z[] = atoms[:,3]    
    end
end

