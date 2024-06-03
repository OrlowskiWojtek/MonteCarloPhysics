# algorytm wyrzażania
using ProgressMeter


function rand_move!(atoms::Vector{Atom},i,p,brenner_params, β)
    atom = atoms[i]
    
    Δr = atom.r*(2*rand()-1)*p.wr
    Δϕ = atom.ϕ*(2*rand()-1)*p.wϕ
    Δθ = atom.θ*(2*rand()-1)*p.wθ

    rn = atom.r + Δr
    ϕn = atom.ϕ + Δϕ
    θn = atom.θ + Δθ

    if(ϕn < 0) ϕn += 2π end
    if(ϕn > 2π) ϕn -= 2π end

    if(θn < 0) θn = atom.θ end
    if(θn > π) θn = atom.θ end    

    new_atom = Atom(0,0,0,θn,ϕn,rn)
    calc_cartesian!(new_atom)
    
    V_old = Vi(atoms,i,brenner_params)

    atoms[i] = new_atom
    V_new = Vi(atoms,i,brenner_params)


    p_acc = min(1,exp(-β*(V_new - V_old)))
    
    if(rand() <= p_acc)
        atoms[i] = new_atom
    else
        atoms[i] = atom
    end
end

function β(it,itmax,β_min,β_max,p)::Float64
    return (β_min + (it/itmax)^p * (β_max - β_min))
end

function random_r!(atoms::Vector{Atom}, Wall, brenner_params, β)
    E_old = V_total(atoms,brenner_params)
    
    U1 = rand()
    
    for i in eachindex(atoms)
        atoms[i].r = atoms[i].r*(1 + Wall * (2*U1 - 1))
        calc_cartesian!(atoms[i])
    end

    E_new = V_total(atoms,brenner_params)

    p_acc = min(1,exp(-β * (E_new - E_old)))

    if( rand() <= p_acc )
        return
    else
        for i in eachindex(atoms)
            atoms[i].r = atoms[i].r/(1 + Wall * (2*U1 - 1))
            calc_cartesian!(atoms[i])
        end 
    end
end

function optimize(atoms, w, brenner_params, params)
    β_vec = Vector{Float64}(undef,Int(params.itmax/1))
    V_vec = Vector{Float64}(undef,Int(params.itmax/1))
    r_vec = Vector{Float64}(undef,Int(params.itmax/1))
    atom_positions = Vector{Vector{Atom}}(undef,Int(params.itmax/1))
    counter::Int64 = 1

    p = Progress(params.itmax; desc="Optimizing structure")
    @showprogress for it in 1:params.itmax
        β_it = β(it,params.itmax,params.β_min,params.β_max,params.p)

        for i in eachindex(atoms)
            rand_move!(atoms, i, w, brenner_params,β_it)
        end

        random_r!(atoms,params.Wall,brenner_params,β_it)

        if(it % 1 == 0)
            V_vec[counter] = V_total(atoms,brenner_params)
            β_vec[counter] = β_it
            r_vec[counter] = calc_r_avg(atoms);
            atom_positions[counter] = copy(atoms)
            counter+=1
        end
    end

    println("Końcowa energia wynosi: $(V_total(atoms,Brenner_params)) eV") 
    println("Średnia odległość od środka wyniosła: $(r_vec[end]) Å")

    return (V_vec, β_vec, r_vec, atom_positions)
end
