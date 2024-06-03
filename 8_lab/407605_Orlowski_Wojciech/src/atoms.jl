mutable struct Atom
    x::Float64
    y::Float64
    z::Float64
    θ::Float64
    ϕ::Float64
    r::Float64
end

function r(a::Atom)
    return sqrt(a.x^2 + a.y^2 + a.z^2)
end

# From spherical to cartesian transform
function x(a::Atom)
    return a.r * sin(a.θ) * cos(a.ϕ)
end

function y(a::Atom)
    return a.r * sin(a.θ) * sin(a.ϕ)
end

function z(a::Atom)
    return a.r * cos(a.θ)
end

function calc_cartesian!(a::Atom)
    a.x = x(a)
    a.y = y(a)
    a.z = z(a)
    a
end

function calc_cartesian(a::Atom)
    return [x(a), y(a), z(a)]
end

function calc_spherical!(a::Atom)
    a.r = r(a)
    a.θ = acos(a.z/a.r)
    a.ϕ = acos(a.x/(a.r * sin(a.θ)))
    a
end

function distance(a::Atom, b::Atom)
    return sqrt((a.x - b.x)^2 + (a.y - b.y)^2 + (a.z - b.z)^2)
end

function r⃗(a::Atom, b::Atom)::Vector{Float64}
    return [b.x - a.x, b.y - a.y, b.z - a.z]
end

function read_atoms(atom_pos::Matrix)
    atoms = Vector{Atom}(undef,0)
    for a in eachrow(atom_pos)
        atom = Atom(a[1],a[2],a[3],0.0,0.0,0.0)
        atom.r = r(atom)
        push!(atoms,atom)
    end
    return atoms
end

function rand_atoms(r_start,n)
    atoms = Vector{Atom}(undef,n)
    
    for i in 1:n
        atoms[i] = Atom(0,0,0,rand()*π, rand()*2π, r_start)
        calc_cartesian!(atoms[i])
    end

    return atoms    
end

function calc_r_avg(atoms::Vector{Atom})
    sum = 0 
    for atom in atoms
        sum += atom.r
    end

    return sum / length(atoms)
end

function save_atoms_avogadro(atoms::Vector{Atom})
    file = open("save_data.xyz", "w")

    write(file, "60 \n \n")

    for atom in atoms
        write(file, "C $(atom.x) $(atom.y) $(atom.z)\n")
    end
    close(file)
end

function save_atoms(atoms::Vector{Atom}, name::String, energy, params,w)
    file = open("atoms_pos/"*name, "w")
    write(file, "# Otrzymana energia: $energy\n")
    write(file, "# $params\n")
    write(file, "# $w\n")

    for atom in atoms
        write(file, "$(atom.x) $(atom.y) $(atom.z)\n")
    end
    close(file)
end

Base.show(io::IO, a::Atom) = print(io,"Atom located at: ($(a.x), $(a.y), $(a.z)) ")

