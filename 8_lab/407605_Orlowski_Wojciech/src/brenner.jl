Brenner_params = (R0 = 1.315, R1 = 1.7, R2 = 2.00,
                De = 6.325, S = 1.29, λ = 1.5, 
                δ = 0.80469, a0 =  0.011304, c0 = 19.,
                d0 = 2.5)

function fcut(r, p)
    if(r > p.R2)
        return 0
    end
    if(r <= p.R1)
        return 1.
    end
    if(r > p.R1 && r <= p.R2)
        return 0.5 * (1. + cos((r - p.R1)/(p.R2 - p.R1) * π))
    end
end

function Vr(r,p)
    return p.De/(p.S - 1.) * exp(-sqrt(2*p.S)*p.λ *(r - p.R0))
end

function Va(r,p)
    return p.De * p.S / (p.S - 1.) * exp(-sqrt(2/p.S)*p.λ *(r - p.R0))
end

function cosθ(atoms::Vector{Atom}, i::Int, j::Int, k::Int)
    #rij = r⃗(atoms[i],atoms[j])
    #rik = r⃗(atoms[i],atoms[k])

    scalar_ij_ik = (atoms[j].x - atoms[i].x) * (atoms[k].x-atoms[i].x) + (atoms[j].y - atoms[i].y) * (atoms[k].y-atoms[i].y) + (atoms[j].z - atoms[i].z) * (atoms[k].z-atoms[i].z)
    rij_dist = sqrt((atoms[j].x - atoms[i].x)^2  + (atoms[j].y - atoms[i].y)^2 + (atoms[j].z - atoms[i].z)^2)
    rik_dist = sqrt((atoms[k].x - atoms[i].x)^2  + (atoms[k].y - atoms[i].y)^2 + (atoms[k].z - atoms[i].z)^2)
    return scalar_ij_ik/(rij_dist*rik_dist)
    #return (rij'*rik)/(sqrt(rij'*rij)*sqrt(rik'*rik))
end

function g(cosθ,p)
    return p.a0 * (1 + p.c0^2/p.d0^2 - p.c0^2/(p.d0^2 + (1 + cosθ)^2))
end

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
        sum += f_cut * g(cosθ(atoms,i,j,k),p)        
    end
    return sum
end

function B(atoms::Vector{Atom},i::Int,j::Int,p)
    return (1. + ζ(atoms,i,j,p))^(-p.δ)
end

function B̄(atoms::Vector{Atom},i::Int, j::Int,p)
    return (B(atoms,i,j,p) + B(atoms,j,i,p))/2
end

function Vi(atoms::Vector{Atom},i,p)
    sums = 0
     
    for j in eachindex(atoms)
        if(j == i)
            continue
        end
        r_ij = distance(atoms[i],atoms[j])
        f_cut = fcut(r_ij,p)
        if(f_cut == 0)
            continue
        end 
        sums += f_cut*(Vr(r_ij,p) - B̄(atoms,i,j,p)*Va(r_ij,p))
    end
    return sums
end

function V_total(atoms::Vector{Atom},p)
    sum = 0
    for i in eachindex(atoms)
        sum+=Vi(atoms,i,p)
    end
    return 0.5 * sum
end

