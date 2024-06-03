function r_avg(atoms::Vector{Atom})
    sum = 0.
    for atom in atoms
        sum += atom.r
    end
    return sum/(length(atoms))
end

function pcf(atoms::Vector{Atom},M)
    ravg = r_avg(atoms)
    r_max = 2.5 * ravg
    Δr = r_max/M
    Ω = 4π * ravg^2 
    pcf = zeros(Float64,M)
    n = length(atoms)

    for i in eachindex(atoms)
        for j in (i+1):n
            r = distance(atoms[i],atoms[j])
            m = Int(floor(r/Δr)+1)
            if m <= M
                pcf[m] += 2*Ω/(n^2 * 2π * r * Δr)
            end
        end
    end

    radius = range(0,r_max,length = M)
    return (pcf,radius, Δr)
end


