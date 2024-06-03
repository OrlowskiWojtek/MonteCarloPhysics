function init_V(Vl,Vb,Vt,nx,ny,Δ)
    A = zeros(nx+1,ny+1)
    V = OffsetArrays.Origin(0)(A)
    
    for i in 0:nx
        V[i,0] = Vb*sin(π*Δ*i/(Δ*nx))
        V[i,ny] = Vt*sin(π*Δ*i/(Δ*nx)) 
    end

    for j in 0:ny
        V[0,j] = Vl * sin(π*Δ*j/(Δ*ny))
    end
    
    return V
end

function relax(V,ω,Δ,ε,Aρ)
    for i in 1:(size(V,1)-2), j in 1:(size(V,2)-2)
        V[i,j] = (1-ω)*V[i,j] + ω/4 * (V[i+1,j] + V[i-1,j] + V[i,j+1] + V[i,j-1] + Δ^2/ε * Aρ[i,j]) 
    end
end

function energy_F(V,Aρ)
    sum = 0
    for i in 1:(size(V,1)-2), j in 1:(size(V,2)-2)
        Ex = (V[i+1,j] - V[i-1,j])/(2*Δ)
        Ey = (V[i,j+1] - V[i,j-1])/(2*Δ)
        sum += (1/2 * (Ex^2 + Ey^2) - Aρ[i,j]*V[i,j])
    end
    return sum
end

function solve(V,Aρ,p)
    itmax = 10_000
    tol = 10^(-6)
    ω = 1.8
    ny = size(V,2)-1
    nx = size(V,1)-1
    F0 = energy_F(V,Aρ)
    F1 = F0
    it = 0
    while(true)
        F1 = F0
        relax(V,ω,p.Δ,p.ε,Aρ)
        for j in 1:ny-1
            V[nx,j] = V[nx-1,j] 
        end
        F0 = energy_F(V,Aρ)
        it+=1
        if(it > itmax || abs((F0-F1)/F0) < tol)
            break
        end
    end
    
    println("Ended on iteration $it with energy change $(abs((F0-F1)/F0))")
end
