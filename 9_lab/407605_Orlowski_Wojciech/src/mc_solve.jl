using ProgressMeter

struct MC_poisson
    nx::Int64
    ny::Int64
    N::Int64
    n::Int64
    Δ::Float64
    block::Bool
    V::OffsetArray{Float64}
    Aρ::OffsetArray{Float64}
    σ::OffsetArray{Float64}
    B::OffsetArray{Bool}
    S::OffsetArray{Float64}
end

function init_rho(mc::MC_poisson,ρmax::Float64)
    xmax = Δ*nx
    ymax = Δ*ny
    ρ_p = xmax/10.
    for i in 0:nx, j in 0:ny
        mc.Aρ[i,j] = ρ(ρmax,ρ_p,[xmax,ymax],[i*Δ,j*Δ])
    end
end

function init_boundary(mc::MC_poisson,Vl,Vb,Vt)
    xmax = mc.Δ * mc.nx
    ymax = mc.Δ * mc.ny
    for i in 0:mc.nx
        mc.V[i,0] = Vb*sin(π*mc.Δ*i/(xmax))
        mc.V[i,mc.ny] = Vt*sin(π*mc.Δ*i/(xmax)) 
        mc.B[i,0] = true
        mc.B[i,mc.ny] = true
    end

    for j in 0:mc.ny
        mc.V[0,j] = Vl * sin(π*mc.Δ*j/(ymax))
        mc.B[0,j] = true
    end
end

function init_mc(nx::Int,ny::Int,N::Int,n::Int,block::Bool,Δ::Float64,ρ_max::Float64, V_boundary)::MC_poisson    
    A = zeros(nx+1,ny+1)
    mc = MC_poisson(nx,ny,N,n,Δ,block,
        OffsetArrays.Origin(0)(copy(A)),
        OffsetArrays.Origin(0)(copy(A)), OffsetArrays.Origin(0)(copy(A)),
        OffsetArrays.Origin(0)(copy(A)), OffsetArrays.Origin(0)(copy(A)))
    
    init_rho(mc, ρ_max)
    init_boundary(mc,V_boundary...)
    return mc
end

function solve(mc::MC_poisson)
    p = Progress(mc.nx-1; desc = "Solving Poisson equation")
    @showprogress for i0 in 1:mc.nx-1
        for j0 in 1:mc.ny-1
            sum_V1 = 0.
            sum_V2 = 0.
            k_chains = 0
            for N in 1:mc.N
                i = i0
                j = j0
                g::Float64 = 0.
                for n in 1:mc.n
                    m = Int(floor(4*rand()))
                    if(m == 0) i -= 1 end
                    if(m == 1) i += 1 end
                    if(m == 2) j -= 1 end
                    if(m == 3) j += 1 end
                    # check for boundary reflection
                    if(i==(mc.nx+1)) i = mc.nx - 1 end
                    # check for boundary absorption
                    if(mc.B[i,j])
                        dV = mc.V[i,j] + g
                        sum_V1 += dV
                        sum_V2 += dV^2
                        k_chains += 1
                        break
                    end
                    g += mc.Aρ[i,j]*(mc.Δ)^2/4. # assumed that ε = 1 TODO change that if needed 
                end
            end
            V1 = sum_V1/k_chains
            V2 = sum_V2/k_chains
            mc.V[i0,j0] = V1
            mc.σ[i0,j0] = sqrt((V2-V1^2)/k_chains)
            if(mc.block) mc.B[i0,j0] = true end
            mc.S[i0,j0] = k_chains/mc.N
        end
    end
    for j in 0:ny
        mc.V[mc.nx,j] = mc.V[mc.nx-1,j] 
    end
end
