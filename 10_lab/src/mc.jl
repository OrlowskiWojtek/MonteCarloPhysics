# monte carlo algorithm for solving integral equation


function mc_solve(npaths::Int,x_start::Float64,t_start::Float64,p)
    μ = p.G / p.C
    c = 1. / sqrt(p.L * p.C)
    
    R0 = sqrt(p.L/p.C)
    
    λ = 0.5 * (p.R/p.L - p.G/p.C)
    
    ζ = R0 / (R0 + p.Rg)
    Γg= (p.Rg - R0)/(p.Rg + R0)
    Γl = (p.Rl - R0)/(p.Rl + R0)
    
    fxt = 0.
    bxt = 0.
    
    for i in 1:2
        suma = 0.
        for n in 1:npaths
            x = copy(x_start)
            t = copy(t_start)
            η = 1.
            sign = (-1)^i
            while(t>0)
                s = -log(rand()) / (λ + μ)
                if(sign == -1)
                    if((x - c*s) > 0)
                        η *= λ/(λ + μ)
                    else
                        s = x/c
                        suma += η*ζ*Vg(t - s)
                        η *= Γg
                    end
                    x -= c*s
                    t -= s
                elseif(sign == 1)
                    if((x+c*s) < p.l)
                        η *= λ/(λ+μ)
                    else
                        s = (p.l-x)/c
                        η *= Γl
                    end
                    x += c*s
                    t -= s
                end
                sign = -1*sign
            end
        end

        if(i == 1)
            fxt = suma/npaths
        elseif (i == 2)
            bxt = suma/npaths
        end
    end
    return (fxt,bxt,fxt+bxt)
end



