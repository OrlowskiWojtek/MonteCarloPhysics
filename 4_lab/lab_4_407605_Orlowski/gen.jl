function circle(Rα,xα, yα)
    u1 = rand()
    u2 = rand()

    x = √(-2*log(u1))*sin(2π*u2)
    y  = √(-2*log(u1))*cos(2π*u2)

    r = √(x^2 + y^2)
    x = x/r
    y = y/r

    q = sqrt(rand())

    x = q*x*Rα + xα
    y = q*y*Rα + yα

    return (x,y)
end

function circle(circle_data::Tuple{Float64,Float64,Float64})
    return circle(circle_data[1],circle_data[2],circle_data[3])
end

function gen_n_in_circle(N::Int, circle_data::Tuple{Float64,Float64,Float64})
    
    x = Vector{Float64}(undef,N)
    y = Vector{Float64}(undef,N)

    for i in 1:N
        (x[i],y[i]) = circle(circle_data)
    end

    return (x,y)
end

function gen_circle_label(circle::Tuple{Float64,Float64,Float64})

    s = "R=$(circle[1]);x=$(circle[2]),y=$(circle[3])"

    return s
end

#here there is no need to check for circle_1, but making it in case of adjustments
function θ(x,y,circle_1::Tuple{Float64,Float64,Float64},circle_2::Tuple{Float64,Float64,Float64})
    in_circle_1::Bool = ((x-circle_1[2])^2 + (y-circle_1[3])^2) < circle_1[1]^2
    in_circle_2::Bool = ((x-circle_2[2])^2 + (y-circle_2[3])^2) < circle_2[1]^2

    if(in_circle_1 && in_circle_2)
        return 1.
    else
        return 0.
    end
end

#same thing given that x,y are from circle_1 at 100% 
function θ(x,y,circle_2::Tuple{Float64,Float64,Float64})
    in_circle_2::Bool = ((x-circle_2[2])^2 + (y-circle_2[3])^2) < circle_2[1]^2

    if(in_circle_2)
        return 1.
    else
        return 0.
    end
end

function calc_common_area(n, circle_1::Tuple{Float64,Float64,Float64},circle_2::Tuple{Float64,Float64,Float64})
    s = 0

    for i in 1:n
        (x,y) = circle(circle_1)
        s +=  θ(x,y,circle_2)
    end

    μ1 = π*circle_1[1]^2 * s/n
    μ2 = π*circle_1[1]^2*μ1 
    σ = √((μ2 - μ1^2)/n)

    return (μ1,μ2,σ)
end

function calc_common_area(circle_1::Tuple{Float64,Float64,Float64},circle_2::Tuple{Float64,Float64,Float64}, save_at)
    s = 0

    n = save_at[end]

    counter = 1
    s_vec = Vector{Float64}(undef,length(save_at))

    for i in 1:n
        (x,y) = circle(circle_1)
        s +=  θ(x,y,circle_2)
        if(i == save_at[counter])
            s_vec[counter] = s
            counter += 1
        end    
    end

    μ1 = π*circle_1[1]^2 .* s_vec ./ save_at
    μ2 = π*circle_1[1]^2 .* μ1 
    σ = sqrt.((μ2 .- μ1.^2) ./ save_at)

    return (μ1,σ)
end

function gen_circles_plot(circle_1, circle_2, n)
    (x,y) = gen_n_in_circle(n,circle_2)
    plt = scatter(x,y, label = "Kβ", aspect_ratio=:equal, legendfontsize = 23)
    
    (x,y) = gen_n_in_circle(n,circle_1)
    scatter!(x,y, label = "Kα")
    return plt
end

