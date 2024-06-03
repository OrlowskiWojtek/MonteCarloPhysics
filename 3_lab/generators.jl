function box_muller()
    U₁ = rand()
    U₂ = rand()

    X = √(- 2 * log(1 - U₁)) * cos(2π*U₂)
    Y = √(- 2 * log(1 - U₁)) * sin(2π*U₂)

    return (X,Y)
end


function gen_n_box_muller(n::Int, radius::Bool = false)
    x = Vector{Float64}(undef,n)
    y = Vector{Float64}(undef,n)
    
    for i in 1:n
        (x[i],y[i]) = box_muller()
    end

    if(!radius)
        return (x,y)
    else
        return sqrt.(x .^2 .+ y .^2)
    end
end

function on_circle()
    (X,Y) = box_muller()

    Xc = X/√(X^2 + Y^2)
    Yc = Y/√(X^2 + Y^2)

    return (Xc,Yc)
end

function gen_n_on_circle(n::Int)
    x = Vector{Float64}(undef,n)
    y = Vector{Float64}(undef,n)
    
    for i in 1:n
        (x[i],y[i]) = on_circle()
    end

    return (x,y)
end

function in_circle()
    R = sqrt(rand())
    (Xc,Yc) = on_circle()
    (Xi,Yi) = (Xc*R, Yc*R)

    return (Xi, Yi)
end

function gen_n_in_circle(n::Int)
    x = Vector{Float64}(undef,n)
    y = Vector{Float64}(undef,n)
    
    for i in 1:n
        (x[i],y[i]) = in_circle()
    end

    return (x,y)
end

function gen_trans_matrix(α,b₁,b₂)
    R = [cos(α) -sin(α);
    sin(α) cos(α)]

    r₁ = R * b₁ * [1, 0]
    r₂ = R * b₂ * [0, 1]

    A = [r₁[1] r₂[1] ; r₁[2] r₂[2]]
    return A
end

function move_to_ellipse(α,b₁,b₂, X, Y)
    A = gen_trans_matrix(α,b₁,b₂)
    result = A*transpose([X Y])
    
    return (result[1,:], result[2,:])
end

function σₓ(X)
    s = 0
    s_squared = 0
    for i in X
        s += i
        s_squared += i^2
    end

    return s_squared/length(X) - (s/length(X))^2
end

function σ(X,Y)
    if(length(X) != length(Y))
        print("length of x and y do not match")
    end
    
    sum_xy = 0
    sum_x = 0
    sum_y = 0

    for i in eachindex(X)
        sum_xy += X[i]*Y[i]
        sum_x += X[i]
        sum_y += Y[i]
    end

    return sum_xy/length(X) - sum_x * sum_y/(length(X)^2)
end

function cov(X,Y)
    σ_XY = σ(X,Y)
    mat = [σₓ(X) σ_XY; σ_XY σₓ(Y)]
    return mat
end

function cov_coeff(cov_matrix::Matrix)
    r = cov_matrix[1,2]/sqrt(cov_matrix[1,1] * cov_matrix[2,2])
    return r
end

function cov_coeff(X::Vector, Y::Vector)
    cov_matrix = cov(X,Y)
    return cov_coeff(cov_matrix)
end

