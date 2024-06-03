# standardowe sprawozdania
# ocena tylko z raportów 
# raport w postaci pdf
# czas na raport to tydzień
# za każdy kolejny rozpoczęty tydzień 20pkt ucinane jest

function expected_value(x::Vector{Float64})
    E = sum(x)/length(x)
    return E
end

function variance(x::Vector{Float64})
    Z_1 = expected_value(x.^2)
    Z_2 = expected_value(x)^2

    return (Z_1 - Z_2)
end

function bernoulli(p)
    u = rand()
    if u <= p
        return 1
    else
        return 0
    end
end

function make_exp_val_plot(p,exp_val)
    exp_val_plot = plot(xlabel = "n", ylabel = "błąd względny <Z>")

    k = [2,3,4,5,6,7]
    n = 10 .^ k

    for i ∈ eachindex(p)
        plot!(n,abs.(exp_val[i,10 .^k] .- p[i]) ./ p[i],
            xaxis = :log, 
            yaxis = :log,
            label = "p = $(p[i])",
            lw = 2,
            minorgrid = true
            )
    end

    display(exp_val_plot)
    savefig(exp_val_plot, "galeria/rel_err_exp_val.png")
end

function make_var_plot(p,var)
    var_plot = plot(xlabel = "n", ylabel = "błąd względny wariancji")

    k = [2,3,4,5,6,7]
    n = 10 .^ k

    for i ∈ eachindex(p)
        plot!(n ,abs.(var[i,10 .^k ] .- (p[i] - p[i]^2)) ./ (p[i] - p[i]^2),
            xaxis = :log, 
            yaxis = :log,
            label = "p = $(p[i])",
            lw = 2,
            minorgrid = true
            )
    end

    display(var_plot)
    savefig(var_plot, "galeria/rel_err_var.png")
end
