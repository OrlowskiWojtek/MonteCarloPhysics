module CppExact
  using CxxWrap
  @wrapmodule(() -> joinpath("u_xt_exact/build/lib/","libuexactlib"))

  function __init__()
    @initcxx
  end
end

function exact(x::Vector{Float64},t::Float64,p)
  uexact = Vector{Float64}(undef, length(x))

  # double u_xt_exact(double x, double t, double t0,  double freq, double sigma, double R, double G, double L, double C, 
  # double Rg, double Rl, double length, int number_nodes, int n_sum_terms){
  R0 = sqrt(p.L/p.C)
  
  for i in eachindex(x)
      uexact[i] = CppExact.u_xt_exact(x[i],t,p.t₀,p.ν,p.σ,p.R,p.G,p.L,p.C,p.Rg,p.Rl,p.l,1000,100)
  end
  return uexact
end

