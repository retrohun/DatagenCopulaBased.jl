#TODO: remove
function logseriesquantile(p::Float64, v::Vector{Float64})
  dist = Logarithmic(p)
  map(x->quantile(dist, x), v)
end

#TODO: remove
function levyel(θ::Union{Int, Float64})
  rand(StableGumbel(θ))
end # rand δ -> μ, γ -> c

#TODO: remove
function levygen(θ::Union{Int, Float64}, u::Vector{Float64})
  rand(StableGumbel(θ), u)
end #u - vector rand()



#TODO remove
function tiltedlevygen(V0::Vector{Float64}, α::Float64)
  rand(ExpTiltedPosLevy(α), V0)
  #=t = length(V0)
  ret = zeros(t)
  for i in 1:t
    x = levyel(α)
    u = rand()
    while exp(-V0[i]^α*x)/(35*exp(-V0[i])) < u
      x = levyel(α)
      u = rand()
    end
    ret[i] = x
  end
  ret.*V0.^α=#
end #rand ExpTiltedPosStable


"""
  Ginv(y::Float64, α::Float64)

Returns Float64, helper for the joe/frank nested copula generator
"""
Ginv(y::Float64, α::Float64) = ((1-y)*gamma(1-α))^(-1/α)

"""
  InvlaJ(n::Int, α::Float64)

Returns Float64, n-th element of the inverse laplacea transform of generator of Joe nested copula
"""
InvlaJ(n::Int, α::Float64) = 1-1/(n*beta(n, 1-α))

"""
  sampleInvlaJ(α::Float64, v::Float64)

Returns Int, a sample of inverce laplacea transform of generator of Joe nested copula,
given a parameter α and a random numver v ∈ [0,1]
"""

function sampleInvlaJ(α::Float64, v::Float64)
  if v <= α
    return 1
  else
    G = Ginv(v, α)
    if G > 2^62
      return 2^62
    else
      return (InvlaJ(floor(Int, G), α) < v)? ceil(Int, G): floor(Int, G)
    end
  end
end

"""
  elInvlaF(θ₁::Float64, θ₀::Float64)

Returns Int, a single sample of the inverse laplacea transform of the generator
of nested Frank copula
"""

function elInvlaF(θ₁::Float64, θ₀::Float64)
  c1 = 1-exp(-θ₁)
  α = θ₀/θ₁
  if θ₀ <= 1
    v = rand()
    X = logseriesquantile(c1, rand(1))[1]
    while v > 1/((X-α)*beta(X, 1-α))
      v = rand()
      X = logseriesquantile(c1, rand(1))[1]
    end
    return X
  else
    v = rand()
    X = sampleInvlaJ(α, rand())
    while v > c1^(X-1)
      X = sampleInvlaJ(α, rand())
      v = rand()
    end
    return X
  end
end

"""
  nestedfrankgen(θ₁::Float64, θ₀::Float64, V0::Vector{Int})

Return vector of int, samples of inverse laplacea trensform of nested
Frak copula given parametes and V0 - vector of samples if invlaplace of perents copula
"""

function nestedfrankgen(θ₁::Float64, θ₀::Float64, V0::Vector{Int})
  if nprocs() == 1
    return map(k -> sum([elInvlaF(θ₁, θ₀) for j in 1:k]), V0)
  end
  u = SharedArray{Float64}(length(V0))
  @sync @parallel for i = 1:length(V0)
    u[i] = sum([elInvlaF(θ₁, θ₀) for j in 1:V0[i]])
  end
  Array(u)
end
