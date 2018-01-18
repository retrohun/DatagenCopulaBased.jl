
import Base:rand

#TODO: fill
abstract type Stable{T<:Real} <: ContinuousUnivariateDistribution end

minimum(d::Stable{T}) where T = d.α < 1 && d.β == one(T) ? d.μ : -Inf
maximum(d::Stable{T}) where T = d.α < 1 && d.β == -one(T) ? d.μ : Inf

#TODO fill
"""
    StableGumbel
α = 1/θ
β = 1
c = cos^θ(π/(2θ))
μ = 0

"""
struct StableGumbel{T<:Real} <: Stable{T}
   θ::T

   function StableGumbel{T}(θ::Real) where T
      1 <= θ || throw(ArgumentError("θ needs to be at least 1"))
      new{T}(θ)
   end
end

function StableGumbel(θ::T) where T
   StableGumbel{T}(θ)
end

stability(d::StableGumbel{T}) where T = one(T)/d.θ
skewness(d::StableGumbel{T}) where T= one(T)
scale(d::StableGumbel{T}) where T = cos(π/2/d.θ)^d.θ
location(d::StableGumbel{T}) where T = zero(T)



function rand(d::StableGumbel)
   ϕ = π*rand()-π/2
   v = quantile.(Exponential(1.), rand())
   γ = scale(d)
   v = ((cos(pi/(2*d.θ)+(1/d.θ-1)*ϕ))/v)^(d.θ-1)
   γ*v*sin(1/d.θ*(pi/2+ϕ))*(cos(pi/(2*d.θ))*cos(ϕ))^(-d.θ)
end

#TODO: correct
"""
Return a Vectof(Floast} of  of pseudo cdf of Levy stable distribution with parameters
α = 1/θ, β = 1, γ = (cos(pi/(2*θ)))^θ and δ = 0, given a vector of Floats - u
"""
function rand(d::StableGumbel, v::Vector{Float64})
   p = invperm(sortperm(v))
   v = rand(d, length(v))
   sort(v)[p]
end
