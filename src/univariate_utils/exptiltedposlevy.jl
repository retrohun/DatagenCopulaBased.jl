#TODO fill
"""
    ExpTiltedPosLevy(α)

"""
struct ExpTiltedPosLevy{T<:Number} <: ContinuousUnivariateDistribution
    α::T

    function ExpTiltedPosLevy{T}(α::Real) where T
        1 <= α || throw(AssertionError("Parameter α needs to be at least 1"))
        new{T}(α)
    end
end

ExpTiltedPosLevy(α::T) where T<:Number = ExpTiltedPosLevy{T}(α)

#TODO: correct
"""
  rand(d::ExpTiltedPosLevy, v::Vector{Float64})


Returns a Vector{Floats} genrated from the expotencialy tilted levy stable pdf
f(x; V0, α) = exp(-V0^α) g(x; α)/exp(-V0), where g(x; α) is a stable Levy pdf
with parameters α = 1/θ, β = 1, γ = (cos(pi/(2*θ)))^θ and δ = 0.


"""
function rand(d::ExpTiltedPosLevy, v::Vector{Float64})
   t = length(v)
   ret = zeros(t)
   for i in 1:t
       x = rand(StableGumbel(d.α))
       u = rand()
       while exp(-v[i]^d.α*x)/(35*exp(-v[i])) < u
           x = rand(StableGumbel(d.α))
           u = rand()
       end
       ret[i] = x
   end
   ret.*v.^d.α
end
