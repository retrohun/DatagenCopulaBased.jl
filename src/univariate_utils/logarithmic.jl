import Base: mean, var, minimum, maximum
import Distributions: minimum, maximum, sampler, pdf, cdf, mean, mode, var, mgf, cf, insupport, quantile

#TODO: fill
"""
    Logarithmic(p)

"""
struct Logarithmic{T<:Real} <: DiscreteUnivariateDistribution
    p::T

    function Logarithmic{T}(p::Real) where T
        zero(T)<p<one(T) || throw(AssertionError("Parameter p needs to be between 0 and 1"))
        new{T}(p)
    end
end

Logarithmic(p::T) where T<:Real = Logarithmic{T}(p)
Logarithmic(p::Integer) = Logarithmic(Float64(p))

minimum(::Logarithmic) = 1
maximum(::Logarithmic) = Inf

pdf(dist::Logarithmic, k::Int) = -dist.p^k/log(1-dist.p)/k
cdf(dist::Logarithmic, k::Int) = 1 + quadgk(t->t^k/(1-t), 0, dist.p)[1]/log(1-dist.p)

mean(dist::Logarithmic) = -dist.p/log(1-dist.p)/(1-dist.p)
mode(dist::Logarithmic) = 1
var(dist::Logarithmic) = -dist.p*(dist.p + log(1-dist.p))/((1-dist.p)*log(1-dist.p))^2

params(dist::Logarithmic) = (dist.p,)

function mgf(dist::Logarithmic, t::Real)
    t < -log(dist.p) || throw(ArgumentError("t needs to be smaller than -ln(p)"))
    log(1-dist.p*exp(t))/log(1-dist.p)
end

cf(dist::Logarithmic, t::Real) = log(1-dist.p*exp(t))/log(1-dist.p)

function pgf(dist::Logarithmic, z::Complex128)
    abs(z) < 1/dist.p || throw(ArgumentError("argument needs to satisfy |z|<1/p"))
    log(1-dist.p*z)/log(1-dist.p)
end

function quantile(dist::Logarithmic, p::T) where T<:Real
    zero(T) <= p <= one(T) || throw(ArgumentError("p needs to be a probability"))
    i = 1
    sum = pdf(dist, minimum(dist))
    while p > sum
        i += 1
        sum += pdf(dist, i)
    end
    i
end
