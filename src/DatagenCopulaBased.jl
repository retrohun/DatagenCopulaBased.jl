module DatagenCopulaBased
  using Distributions
  using NLsolve
  using Combinatorics
  using HypothesisTests
  using Cubature
  using PyCall
  @pyimport scipy.cluster.hierarchy as sch


  # new
  include("univariate_utils/logarithmic.jl")
  include("univariate_utils/stable.jl")
  include("univariate_utils/exptiltedposlevy.jl")

  export rand
  # old
  include("sampleunivdists.jl")
  include("archcopcorrelations.jl")
  include("archcopulagendat.jl")
  include("nestedarchcopulagendat.jl")
  include("chaincopulagendat.jl")
  include("subcopulasgendat.jl")
  include("copulagendat.jl")
  include("marshalolkincopcor.jl")

  export archcopulagen, chaincopulagen, nestedarchcopulagen
  export cormatgen, convertmarg!
  export tstudentcopulagen, gausscopulagen
  export frechetcopulagen, marshalolkincopulagen, chainfrechetcopulagen
  export copulamix
end
