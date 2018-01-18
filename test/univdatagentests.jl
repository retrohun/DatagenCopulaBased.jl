
@testset "logseries dist" begin
  @test map(x->cdf(Logarithmic(0.01), x), 1:2) ≈ [0.994992, 0.999967] atol=1.0e-5

  dist = Logarithmic(0.9)
  @test map(x->quantile(dist, x), [0.25, 0.5, 0.75]) == [1, 2, 5]

  srand(43)
  d = Logarithmic(0.4)
  @test mean(d) ≈ 1.304 atol=1.0e-2
  @test std(d) ≈ 0.687 atol=1.0e-2

  v = rand(d, 1000000)
  @test skewness(v) ≈ 3.1 atol=1.0e-2
  @test kurtosis(v) ≈ 13.5 atol=1.0

end
@testset "stable levy dist" begin
  srand(43)
  @test levyel(2.) ≈ 0.18170339379413047
  srand(43)
  @test levygen(2., [0.2, 0.4, 0.6, 0.8]) ≈ [0.159748, 0.181703, 3.20539, 6.91497] atol=1.0e-4
  srand(43)
  @test tiltedlevygen([0.2, 0.4, 0.6], 2.) ≈ [0.271861, 0.0285145, 1.11889] atol=1.0e-5
end
@testset "nested copulas data generators" begin
  @test Ginv(0.5, 0.5) ≈ 1.2732395447351625
  @test InvlaJ(4, 0.5) ≈ 0.7265625
  @test sampleInvlaJ(0.5, 0.5) == 1
  @test sampleInvlaJ(0.5, 0.8) == 8
  srand(43)
  @test elInvlaF(4., 2.) == 7
  srand(43)
  @test elInvlaF(4., .5) == 16
  srand(43)
  @test nestedfrankgen(5., 3., [1, 1, 2]) == [3, 1, 49]
end
