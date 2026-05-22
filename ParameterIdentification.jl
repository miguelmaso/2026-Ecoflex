using JLD2
using Plots
using HyperFEM
using Optimization
using Optim, OptimizationOptimJL
using Metaheuristics, OptimizationMetaheuristics
# using ParallelParticleSwarms, KernelAbstractions

include("src/Calibration.jl")
# using .Calibration


## Load data

data = begin
  d = load(abspath(dirname(@__FILE__),"data/20260521/GoodSample2&1_Cyclic_Relax_Quasi_TestData.jld2"))
  d["Test_Data"]
end

qs_tests = Vector{QuasiStaticTest}()
creep_tests = Vector{CreepTest}()
onecycle_tests = Vector{OneCycleTest}()

for (k,v) in data
  if occursin(r"^Test = [34].*", k)
    test = OneCycleTest(v)
    push!(onecycle_tests, test)
  elseif occursin(r"^Test = Relax*", k)
    test = CreepTest(v)
    push!(creep_tests, test)
  elseif occursin(r"^Test = Quasi*", k)
    test = QuasiStaticTest(v)
    push!(qs_tests, test)
  end
end


## Data inspection

println(onecycle_tests)
println(creep_tests)
println(qs_tests)

test_1 = getfirst(t -> t.v≈0.025 && t.λ_max≈1.6, onecycle_tests)
display(plot(test_1.λ, test_1.σ, label=nothing))
test_2 = qs_tests[1]
display(plot(test_2.λ, test_2.σ, label=false))


## Step 1: Hyperelastic characterization

build_long(C1, C2, C3) = Yeoh3D(λ=0.0, C10=C1, C20=C2, C30=C3)
pn = ["C10",  "C20",  "C30"]  # Parameter names
p0 = [  3e4,   -2e2,    3e0]  # Initial seed
lb = [1.0e3, -2.0e3,  0.0e0]  # Minimum search limits
ub = [2.0e5,  2.0e3,  2.0e2]  # Maximum search limits

opt_func = OptimizationFunction((p,d) -> loss(build_long, p, d))
opt_prob = OptimizationProblem(opt_func, p0, qs_tests, lb=lb, ub=ub)
opt_long = solve(opt_prob, ParticleSwarm(lower=lb, upper=ub, n_particles=1000), maxiters=1000, maxtime=60)
opt_prob = OptimizationProblem(opt_func, opt_long.u, qs_tests)
opt_long = solve(opt_prob, Optim.NelderMead(), maxiters=100, maxtime=30)
sol_long = opt_long.u


## Visualization of the long-term component

long_term = build_long(sol_long...)
r2 = stats(build_long, sol_long, qs_tests, pn)
p = plot(xlabel="Stretch [-]", ylabel="Stress [KPa]")
plot_experiment!(long_term, getfirst(r -> true, qs_tests))
annotate_r2!(r2, 0.7)
display(p);


## Step 2: Viscoelastic characterization

build_branch(μ, t) = ViscousIncompressible(IncompressibleNeoHookean3D(λ=0.0, μ=μ), τ=exp10(t))
build_branches(p...) = map(splat(build_branch), Iterators.partition(p,2))
build_visco(p...) = GeneralizedMaxwell(build_longterm(sol_long...), build_branches(p...)...)
n_branches = 2
pn = reduce(vcat, ["μ$i", "t$i"] for i in 1:n_branches)  # Parameter names
p0 = reduce(vcat, [  1e4,   1.0] for _ in 1:n_branches)  # Initial seed
lb = reduce(vcat, [  1e3,  -1.0] for _ in 1:n_branches)  # Lower search limits
ub = reduce(vcat, [  1e5,   4.0] for _ in 1:n_branches)  # Upper search limits
