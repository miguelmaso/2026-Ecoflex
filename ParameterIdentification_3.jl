using JLD2
using CSV
using DataFrames
using Plots
using HyperFEM
using Optimization
using Optim, OptimizationOptimJL
using Metaheuristics, OptimizationMetaheuristics
include("src/Calibration.jl")


## MARK: Uniaxial Load-Unload

IL = 70.0 # mm
IW = 10.5 # mm
IT = 2.43 # mm
IA = IW*IT # mm^2

λ_max_list = [60,40,20]
v_list = [0.05,0.025,0.01]
Test_Num_list = [3,4,5]
Δt = 0.2
npoints_ = 100
OneCycleTest__list = Vector{ExperimentData}()
for λ_max in λ_max_list, v in v_list
    Raw_onecycle_tests = Vector{RawOneCycleTest}()
    for Test_Num in Test_Num_list
        push!(Raw_onecycle_tests,ReadData_OneCycle(λ_max,v,IL,IA,Test_Num))
    end
    push!(OneCycleTest__list,OneCycleTest(Raw_onecycle_tests, npoints_))
end
println(OneCycleTest__list)
length(OneCycleTest__list)
onecycle_tests = OneCycleTest__list

## MARK: Biaxial Relaxation

λ_max,v,IL,IA = 75, 0.05, 68, 68*2.9
RawEquiBiaxialRelaxTest_ = ReadData_EquiBiaxialRelax(λ_max,v,IL,IA)

t_max = 1000
npoints_ = 200
EquiBiaxialRelaxTest_ = EquiBiaxialRelaxTest(RawEquiBiaxialRelaxTest_, t_max, npoints_);

push!(onecycle_tests,EquiBiaxialRelaxTest_)

##  MARK: Uniaxial Load unload Plots

plots = []
for OneCycleTest_ in OneCycleTest__list
    push!(plots, plot(OneCycleTest_.λ,OneCycleTest_.σ./1000,marker=2,xlabel="λ",ylabel="P (kPa)",label = false, title = "λ_max = $(OneCycleTest_.λ_max) v = $(OneCycleTest_.v)"))
end
grid_plot = plot(plots..., layout=(3, 3), size=(900, 900))

p = plot(size=(1000,1000))
for OneCycleTest_ in OneCycleTest__list
    p = plot!(OneCycleTest_.λ,OneCycleTest_.σ./1000,marker=2,xlabel="λ",ylabel="P (kPa)",label = "λ_max = $(OneCycleTest_.λ_max) v = $(OneCycleTest_.v)",legend = :outerbottom)
end
display(p)
savefig("C:/Users/mjbarillas/Documents/LaTeX/Elastomer-X_DEA/Load_Unload_9Cycles_Data.pdf")

## MARK: uniaxial QuasiStatic
npoints_ = 100
window = 4
v = 0.0001
λ_max = 60
IL_ = 70
Raw_Quasi_list = Vector{RawQuasiTest}()
Raw_Quasi = ReadData_Quasi(λ_max,v,IL_,IA)
push!(Raw_Quasi_list,Raw_Quasi)
Quasi_test =  QuasiStaticTest(Raw_Quasi, npoints_, window)
length(Quasi_test.λ)
λ_max = 70
IL_ = 65
Raw_Quasi = ReadData_Quasi(λ_max,v,IL_,IA)
push!(Raw_Quasi_list,Raw_Quasi)
Quasi_test =  QuasiStaticTest(Raw_Quasi, npoints_, window)
Quasi_test =  QuasiStaticTest(Raw_Quasi_list, npoints_, window)
qs_tests = Vector{QuasiStaticTest}()
push!(qs_tests,Quasi_test)

plot(qs_tests[1].λ,qs_tests[1].σ./1000,marker=2,xlabel="λ",ylabel="P (kPa)",label=false)
# savefig("C:/Users/mjbarillas/Documents/LaTeX/Elastomer-X_DEA/Quasi-Static_Data.pdf")


## MARK: Calibration Long term


build_longterm(C1, C2, C3) = Yeoh3D(λ=0.0, C10=C1, C20=C2, C30=C3)
pn = ["C10",  "C20",  "C30"]  # Parameter names
p0 = [  3e4,   -2e2,    3e0]  # Initial seed
lb = [1.0e3, -6.0e3,  0.0e0]  # Minimum search limits
ub = [3.0e5,  4.0e3,  1.0e3]  # Maximum search limits

build_longterm(μ) = IncompressibleNeoHookean2D(λ=0.0,μ=μ)
pn = ["μ"]  # Parameter names
p0 = [  3e4]  # Initial seed
lb = [1.0e3]  # Minimum search limits
ub = [1.0e6]  # Maximum search limits

build_longterm(μ_1,μ_2) = MooneyRivlin3D(λ=0.0, μ1=μ_1, μ2=μ_2)
pn = ["μ_1","μ_2"]  # Parameter names
p0 = [  4.0e5,   4.0e5]  # Initial seed
lb = [1.0e3, -1.0e3]  # Minimum search limits
ub = [1.0e7,  1.0e7]  # Maximum search limits


opt_func = OptimizationFunction((p,d) -> loss(build_longterm, p, d))
opt_prob = OptimizationProblem(opt_func, p0, qs_tests, lb=lb, ub=ub)
opt_long = solve(opt_prob, ParticleSwarm(lower=lb, upper=ub, n_particles=1000), maxiters=1000, maxtime=180)
opt_prob = OptimizationProblem(opt_func, opt_long.u, qs_tests)
opt_long = solve(opt_prob, Optim.NelderMead(), maxiters=100, maxtime=30)
sol_long = opt_long.u

long_term = build_longterm(sol_long...)
r2 = stats(build_longterm, sol_long, qs_tests, pn)

## MARK: Plot Calibration Long term
p = plot(xlabel="Stretch [-]", ylabel="P [KPa]")
plot_experiment!(long_term, getfirst(r -> true, qs_tests))
annotate_r2!(r2, 0.7)
savefig("C:/Users/mjbarillas/Documents/LaTeX/Elastomer-X_DEA/NeoHookean_Longterm.pdf")



## MARK: Calibration Visco

build_branch(μ, t) = ViscousIncompressible(IncompressibleNeoHookean3D(λ=0.0, μ=μ), τ=exp10(t))
# build_branch(C1, C2, C3, t) = ViscousIncompressible(Yeoh3D(λ=0.0, C10=C1, C20=C2, C30=C3), τ=exp10(t))
build_branches(p...) = map(splat(build_branch), Iterators.partition(p,2))
build_visco(p...) = GeneralizedMaxwell(build_longterm(sol_long...), build_branches(p...)...)
n_branches = 2
pn = reduce(vcat, ["\\(\\mu_$i\\)", "\\(\\text{log}(\\tau_$i)\\)"] for i in 1:n_branches)  # Parameter names
p0 = reduce(vcat, [  1e4,      1.0] for _ in 1:n_branches)  # Initial seed
# p0 = [
#     201388.9322396927,
#     1.9443069865791058,
#     4.418714577779876e6,
#     -0.8227763529584671,
#     440261.76718093106,
#     0.4633834793083379
# ]
lb = reduce(vcat, [  1e3,     -4.0] for _ in 1:n_branches)  # Lower search limits
ub = reduce(vcat, [  1e9,      4.0] for _ in 1:n_branches)  # Upper search limits
# p0 = reduce(vcat, [    3e4,   -2e2,    3e0,      1.0] for _ in 1:n_branches)  # Initial seed
# lb = reduce(vcat, [  1.0e3, -6.0e3,  0.0e0,     -1.0] for _ in 1:n_branches)  # Lower search limits
# ub = reduce(vcat, [  3.0e5,  4.0e3,  1.0e3,      4.0] for _ in 1:n_branches)  # Upper search limits

opt_func = OptimizationFunction((p,d) -> loss(build_visco, p, d))
opt_prob = OptimizationProblem(opt_func, p0, onecycle_tests, lb=lb, ub=ub)
opt_visco_ = solve(opt_prob, ParticleSwarm(lower=lb, upper=ub, n_particles=1000), maxiters=5000, maxtime=60*60*3)
opt_prob = OptimizationProblem(opt_func, opt_visco_.u, onecycle_tests)
opt_visco = solve(opt_prob, Optim.NelderMead(parameters = Optim.FixedParameters()), maxiters=500, maxtime=180)
sol_visco = opt_visco.u


long_term = build_longterm(sol_long...)
model = build_visco(sol_visco...)
r2 = stats(build_visco, sol_visco, onecycle_tests, pn)

##
p = plot(xlabel="λ [-]", ylabel="P [KPa]")
plot_experiments(model, filter(r -> r.λ_max ≈ 1.6, onecycle_tests), stretch_label, vel_label, "λ [-]", "P [KPa]")
annotate_r2!(r2, 0.5)
savefig("C:/Users/mjbarillas/Documents/LaTeX/Elastomer-X_DEA/$(length(onecycle_tests))Tests_$(n_branches)_NoNM_Branches_LambdaMax_1_6.pdf")

##
plot(onecycle_tests[10].λ_1_1,onecycle_tests[10].σ_1_1./1000.0,label= "Experiment",xlabel = "λ",ylabel = "P (kPa)")
plot!(onecycle_tests[10].λ_1_1,evaluate_stress(model, onecycle_tests[10].Δt_1, onecycle_tests[10].λ_1_1, onecycle_tests[10].λ_2_1)[1]./1000.0,label= "Model",xlabel = "λ",ylabel = "P (kPa)")
savefig("C:/Users/mjbarillas/Documents/LaTeX/Elastomer-X_DEA/$(length(onecycle_tests))Tests_$(n_branches)Branches_BiaxialLoad.pdf")

##

plot(onecycle_tests[10].λ_1_1,onecycle_tests[10].σ_1_1./1000.0,label= "Experiment",xlabel = "λ",ylabel = "P (kPa)",title = "Load Biaxial Relaxation")
plot!(onecycle_tests[10].λ_1_1,evaluate_stress(model, onecycle_tests[10].Δt_1, onecycle_tests[10].λ_1_1, onecycle_tests[10].λ_2_1)[1]./1000.0,label= "Model",xlabel = "λ",ylabel = "P (kPa)")
savefig("C:/Users/mjbarillas/Documents/LaTeX/Elastomer-X_DEA/$(length(onecycle_tests))Tests_$(n_branches)Branches_BiaxialLoad.pdf")

plot([(i-1)*onecycle_tests[10].Δt_2 for i in  1:lastindex(onecycle_tests[10].λ_1_2)],onecycle_tests[10].σ_1_2./1000.0,label= "Experiment",xlabel = "Time (s)",ylabel = "P (kPa)",title = "Relaxation")
plot!([(i-1)*onecycle_tests[10].Δt_2 for i in  1:lastindex(onecycle_tests[10].λ_1_2)],evaluate_stress(model,onecycle_tests[10].Δt_1,onecycle_tests[10].v,onecycle_tests[10].npoints_load,onecycle_tests[10].Δt_2,onecycle_tests[10].λ_max,onecycle_tests[10].t_max)[3]./1000.0,label= "Model",xlabel = "λ",ylabel = "P (kPa)")
savefig("C:/Users/mjbarillas/Documents/LaTeX/Elastomer-X_DEA/$(length(onecycle_tests))Tests_$(n_branches)Branches_BiaxialRelaxation.pdf")
