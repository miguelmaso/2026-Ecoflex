using JLD2
using CSV
using DataFrames
using Plots
using HyperFEM
using Optimization
using Optim, OptimizationOptimJL
using Metaheuristics, OptimizationMetaheuristics
# using ParallelParticleSwarms, KernelAbstractions

include("src/Calibration.jl")

## MARK: Initial Data
IL = 70.0 # mm
IW = 10.5 # mm
IT = 2.43 # mm
IA = IW*IT # mm^2

## MARK: Cyclic

λ_max,v,Test_Num = 60, 0.05, 3
Raw_OneCycle = ReadData_OneCycle(λ_max,v,IL,IA,Test_Num)
println(Raw_OneCycle)
plot(Raw_OneCycle.λ,Raw_OneCycle.P./1000,marker=2,xlabel="λ",ylabel="P (kPa)",label = "λ_max = $(Raw_OneCycle.λ_max) v = $(Raw_OneCycle.v)")

λ_max_list = [60,40,20]
v_list = [0.05,0.025,0.01]
Test_Num_list = [3,4,5,6]
Raw_onecycle_tests = Vector{RawOneCycleTest}()
for λ_max in λ_max_list, v in v_list, Test_Num in Test_Num_list
    push!(Raw_onecycle_tests,ReadData_OneCycle(λ_max,v,IL,IA,Test_Num))
end
println(Raw_onecycle_tests)

IL_list = [70,66]
plots = []
for λ_max in λ_max_list, v in v_list
    p = plot()
    for Test_Num in Test_Num_list, IL_ in IL_list
        try
                    Raw_OneCycle = ReadData_OneCycle(λ_max,v,IL_,IA,Test_Num)
                    p = plot!(Raw_OneCycle.λ,Raw_OneCycle.P./1000,marker=2,xlabel="λ",ylabel="P (kPa)",label = false, title = "λ_max = $(Raw_OneCycle.λ_max) v = $(Raw_OneCycle.v)")
        catch
        end
    end
    push!(plots,p)
end
grid_plot = plot(plots..., layout=(3, 3), size=(900, 900))

IL_ = 65
λ_max,v = 60, 0.05
Raw_Relax = ReadData_Relax(λ_max,v,IL_,IA)
p = plot()
i = 0
for λ_max in [60,40,20]
    Raw_Relax = ReadData_Relax(λ_max,v,IL_,IA)
    p = plot!(Raw_Relax.Time.+(i*100.0),Raw_Relax.P./1000,marker=2,xlabel="λ",ylabel="P (kPa)",label = "λ_max = $(Raw_Relax.λ_max) v = $(Raw_Relax.v)", legend=:topright)
    i+=1
end
display(p)

v = 0.0001
λ_max = 60
IL_ = 70
Raw_Quasi = ReadData_Quasi(λ_max,v,IL_,IA)
p = plot()
p = plot!(Raw_Quasi.λ,Raw_Quasi.P./1000,marker=2,xlabel="λ",ylabel="P (kPa)",label = "λ_max = $(Raw_Quasi.λ_max) v = $(Raw_Quasi.v)", legend=:topleft)
λ_max = 70
IL_ = 65
Raw_Quasi = ReadData_Quasi(λ_max,v,IL_,IA)
p = plot!(Raw_Quasi.λ,Raw_Quasi.P./1000,marker=2,xlabel="λ",ylabel="P (kPa)",label = "λ_max = $(Raw_Quasi.λ_max) v = $(Raw_Quasi.v)", legend=:topleft)

npoints_ = 100
λ_max, Δt, λ, P, times = NoUnloadBuckling_And_Uniform_TimeInc(npoints_,Raw_OneCycle.λ,Raw_OneCycle.P,Raw_OneCycle.Time)
plot(λ,P./1000,markers=2)

plot(Raw_OneCycle.Time)

OneCycleTest_ = OneCycleTest(Raw_OneCycle, npoints_)
plot(OneCycleTest_.λ,OneCycleTest_.σ./1000,marker=2,xlabel="λ",ylabel="P (kPa)",label = "λ_max = $(OneCycleTest_.λ_max) v = $(OneCycleTest_.v)")

λ_max_list = [60,40,20]
v_list = [0.05,0.025,0.01]
Test_Num_list = [3,4,5]
Δt = 0.2
npoints_ = 100
OneCycleTest__list = Vector{OneCycleTest}()
for λ_max in λ_max_list, v in v_list
    Raw_onecycle_tests = Vector{RawOneCycleTest}()
    for Test_Num in Test_Num_list
        push!(Raw_onecycle_tests,ReadData_OneCycle(λ_max,v,IL,IA,Test_Num))
    end
    push!(OneCycleTest__list,OneCycleTest(Raw_onecycle_tests, npoints_))
end
println(OneCycleTest__list)
length(OneCycleTest__list)

plots = []
for OneCycleTest_ in OneCycleTest__list
    push!(plots, plot(OneCycleTest_.λ,OneCycleTest_.σ./1000,marker=2,xlabel="λ",ylabel="P (kPa)",label = false, title = "λ_max = $(OneCycleTest_.λ_max) v = $(OneCycleTest_.v)"))
end
grid_plot = plot(plots..., layout=(3, 3), size=(900, 900))

p = plot()
for OneCycleTest_ in OneCycleTest__list
    p = plot!(OneCycleTest_.λ,OneCycleTest_.σ./1000,marker=2,xlabel="λ",ylabel="P (kPa)",label = "λ_max = $(OneCycleTest_.λ_max) v = $(OneCycleTest_.v)",legend = :outerbottom,size=(800,1000))
end
display(p)

## MARK: Relaxation


RelaxTests = Vector{CreepTest}()
npoints_ = 200
IL_ = 65
λ_max = 60
for λ_max in [60,40,20]
    push!(RelaxTests,CreepTest(ReadData_Relax(λ_max,v,IL_,IA),1000,npoints_))
end
println(RelaxTests)

findmin((ReadData_Relax(λ_max,v,IL_,IA).Time.-1000).^2)
findmax(ReadData_Relax(λ_max,v,IL_,IA).λ)
p = plot()
i = 0
for RelaxTests_ in RelaxTests
    p = plot!(RelaxTests_.t2.+(i*100.0),RelaxTests_.σ2./1000,marker=2,xlabel="Time (s)",ylabel="P (kPa)",label = "λ_max = $(RelaxTests_.λ_max)",legend = :outerbottom,size=(600,600))
    i+=1
end
display(p)

p = plot()
i = 0
for RelaxTests_ in RelaxTests
    p = plot!(RelaxTests_.λ1,RelaxTests_.σ1./1000,marker=2,xlabel="Time (s)",ylabel="P (kPa)",label = "λ_max = $(RelaxTests_.λ_max)",legend = :outerbottom,size=(600,600))
    i+=1
end
display(p)

length(RelaxTests[1].λ2)

## MARK: Quasi-Static


npoints_ = 100
window = 4
v = 0.0001
λ_max = 60
IL_ = 70
Raw_Quasi_list = Vector{RawQuasiTest}()
Raw_Quasi = ReadData_Quasi(λ_max,v,IL_,IA)
p = plot()
p = plot!(Raw_Quasi.λ,Raw_Quasi.P./1000,marker=2,xlabel="λ",ylabel="P (kPa)",label = "λ_max = $(Raw_Quasi.λ_max) v = $(Raw_Quasi.v)", legend=:topleft)
Quasi_test =  QuasiStaticTest(Raw_Quasi, npoints_, window)
push!(Raw_Quasi_list,Raw_Quasi)
length(Quasi_test.λ)
p = plot!(Quasi_test.λ,Quasi_test.σ./1000,marker=2,xlabel="λ",ylabel="P (kPa)",label = "λ_max = $(Raw_Quasi.λ_max) v = $(Raw_Quasi.v)", legend=:topleft)
λ_max = 70
IL_ = 65
Raw_Quasi = ReadData_Quasi(λ_max,v,IL_,IA)
push!(Raw_Quasi_list,Raw_Quasi)
Quasi_test =  QuasiStaticTest(Raw_Quasi, npoints_, window)
p = plot!(Quasi_test.λ,Quasi_test.σ./1000,marker=2,xlabel="λ",ylabel="P (kPa)",label = "λ_max = $(Raw_Quasi.λ_max) v = $(Raw_Quasi.v)", legend=:topleft)

p = plot!(Raw_Quasi.λ,Raw_Quasi.P./1000,marker=2,xlabel="λ",ylabel="P (kPa)",label = "λ_max = $(Raw_Quasi.λ_max) v = $(Raw_Quasi.v)", legend=:topleft)

Quasi_test =  QuasiStaticTest(Raw_Quasi_list, npoints_, window)
p = plot!(Quasi_test.λ,Quasi_test.σ./1000,marker=2,xlabel="λ",ylabel="P (kPa)",label = "λ_max = $(Raw_Quasi.λ_max) v = $(Raw_Quasi.v)", legend=:topleft)

## MARK: Calibration

qs_tests = Vector{QuasiStaticTest}()
push!(qs_tests,Quasi_test)
creep_tests = RelaxTests
onecycle_tests = OneCycleTest__list


build_longterm(C1, C2, C3) = Yeoh3D(λ=0.0, C10=C1, C20=C2, C30=C3)
pn = ["C10",  "C20",  "C30"]  # Parameter names
p0 = [  3e4,   -2e2,    3e0]  # Initial seed
lb = [1.0e3, -6.0e3,  0.0e0]  # Minimum search limits
ub = [3.0e5,  4.0e3,  1.0e3]  # Maximum search limits

opt_func = OptimizationFunction((p,d) -> loss(build_longterm, p, d))
opt_prob = OptimizationProblem(opt_func, p0, qs_tests, lb=lb, ub=ub)
opt_long = solve(opt_prob, ParticleSwarm(lower=lb, upper=ub, n_particles=1000), maxiters=1000, maxtime=60)
opt_prob = OptimizationProblem(opt_func, opt_long.u, qs_tests)
opt_long = solve(opt_prob, Optim.NelderMead(), maxiters=100, maxtime=30)
sol_long = opt_long.u


build_branch(μ, t) = ViscousIncompressible(IncompressibleNeoHookean3D(λ=0.0, μ=μ), τ=exp10(t))
build_branches(p...) = map(splat(build_branch), Iterators.partition(p,2))
build_visco(p...) = GeneralizedMaxwell(build_longterm(sol_long...), build_branches(p...)...)
n_branches = 2
pn = reduce(vcat, ["μ$i", "logτ$i"] for i in 1:n_branches)  # Parameter names
p0 = reduce(vcat, [  1e4,      1.0] for _ in 1:n_branches)  # Initial seed
lb = reduce(vcat, [  1e3,     -1.0] for _ in 1:n_branches)  # Lower search limits
ub = reduce(vcat, [  1e7,      4.0] for _ in 1:n_branches)  # Upper search limits

opt_func = OptimizationFunction((p,d) -> loss(build_visco, p, d))
opt_prob = OptimizationProblem(opt_func, p0, onecycle_tests, lb=lb, ub=ub)
opt_visco = solve(opt_prob, ParticleSwarm(lower=lb, upper=ub, n_particles=10000), maxiters=1000, maxtime=60)
opt_prob = OptimizationProblem(opt_func, opt_visco.u, onecycle_tests)
opt_visco = solve(opt_prob, Optim.NelderMead(), maxiters=100, maxtime=30)
sol_visco = opt_visco.u



plot()
λ_max,v,IL,IA,Test_Num = 20,0.05,68,68*2.9,3
RawEquiBiaxialOneCycleTest_ = ReadData_EquiBiaxialOneCycle(λ_max,v,IL,IA,Test_Num)
plot!(RawEquiBiaxialOneCycleTest_.λ_1,RawEquiBiaxialOneCycleTest_.P_1./1000,marker=2,xlabel="λ",ylabel="P (kPa)",label = "λ_max = $(RawEquiBiaxialOneCycleTest_.λ_max) v = $(RawEquiBiaxialOneCycleTest_.v)")
plot!(RawEquiBiaxialOneCycleTest_.λ_2,RawEquiBiaxialOneCycleTest_.P_2./1000,marker=2,xlabel="λ",ylabel="P (kPa)",label = "λ_max = $(RawEquiBiaxialOneCycleTest_.λ_max) v = $(RawEquiBiaxialOneCycleTest_.v)")

λ_max = 75
RawEquiBiaxialRelaxTest_ = ReadData_EquiBiaxialRelax(λ_max,v,IL,IA)
plot(size=(700,700))
plot!(RawEquiBiaxialRelaxTest_.Time,RawEquiBiaxialRelaxTest_.P_1./1000,marker=2,xlabel="Time (s)",ylabel="P (kPa)",label = "λ_max = $(RawEquiBiaxialRelaxTest_.λ_max) v = $(RawEquiBiaxialRelaxTest_.v)", legend = :outerbottom)
plot!(RawEquiBiaxialRelaxTest_.Time.+50,RawEquiBiaxialRelaxTest_.P_2./1000,marker=2,xlabel="Time (s)",ylabel="P (kPa)",label = "λ_max = $(RawEquiBiaxialRelaxTest_.λ_max) v = $(RawEquiBiaxialRelaxTest_.v)")

## MARK:Biaxial
λ_max,v,IL,IA,Test_Num = 20, 0.05, 68, 68*3, 3
npoints_ = 100

plot()
λ_max,v,IL,IA,Test_Num = 20,0.05,68,68*2.9,3
RawEquiBiaxialOneCycleTest_ = ReadData_EquiBiaxialOneCycle(λ_max,v,IL,IA,Test_Num)
plot!(RawEquiBiaxialOneCycleTest_.λ_1,RawEquiBiaxialOneCycleTest_.P_1./1000,marker=2,xlabel="λ",ylabel="P (kPa)",label = "λ_max = $(RawEquiBiaxialOneCycleTest_.λ_max) v = $(RawEquiBiaxialOneCycleTest_.v)")
plot!(RawEquiBiaxialOneCycleTest_.λ_2,RawEquiBiaxialOneCycleTest_.P_2./1000,marker=2,xlabel="λ",ylabel="P (kPa)",label = "λ_max = $(RawEquiBiaxialOneCycleTest_.λ_max) v = $(RawEquiBiaxialOneCycleTest_.v)")

EquiBiaxialOneCycleTest_ = EquiBiaxialOneCycle(RawEquiBiaxialOneCycleTest_, npoints_)
plot(EquiBiaxialOneCycleTest_.λ_1,EquiBiaxialOneCycleTest_.σ_1./1000,marker=2,xlabel="λ",ylabel="P (kPa)",label = "λ_max = $(EquiBiaxialOneCycleTest_.λ_max) v = $(EquiBiaxialOneCycleTest_.v)")
plot!(EquiBiaxialOneCycleTest_.λ_2,EquiBiaxialOneCycleTest_.σ_2./1000,marker=2,xlabel="λ",ylabel="P (kPa)",label = "λ_max = $(EquiBiaxialOneCycleTest_.λ_max) v = $(EquiBiaxialOneCycleTest_.v)")

scatter(EquiBiaxialOneCycleTest_.λ_1,EquiBiaxialOneCycleTest_.λ_2)

EquiBiaxialOneCycleTest_.λ_1./EquiBiaxialOneCycleTest_.λ_2

build_longterm(μ) = IncompressibleNeoHookean2D(λ=0.0,μ=μ)

build_branch(μ, t) = ViscousIncompressible(IncompressibleNeoHookean3D(λ=0.0, μ=μ), τ=exp10(t))
# build_branch(C1, C2, C3, t) = ViscousIncompressible(Yeoh3D(λ=0.0, C10=C1, C20=C2, C30=C3), τ=exp10(t))
build_branches(p...) = map(splat(build_branch), Iterators.partition(p,2))
build_visco(p...) = GeneralizedMaxwell(build_longterm(sol_long...), build_branches(p...)...)

sol_long = [61347.76647907597]
# Sol visco no equibiaxial relaxation 3 branches
sol_visco = [
    201388.9322396927,
    1.9443069865791058,
    4.418714577779876e6,
    -0.8227763529584671,
    440261.76718093106,
    0.4633834793083379
]

# Sol visco no equibiaxial relaxation 2 branches
sol_visco = [
          1.5644287269978845e6,
     -0.02256351676659682,
 239754.3995972594,
      1.7868023645518827
]

# Sol visco with equibiaxial relaxation 2 branches
sol_visco = [
     668717.4045992043,
      2.3287997978705994,
      5.287834850518206e6,
     -0.6503831837395962
]

long_term = build_longterm(sol_long...)
model = build_visco(sol_visco...)
P1, P2 = evaluate_stress(model, EquiBiaxialOneCycleTest_.Δt, EquiBiaxialOneCycleTest_.λ_1, EquiBiaxialOneCycleTest_.λ_2)

plot!(EquiBiaxialOneCycleTest_.λ_1,P1./1000,marker=2,xlabel="λ",ylabel="P (kPa)",label = "λ_max = $(EquiBiaxialOneCycleTest_.λ_max) v = $(EquiBiaxialOneCycleTest_.v) Pred")
plot!(EquiBiaxialOneCycleTest_.λ_2,P2./1000,marker=2,xlabel="λ",ylabel="P (kPa)",label = "λ_max = $(EquiBiaxialOneCycleTest_.λ_max) v = $(EquiBiaxialOneCycleTest_.v) Pred")

λ_max = 75
RawEquiBiaxialRelaxTest_ = ReadData_EquiBiaxialRelax(λ_max,v,IL,IA)
(RawEquiBiaxialRelaxTest_.λ_1[10]-RawEquiBiaxialRelaxTest_.λ_1[11])/((RawEquiBiaxialRelaxTest_.Time[10]-RawEquiBiaxialRelaxTest_.Time[11]))
plot(size=(700,700))
plot!(RawEquiBiaxialRelaxTest_.Time,RawEquiBiaxialRelaxTest_.P_1./1000,marker=2,xlabel="Time (s)",ylabel="P (kPa)",label = "λ_max = $(RawEquiBiaxialRelaxTest_.λ_max) v = $(RawEquiBiaxialRelaxTest_.v)", legend = :outerbottom)
plot!(RawEquiBiaxialRelaxTest_.Time.+50,RawEquiBiaxialRelaxTest_.P_2./1000,marker=2,xlabel="Time (s)",ylabel="P (kPa)",label = "λ_max = $(RawEquiBiaxialRelaxTest_.λ_max) v = $(RawEquiBiaxialRelaxTest_.v)")

t_max = 1000
npoints_ = 200
EquiBiaxialRelaxTest_ = EquiBiaxialRelaxTest(RawEquiBiaxialRelaxTest_, t_max, npoints_);
plot([EquiBiaxialRelaxTest_.λ_1_1,EquiBiaxialRelaxTest_.λ_2_1],[EquiBiaxialRelaxTest_.σ_1_1./1000,EquiBiaxialRelaxTest_.σ_2_1./1000],marker=2)
plot([[(i-1)*EquiBiaxialRelaxTest_.Δt_2 for i in 1:lastindex(EquiBiaxialRelaxTest_.σ_1_2)],[(i-1)*EquiBiaxialRelaxTest_.Δt_2 for i in 1:lastindex(EquiBiaxialRelaxTest_.σ_2_2)]],[EquiBiaxialRelaxTest_.σ_1_2,EquiBiaxialRelaxTest_.σ_2_2],marker=2)
EquiBiaxialRelaxTest_.λ_1_1[1]

(EquiBiaxialRelaxTest_.λ_max-1.0)/(EquiBiaxialRelaxTest_.Δt_1*EquiBiaxialRelaxTest_.v)

P1,P2,P3,P4 = evaluate_stress(model,EquiBiaxialRelaxTest_.Δt_1, EquiBiaxialRelaxTest_.v,EquiBiaxialRelaxTest_.npoints_load,EquiBiaxialRelaxTest_.Δt_2,EquiBiaxialRelaxTest_.λ_max,EquiBiaxialRelaxTest_.t_max)
plot([EquiBiaxialRelaxTest_.λ_1_1,EquiBiaxialRelaxTest_.λ_2_1,EquiBiaxialRelaxTest_.λ_1_1,EquiBiaxialRelaxTest_.λ_2_1],
[EquiBiaxialRelaxTest_.σ_1_1,EquiBiaxialRelaxTest_.σ_2_1,P1,P2],marker=2)
plot([[(i-1)*EquiBiaxialRelaxTest_.Δt_2 for i in 1:lastindex(EquiBiaxialRelaxTest_.σ_1_2)] for _ in 1:4],[EquiBiaxialRelaxTest_.σ_1_2,EquiBiaxialRelaxTest_.σ_2_2,P3,P4],marker=2)
EquiBiaxialRelaxTest_.λ_1_1[1]

y_true, y_pred = experiment_prediction(model,EquiBiaxialRelaxTest_)

plot(abs.(y_true-y_pred)./maximum(y_true))

length(y_true)
length(y_pred)

## MARK: Relaxation data test

P1,P2 = evaluate_stress(model,RelaxTests[1].Δt1,RelaxTests[1].npoints_load,RelaxTests[1].Δt2,RelaxTests[1].λ_max,1000);
length(P1)
length(P2)

v =((RelaxTests[1].λ_max-1.0)/(RelaxTests[1].npoints_load-1))/RelaxTests[1].Δt1
98*RelaxTests[1].Δt1*v+1.0<RelaxTests[1].λ_max

p = plot()
i = 0
for RelaxTests_ in RelaxTests
    p = plot!(RelaxTests_.t2.+(i*100.0),RelaxTests_.σ2./1000,marker=2,xlabel="Time (s)",ylabel="P (kPa)",label = "λ_max = $(RelaxTests_.λ_max)",legend = :outerbottom,size=(600,600))
    i+=1
end
display(p)

i = 0
for RelaxTests_ in RelaxTests
    P1,P2 = evaluate_stress(model,RelaxTests_.Δt1,RelaxTests_.npoints_load,RelaxTests_.Δt2,RelaxTests_.λ_max,1000);

    p = plot!(RelaxTests_.t2.+(i*100.0),P2./1000,marker=2,xlabel="Time (s)",ylabel="P (kPa)",label = "λ_max = $(RelaxTests_.λ_max)",legend = :outerbottom,size=(600,600))
    i+=1
end
display(p)

p = plot()
i = 0
for RelaxTests_ in RelaxTests
    p = plot!(RelaxTests_.λ1,RelaxTests_.σ1./1000,marker=2,xlabel="Time (s)",ylabel="P (kPa)",label = "λ_max = $(RelaxTests_.λ_max)",legend = :outerbottom,size=(600,600))
    i+=1
end
display(p)

plot!(size=(600,1600))

# i = 0
for RelaxTests_ in RelaxTests
    # P1,P2 = evaluate_stress(model,RelaxTests_.Δt1,RelaxTests_.npoints_load,RelaxTests_.Δt2,RelaxTests_.λ_max,1000);
    P1 = evaluate_stress(model, RelaxTests_.Δt1, RelaxTests_.λ1)
    p = plot!(RelaxTests_.λ1,P1./1000,marker=2,xlabel="Time (s)",ylabel="P (kPa)",label = "λ_max = $(RelaxTests_.λ_max)",legend = :outerbottom,size=(600,600))
    i+=1
end
display(p)