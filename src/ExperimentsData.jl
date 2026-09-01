using DataInterpolations
using Loess
using Printf

abstract type ExperimentData end

#region Structs

"""
One-cycle loading-unloading test.
The constructor takes a dictionary with keys 'Time', 'λ' and 'P',
it perform some sanity checks, remove non-significant points and
generates a uniform time-step series.
"""
struct OneCycleTest <: ExperimentData
  v::Float64
  Δt::Float64
  λ_max::Float64
  λ::Vector{Float64}
  σ::Vector{Float64}
  weight::Float64
end

function OneCycleTest(Raw_OneCycle::RawOneCycleTest, npoints_, weight=1.0)

  # Get the peak stress to copmute the loading rate
  _, max_pos = findmax(Raw_OneCycle.P)
  vel = (Raw_OneCycle.λ[max_pos-1]-Raw_OneCycle.λ[1]) / Raw_OneCycle.Time[max_pos-1]
  vel_tol = 0.005
  vel = round(vel / vel_tol) * vel_tol

  λ_max, Δt, λ, P, times = NoUnloadBuckling_And_Uniform_TimeInc(npoints_,Raw_OneCycle.λ,Raw_OneCycle.P,Raw_OneCycle.Time)
  OneCycleTest(vel, Δt, λ_max, λ, P, weight)
end

function OneCycleTest(Raw_OneCycles::Vector{RawOneCycleTest}, npoints_, weight=1.0)

  # Get the peak stress to copmute the loading rate
  v_list = []
  for Raw_OneCycle in Raw_OneCycles
    push!(v_list,Raw_OneCycle.v)
  end
  if all(x -> x == v_list[1], v_list)
  else
    error("Data does not com from the same strain rate test")
  end
  λ_max_list = []
  for Raw_OneCycle in Raw_OneCycles
    push!(λ_max_list,Raw_OneCycle.λ_max)
  end
  if all(x -> x == λ_max_list[1], λ_max_list)
  else
    error("Data does not com from the same λ_max test")
  end

   vel_, λ_max_, Δt_, λ_, P_, times_ = [], [], [], [], [], []
  for Raw_OneCycle in Raw_OneCycles
    _, max_pos = findmax(Raw_OneCycle.λ)
    vel = (Raw_OneCycle.λ[max_pos-1]-Raw_OneCycle.λ[1]) / Raw_OneCycle.Time[max_pos-1]
    vel_tol = 0.005
    vel = round(vel / vel_tol) * vel_tol
    λ_max, Δt, λ, P, times = NoUnloadBuckling_And_Uniform_TimeInc(npoints_,Raw_OneCycle.λ,Raw_OneCycle.P,Raw_OneCycle.Time)
    push!(vel_,vel)
    push!(λ_max_,λ_max)
    push!(Δt_,Δt)
    push!(λ_,λ)
    push!(P_,P)
    push!(times_,times)
  end
  vel = round(sum(vel_)/length(vel_),digits = 3)
  λ_max = round(sum(λ_max_)/length(λ_max_), digits = 2)
  Δt = sum(Δt_)/length(Δt_)
  λ = sum(λ_)./length(λ_)
  P = sum(P_)./length(P_)
  times = sum(times_)./length(times_)

  OneCycleTest(vel, Δt, λ_max, λ, P, weight)
end

"""
Creep or relaxation test.
The constructor takes a dictionary with keys 'Time', 'λ' and 'P'
"""
struct CreepTest <: ExperimentData
  λ_max::Float64
  Δt1::Float64
  npoints_load::Int32
  λ1::Vector{Float64}
  σ1::Vector{Float64}
  t1::Vector{Float64}
  Δt2::Float64
  λ2::Vector{Float64}
  σ2::Vector{Float64}
  t2::Vector{Float64}
  weight::Float64
end

function CreepTest(Raw_Relax::RawRelaxTest, t_max, npoints_, weight=1.0)
  _, maxpos = findmax(Raw_Relax.P)
  npoints_load = Int32(round(npoints_/2,digits=0))
  λ_max, Δt1, λ1, P1, times1 = NoUnloadBuckling_And_Uniform_TimeInc(npoints_load,Raw_Relax.λ[1:maxpos],Raw_Relax.P[1:maxpos],Raw_Relax.Time[1:maxpos])
  _, min_pos = findmin((Raw_Relax.Time.-t_max).^2)
  _, Δt2, λ2, P2, times2 = NoUnloadBuckling_And_Uniform_TimeInc(Int32(round(npoints_/2,digits=0)),Raw_Relax.λ[maxpos:min_pos],Raw_Relax.P[maxpos:min_pos],Raw_Relax.Time[maxpos:min_pos])
  println(maximum(λ1))
  CreepTest(maximum(λ1), Δt1, npoints_load, λ1, P1, times1, Δt2, λ2, P2, times2, weight)
end


"""
Quasi-static test.
The constructor takes a dictionary with keys 'λ' and 'P'.
"""
struct QuasiStaticTest <: ExperimentData
  λ::Vector{Float64}
  σ::Vector{Float64}
  weight::Float64
end

function QuasiStaticTest(Raw_Quasi::RawQuasiTest, npoints_, window, weight=1.0)
  λ_max, Δt, λ, P, times = OnlyLoad_OnlyIncreasingDispl_MovingAvg_UniformTimeInc(npoints_,window,Raw_Quasi.λ,Raw_Quasi.P,Raw_Quasi.Time)

  # Generate the new struct
  QuasiStaticTest(λ, P, weight)
end

function QuasiStaticTest(Raw_Quasis::Vector{RawQuasiTest}, npoints_, window, weight=1.0)
  λ_max_, Δt_, λ_, P_, times_ = [], [], [], [], [], []
  for Raw_Quasi in Raw_Quasis
    λ_max, Δt, λ, P, times = OnlyLoad_OnlyIncreasingDispl_MovingAvg_UniformTimeInc(npoints_,window,Raw_Quasi.λ,Raw_Quasi.P,Raw_Quasi.Time)
    push!(λ_max_,λ_max)
    push!(Δt_,Δt)
    push!(λ_,λ)
    push!(P_,P)
    push!(times_,times)
  end
  λ_max = round(sum(λ_max_)/length(λ_max_), digits = 2)
  Δt = sum(Δt_)/length(Δt_)
  λ = sum(λ_)./length(λ_)
  P = sum(P_)./length(P_)
  times = sum(times_)./length(times_)

  # Generate the new struct
  QuasiStaticTest(λ, P, weight)
end

"""
EquiBiaxialOneCycle
The constructor takes a dictionary of   λ_max, v, Δt, λ_1, P_1, λ_2, P_2
"""

struct EquiBiaxialOneCycle <: ExperimentData
  λ_max::Float64
  v::Float64
  Δt::Float64
  λ_1::Vector{Float64}
  σ_1::Vector{Float64}
  λ_2::Vector{Float64}
  σ_2::Vector{Float64}
  weight::Float64
end

function EquiBiaxialOneCycle(Raw_EquiBiaxialOneCycleTest::RawEquiBiaxialOneCycleTest, npoints_, weight=1.0)

  # Get the peak stress to copmute the loading rate
  _, max_pos = findmax(Raw_EquiBiaxialOneCycleTest.P_1)
  vel_1 = (Raw_EquiBiaxialOneCycleTest.λ_1[max_pos-1]-Raw_EquiBiaxialOneCycleTest.λ_1[1]) / Raw_EquiBiaxialOneCycleTest.Time[max_pos-1]
  _, max_pos = findmax(Raw_EquiBiaxialOneCycleTest.P_2)
  vel_2 = (Raw_EquiBiaxialOneCycleTest.λ_2[max_pos-1]-Raw_EquiBiaxialOneCycleTest.λ_2[1]) / Raw_EquiBiaxialOneCycleTest.Time[max_pos-1]
  vel_tol = 0.005
  vel = (vel_1+vel_2)/2
  vel = round(vel / vel_tol) * vel_tol

  λ_max, Δt, λ_1, P_1, _ = NoUnloadBuckling_And_Uniform_TimeInc(npoints_,Raw_EquiBiaxialOneCycleTest.λ_1,Raw_EquiBiaxialOneCycleTest.P_1,Raw_EquiBiaxialOneCycleTest.Time)
  λ_max, Δt, λ_2, P_2, _ = NoUnloadBuckling_And_Uniform_TimeInc(npoints_,Raw_EquiBiaxialOneCycleTest.λ_2,Raw_EquiBiaxialOneCycleTest.P_2,Raw_EquiBiaxialOneCycleTest.Time)
  EquiBiaxialOneCycle(λ_max, vel, Δt, λ_1, P_1, λ_2, P_2, weight)
end

"""
"""

struct EquiBiaxialRelaxTest <: ExperimentData
  λ_max::Float64
  Δt_1::Float64
  npoints_load::Int32
  v::Float64
  λ_1_1::Vector{Float64}
  σ_1_1::Vector{Float64}
  λ_2_1::Vector{Float64}
  σ_2_1::Vector{Float64}
  Δt_2::Float64
  λ_1_2::Vector{Float64}
  σ_1_2::Vector{Float64}
  λ_2_2::Vector{Float64}
  σ_2_2::Vector{Float64}
  t_max::Float64
  weight::Float64
end

function EquiBiaxialRelaxTest(RawEquiBiaxialRelaxTest_::RawEquiBiaxialRelaxTest, t_max, npoints_, weight=1.0)
  _, min_pos = findmin((RawEquiBiaxialRelaxTest_.Time.-t_max).^2)
  Relax1 = CreepTest(RawRelaxTest(RawEquiBiaxialRelaxTest_.λ_max,RawEquiBiaxialRelaxTest_.IL,RawEquiBiaxialRelaxTest_.IA,RawEquiBiaxialRelaxTest_.v,RawEquiBiaxialRelaxTest_.λ_1[1:min_pos],RawEquiBiaxialRelaxTest_.P_1[1:min_pos],RawEquiBiaxialRelaxTest_.Time[1:min_pos],RawEquiBiaxialRelaxTest_.Displacement_1[1:min_pos],RawEquiBiaxialRelaxTest_.Load_1[1:min_pos]), t_max, npoints_, weight)
  # v1 = (Relax1.λ1[11]-Relax1.λ1[1])/(Relax1.t1[11]-Relax1.t1[1])

  Relax2 = CreepTest(RawRelaxTest(RawEquiBiaxialRelaxTest_.λ_max,RawEquiBiaxialRelaxTest_.IL,RawEquiBiaxialRelaxTest_.IA,RawEquiBiaxialRelaxTest_.v,RawEquiBiaxialRelaxTest_.λ_2[1:min_pos],RawEquiBiaxialRelaxTest_.P_2[1:min_pos],RawEquiBiaxialRelaxTest_.Time[1:min_pos],RawEquiBiaxialRelaxTest_.Displacement_2[1:min_pos],RawEquiBiaxialRelaxTest_.Load_2[1:min_pos]), t_max, npoints_, weight)
  # v2 = (Relax2.λ1[11]-Relax2.λ1[1])/(Relax2.t1[11]-Relax2.t1[1])
  # v = (v1+v2)/2
  # println((v1,v2,v))
  EquiBiaxialRelaxTest(Relax1.λ_max,Relax1.Δt1,Relax1.npoints_load,v,Relax1.λ1,Relax1.σ1,Relax2.λ1,Relax2.σ1,Relax1.Δt2,Relax1.λ2,Relax1.σ2,Relax2.λ2,Relax2.σ2,t_max,weight)
end
#endregion
#region Functions

npoints(test::OneCycleTest) = length(test.λ)

npoints(test::CreepTest) = length(test.t)

npoints(test::QuasiStaticTest) = length(test.λ)

npoints(test::EquiBiaxialRelaxTest) = length(test.σ_1_1) + length(test.σ_2_1) + length(test.σ_1_2) + length(test.σ_2_2)


npoints(tests::Vector{<:ExperimentData}) = sum(npoints, tests)

Base.maximum(test::ExperimentData) = maximum(test.σ)

Base.maximum(test::EquiBiaxialRelaxTest) = maximum(test.σ_1_1)


function Base.print(data::Vector{OneCycleTest})
  println("Set of $(length(data)) $(OneCycleTest)")
  println("__λ__|___v___|__w_")
  foreach(r -> @printf(" %.1f | %.3f | %.1f\n", r.λ_max, r.v, r.weight), data)
end

function Base.print(data::Vector{CreepTest})
  println("Set of $(length(data)) $(CreepTest)")
  println("__λ__|__w_")
  foreach(r -> @printf(" %.1f | %.1f\n", r.λ_max, r.weight), data)
end

function Base.print(data::Vector{QuasiStaticTest})
  println("Set of $(length(data)) $(QuasiStaticTest)")
  println("__w_")
  foreach(r -> @printf(" %.1f\n", r.weight), data)
end

function Base.println(data::Vector{<:ExperimentData})
  print(data)
  print("\n")
end



getfirst(pred,itr) = first(Iterators.filter(pred,itr))

#endregion
