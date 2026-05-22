using DataInterpolations
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

function OneCycleTest(df, weight=1.0)
  raw_times = df["Time"]
  raw_λ = df["λ"]
  raw_σ = df["P"]

  # Sanity check
  @assert length(raw_times) == length(raw_λ) == length(raw_σ)
  raw_times .-= raw_times[1]

  # Get the peak stress to copmute the loading rate
  σ_max, max_pos = findmax(raw_σ)
  vel = raw_λ[max_pos-1] / raw_times[max_pos-1]

  # Remove the near zero/negative stresses after unloading
  last_measure = findlast(σ -> σ > 0.1*σ_max, raw_σ)
  last_time = raw_times[last_measure]

  # Generate a new time series and interpolate data
  npoints = 100
  times = range(0, last_time, length=npoints)
  interp_λ = LinearInterpolation(raw_λ, raw_times)
  interp_σ = LinearInterpolation(raw_σ, raw_times)
  λ = interp_λ(times)
  σ = interp_σ(times)

  # Generate the new struct
  λ_max = round(maximum(λ); digits=1)
  Δt = times[2] - times[1]
  OneCycleTest(vel, Δt, λ_max, λ, σ, weight)
end


"""
Creep or relaxation test.
The constructor takes a dictionary with keys 'Time', 'λ' and 'P'
"""
struct CreepTest <: ExperimentData
  Δt::Float64
  λ_max::Float64
  t::Vector{Float64}
  λ::Vector{Float64}
  σ::Vector{Float64}
  weight::Float64
end

function CreepTest(df, weight=1.0)
  raw_times = df["Time"]
  raw_λ = df["λ"]
  raw_σ = df["P"]

  # Sanity check
  @assert length(raw_times) == length(raw_λ) == length(raw_σ)
  raw_times .-= raw_times[1]

  # Get the maximum stress and remove near-zero measures
  σ_max = maximum(raw_σ)
  last_measure = findlast(σ -> σ > 0.1*σ_max, raw_σ)
  last_time = raw_times[last_measure]

  # Generate a new time series and interpolate data
  npoints = 100
  times = range(0, last_time, length=npoints)
  interp_λ = LinearInterpolation(raw_λ, raw_times)
  interp_σ = LinearInterpolation(raw_σ, raw_times)
  λ = interp_λ(times)
  σ = interp_σ(times)

  # Generate the new struct
  λ_max = round(maximum(λ); digits=1)
  Δt = times[2] - times[1]
  CreepTest(Δt, λ_max, times, λ, σ, weight)
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

function QuasiStaticTest(df, weight=1.0)
  raw_λ = df["λ"]
  raw_σ = df["P"]

  @assert length(raw_λ) == length(raw_σ)

  QuasiStaticTest(raw_λ, raw_σ, weight)
end

#endregion
#region Functions

npoints(test::OneCycleTest) = length(test.λ)

npoints(test::CreepTest) = length(test.t)

npoints(test::QuasiStaticTest) = length(test.λ)

Base.maximum(test::ExperimentData) = maximum(test.σ)

function Base.print(data::Vector{OneCycleTest})
  println("Set of $(length(data)) $(OneCycleTest)")
  println("__λ__|___v__|__w_")
  foreach(r -> @printf(" %.1f | %.2f | %.1f\n", r.λ_max, r.v, r.weight), data)
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
