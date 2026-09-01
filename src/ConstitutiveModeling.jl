using Gridap.TensorValues
using HyperFEM.PhysicalModels, HyperFEM.TensorAlgebra

#region Kinematics

function F_iso(λ::Float64)
  F_vol(λ, 1.0)
end

function F_iso(λ_1::Float64,λ_2::Float64)
  F_vol(λ_1,λ_2, 1.0)
end

function F_vol(J::Float64)
  λ = J^(-1/3)
  TensorValue(λ, 0, 0, 0, λ, 0, 0, 0, λ)
end

function F_vol(λ::Float64, J::Float64)
  TensorValue(λ, 0, 0, 0, λ^(-1/2), 0, 0, 0, λ^(-1/2)) .* J^(1/3)
end

function F_vol(λ::Float64, λ2::Float64, J::Float64)
  TensorValue(λ, 0, 0, 0, λ2, 0, 0, 0, J/(λ*λ2))
end

#endregion
#region Dynamics

function new_state(model::ViscoElastic, F, Fn, A...)
  map(model.branches, A) do b, Ai
    _, Se, ∂Se∂Ce = SecondPiola(b.elasto)
    HyperFEM.PhysicalModels.ReturnMapping(b, Se, ∂Se∂Ce, F, Fn, Ai)[2]
  end
end

function evaluate_stress(model::Elasto, λ_values)
  P_func = model()[2]
  map(λ_values) do λ
    F = F_iso(λ)
    P = P_func(F)
    p = P[2,2] * F[2,2]  # Volumetric pressure term
    return P[1] - p / F[1]
  end
end

function evaluate_stress(model::ViscoElastic, Δt, λ_values)
  update_time_step!(model, Δt)
  P_func = model()[2]
  n  = length(model.branches)
  A  = ntuple(_ -> VectorValue(I3..., 0.0), Val(n))
  Fn = F_iso(1.0)
  map(λ_values) do λ
    F = F_iso(λ)
    P = try P_func(F, Fn, A...) catch; zeros(3,3) end
    A = try new_state(model, F, Fn, A...) catch; A end
    Fn = F
    p = P[2,2] * F[2,2]  # Volumetric pressure term
    return P[1] - p / F[1]
  end
end

function evaluate_stress(model::ViscoElastic, Δt, λ_1_values, λ_2_values)
  update_time_step!(model, Δt)
  P_func = model()[2]
  n  = length(model.branches)
  A  = ntuple(_ -> VectorValue(I3..., 0.0), Val(n))
  Fn = F_iso(1.0)
  P = map(λ_1_values,λ_2_values) do λ_1,λ_2
    F = F_iso(λ_1,λ_2)
    P = try P_func(F, Fn, A...) catch; zeros(3,3) end
    A = try new_state(model, F, Fn, A...) catch; A end
    Fn = F
    p = P[3,3] * F[3,3]  # Volumetric pressure term
    return [Float64(P[1,1] - p / F[1,1]), Float64(P[2,2] - p / F[2,2])]
  end
  P = reduce(hcat,P)
  return P[1,:],P[2,:]
end

function evaluate_stress(model::ViscoElastic,Δt_1,v,npoints_load,Δt_2,λ_max,t_max)
  t = 0
  update_time_step!(model, Δt_1)
  P_func = model()[2]
  n  = length(model.branches)
  A  = ntuple(_ -> VectorValue(I3..., 0.0), Val(n))
  λ = 1.0
  Fn = F_iso(λ,λ)
  F = F_iso(λ,λ)
  P1 = Vector{Float64}()
  P2 = Vector{Float64}()
  P3 = Vector{Float64}()
  P4 = Vector{Float64}()
  P = P_func(F, Fn, A...)
  p = P[3,3] * F[3,3]
  push!(P1,P[1,1] - p / F[1,1])
  push!(P2,P[2,2] - p / F[2,2])
  v = ((λ_max-1.0)/(npoints_load-1))/Δt_1
  while t<=t_max
    if λ<λ_max
      t+=Δt_1
      λ+=Δt_1*v
      F = F_iso(λ,λ)
      P = try P_func(F, Fn, A...) catch; zeros(3,3) end
      A = try new_state(model, F, Fn, A...) catch; A end
      Fn = F
      p = P[3,3] * F[3,3]
      push!(P1,P[1,1] - p / F[1,1])
      push!(P2,P[2,2] - p / F[2,2])
      if λ==λ_max
        update_time_step!(model, Δt_2)
      end
    else
      t+=Δt_2
      F = F_iso(λ,λ)
      P = try P_func(F, Fn, A...) catch; zeros(3,3) end
      A = try new_state(model, F, Fn, A...) catch; A end
      Fn = F
      p = P[3,3] * F[3,3]
      push!(P3,P[1,1] - p / F[1,1])
      push!(P4,P[2,2] - p / F[2,2])
    end
  end
  return P1,P2,P3,P4
end

function evaluate_stress(model::ViscoElastic,Δt_1,npoints_load,Δt_2,λ_max,t_max)
  t = 0
  update_time_step!(model, Δt_1)
  P_func = model()[2]
  n  = length(model.branches)
  A  = ntuple(_ -> VectorValue(I3..., 0.0), Val(n))
  λ = 1.0
  Fn = F_iso(λ)
  F = F_iso(λ)
  P1 = Vector{Float64}()
  P2 = Vector{Float64}()
  P = P_func(F, Fn, A...)
  p = P[3,3] * F[3,3]
  push!(P1,P[1,1] - p / F[1,1])
  v = ((λ_max-1.0)/(npoints_load-1))/Δt_1
  println(v)
  i = 0
  while t<=t_max
    if λ<λ_max-Δt_1*v
      i += 1
      t+=Δt_1
      λ+=Δt_1*v
      F = F_iso(λ)
      P = try P_func(F, Fn, A...) catch; zeros(3,3) end
      A = try new_state(model, F, Fn, A...) catch; A end
      Fn = F
      p = P[3,3] * F[3,3]
      push!(P1,P[1,1] - p / F[1,1])
      if λ≈λ_max
        println(λ)
        println(i)
        println(length(P1))
        update_time_step!(model, Δt_2)
      end
    else
      t+=Δt_2
      F = F_iso(λ)
      P = try P_func(F, Fn, A...) catch; zeros(3,3) end
      A = try new_state(model, F, Fn, A...) catch; A end
      Fn = F
      p = P[3,3] * F[3,3]
      push!(P2,P[1,1] - p / F[1,1])
    end
  end
  return P1,P2
end

#endregion
