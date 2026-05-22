using Gridap.TensorValues
using HyperFEM.PhysicalModels, HyperFEM.TensorAlgebra

#region Kinematics

function F_iso(λ::Float64)
  F_vol(λ, 1.0)
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

#endregion
