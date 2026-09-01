using Printf

abstract type RawExperimentData end

## Read and prepare machine output
struct RawOneCycleTest <: RawExperimentData
  λ_max::Float64
  IL::Float64
  IA::Float64
  v::Float64
  λ::Vector{Float64}
  P::Vector{Float64}
  Time::Vector{Float64}
  Displacement::Vector{Float64}
  Load::Vector{Float64}
end

struct RawRelaxTest <: RawExperimentData
  λ_max::Float64
  IL::Float64
  IA::Float64
  v::Float64
  λ::Vector{Float64}
  P::Vector{Float64}
  Time::Vector{Float64}
  Displacement::Vector{Float64}
  Load::Vector{Float64}
end

struct RawQuasiTest <: RawExperimentData
  λ_max::Float64
  IL::Float64
  IA::Float64
  v::Float64
  λ::Vector{Float64}
  P::Vector{Float64}
  Time::Vector{Float64}
  Displacement::Vector{Float64}
  Load::Vector{Float64}
end

struct RawEquiBiaxialOneCycleTest <: RawExperimentData
  λ_max::Float64
  IL::Float64
  IA::Float64
  v::Float64
  Time::Vector{Float64}
  λ_1::Vector{Float64}
  P_1::Vector{Float64}
  Displacement_1::Vector{Float64}
  Load_1::Vector{Float64}
  λ_2::Vector{Float64}
  P_2::Vector{Float64}
  Displacement_2::Vector{Float64}
  Load_2::Vector{Float64}
end

struct RawEquiBiaxialRelaxTest <: RawExperimentData
  λ_max::Float64
  IL::Float64
  IA::Float64
  v::Float64
  Time::Vector{Float64}
  λ_1::Vector{Float64}
  P_1::Vector{Float64}
  Displacement_1::Vector{Float64}
  Load_1::Vector{Float64}
  λ_2::Vector{Float64}
  P_2::Vector{Float64}
  Displacement_2::Vector{Float64}
  Load_2::Vector{Float64}
end

function ReadData_OneCycle(λ_max,v,IL,IA,Test_Num)
    V = v*IL
    V = round(V,digits=2)
    Data = CSV.File("data/Raw/GoodSample2/First_test__v$(replace(string(v), "." => "_"))_V$(replace(string(V), "." => "_"))_LambdaMax$(λ_max)_Pretested$(Test_Num)Sample.txt") |> Tables.matrix
    Undef_Lenght = IL*1e-3
    Undef_Area = IA*1e-6
    Time = Data[:,1]
    Displacement = Data[:,3].*1e-3
    Load = Data[:,4]
    Time = Time.-Time[1]
    Time = Time./1000.0
    P = Load./Undef_Area
    λ  = (Displacement./Undef_Lenght).+1
    RawOneCycleTest(round(maximum(λ);digits=2),IL,IA,v,λ,P,Time,Displacement,Load)
end

function ReadData_Relax(λ_max,v,IL,IA)
    V = v*IL
    V = round(V,digits=2)
    Data = CSV.File("data/Raw/GoodSample2/First_test__v$(replace(string(v), "." => "_"))_V$(replace(string(V), "." => "_"))_LambdaMax$(λ_max)_Relax_1800_PretestedSample2.txt") |> Tables.matrix
    Undef_Lenght = IL*1e-3
    Undef_Area = IA*1e-6
    Time = Data[:,1]
    Displacement = Data[:,3].*1e-3
    Load = Data[:,4]
    Time = Time.-Time[1]
    Time = Time./1000.0
    P = Load./Undef_Area
    λ  = (Displacement./Undef_Lenght).+1
    RawRelaxTest(round(maximum(λ);digits=2),IL,IA,v,λ,P,Time,Displacement,Load)
end

function ReadData_Quasi(λ_max,v,IL,IA)
    V = v*IL
    V = round(V,digits=4)
    Data = CSV.File("data/Raw/GoodSample2/First_test__v$(replace(string(v), "." => "_"))_V$(replace(string(V), "." => "_"))_LambdaMax$(λ_max)_PretestedSample2.txt") |> Tables.matrix
    Undef_Lenght = IL*1e-3
    Undef_Area = IA*1e-6
    Time = Data[:,1]
    Displacement = Data[:,3].*1e-3
    Load = Data[:,4]
    Time = Time.-Time[1]
    Time = Time./1000.0
    P = Load./Undef_Area
    λ  = (Displacement./Undef_Lenght).+1
    RawQuasiTest(round(maximum(λ);digits=2),IL,IA,v,λ,P,Time,Displacement,Load)
end

function ReadData_EquiBiaxialOneCycle(λ_max,v,IL,IA,Test_Num)
    V = v*IL
    V = round(V,digits=2)
    Data = CSV.File("data/Raw/Biaxial/MaxEB$(λ_max)%_Pre$(Test_Num).csv") |> Tables.matrix
    Undef_Lenght = IL*1e-3
    Undef_Area = IA*1e-6
    Time = Data[:,1]
    Displacement_1 = Data[:,3].*1e-3
    Load_1 = Data[:,2]
    Displacement_2 = Data[:,6].*1e-3
    Load_2 = Data[:,5]
    Time = Time.-Time[1]
    P_1 = Load_1./Undef_Area
    λ_1  = (Displacement_1./Undef_Lenght).+1
    P_2 = Load_2./Undef_Area
    λ_2  = (Displacement_2./Undef_Lenght).+1
    RawEquiBiaxialOneCycleTest(round(maximum(λ_1);digits=1),IL,IA,v,Time,λ_1,P_1,Displacement_1,Load_1,λ_2,P_2,Displacement_2,Load_2)
end

function ReadData_EquiBiaxialRelax(λ_max,v,IL,IA)
    V = v*IL
    V = round(V,digits=2)
    Data = CSV.File("data/Raw/Biaxial/MaxEB2SRe$(λ_max).csv") |> Tables.matrix
    Undef_Lenght = IL*1e-3
    Undef_Area = IA*1e-6
    Time = Data[:,1]
    Displacement_1 = Data[:,3].*1e-3
    Load_1 = Data[:,2]
    Displacement_2 = Data[:,6].*1e-3
    Load_2 = Data[:,5]
    Time = Time.-Time[1]
    P_1 = Load_1./Undef_Area
    λ_1  = (Displacement_1./Undef_Lenght).+1
    P_2 = Load_2./Undef_Area
    λ_2  = (Displacement_2./Undef_Lenght).+1
    RawEquiBiaxialRelaxTest(round(maximum(λ_1);digits=1),IL,IA,v,Time,λ_1,P_1,Displacement_1,Load_1,λ_2,P_2,Displacement_2,Load_2)
end
##
function Base.print(data::Vector{RawOneCycleTest})
  println("Set of $(length(data)) $(RawOneCycleTest)")
  println("__λ__|___v___")
  foreach(r -> @printf(" %.1f | %.3f\n", r.λ_max, r.v), data)
end

function Base.print(data::Vector{RawRelaxTest})
  println("Set of $(length(data)) $(CreepTest)")
  println("__λ__|__w_")
  foreach(r -> @printf(" %.1f | %.1f\n", r.λ_max, r.weight), data)
end

function Base.print(data::Vector{RawQuasiTest})
  println("Set of $(length(data)) $(QuasiStaticTest)")
  println("__w_")
  foreach(r -> @printf(" %.1f\n", r.weight), data)
end

function Base.println(data::Vector{<:RawExperimentData})
  print(data)
  print("\n")
end