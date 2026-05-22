using JLD2

include("src/Calibration.jl")
using .Calibration

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

print(onecycle_tests)
print(creep_tests)
print(qs_tests)
