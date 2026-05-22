using Plots
using JLD2

Data = load(abspath(dirname(@__FILE__),"20260521/GoodSample2&1_Cyclic_Relax_Quasi_TestData.jld2"))

Data = Data["Test_Data"]
test_list = [3,4]
λ_max_list = [60,40,20]
v_list = [0.05,0.025,0.01]
p = scatter()
for test in test_list
    for λ_max in λ_max_list
        for v in v_list
            p = scatter!(Data["Test = $test ; λ_max = $λ_max ; v = $v"]["λ"],Data["Test = $test ; λ_max = $λ_max ; v = $v"]["P"],xlabel="λ",ylabel="P(Pa)",label="λ_max=$(λ_max/100.0),v=$v")
        end
    end 
end
display(p)

e1 = Data["Test = 3 ; λ_max = 20 ; v = 0.01"]
[keys(e1)...]


p = plot()
P_relax = []
λ_relax = []
for λ_max in λ_max_list
    p = plot!(Data["Test = Relax ; λ_max = $λ_max"]["Time"],Data["Test = Relax ; λ_max = $λ_max"]["P"], label = "λ = $λ_max",ylabel = "P (Pa)", xlabel = "Time (s)")
    loc_relax = argmin((Data["Test = Relax ; λ_max = $λ_max"]["Time"].-1800).^2)
    push!(P_relax,Data["Test = Relax ; λ_max = $λ_max"]["P"][loc_relax])
    push!(λ_relax,λ_max/100.0)
end
display(p)

v = 1e-4
P_max = maximum(Data["Test = Quasi ; v = $v"]["P"])
P_max_loc = argmax(Data["Test = Quasi ; v = $v"]["P"])
λ_P_max = Data["Test = Quasi ; v = $v"]["λ"][P_max_loc]
Load_max = maximum(Data["Test = Quasi ; v = $v"]["Load"])

p = plot(Data["Test = Quasi ; v = $v"]["λ"][[1:P_max_loc...]],Data["Test = Quasi ; v = $v"]["P"][[1:P_max_loc...]],  xlabel = "λ", ylabel = "P", label = "v = $(Data["Test = Quasi ; v = $v"]["v"]) ; P_max = $(round(P_max,digits = 0)) ; λ_max = $(round(λ_P_max,digits = 3))")
p = scatter!(λ_relax,P_relax, label = "Relaxed Stresses")
display(p)
