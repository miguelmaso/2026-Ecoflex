using Printf
using Plots
using Optim

# Your data variables
mu_index = 1
val = 4.7e05
err = 1.1e03
pct = 0.2
total = 14540.5

# The template string
# \u03BC is the Unicode for the Greek letter mu (μ)
@printf("\\(\\mu\\)%d    &  %.1e \\(\\pm\\) %.1e (  %.1f%%) & %.1f\n", 
        mu_index, val, err, pct, total)

## MARK: Density of tensile sample with excess material in the border

W = 11.3e-3
T = 2.43e-3
L = 103.5e-3
V = W*T*L

ΔV = 1.1e-3*(14.7e-3-W)*(0.25*L)

m = 3.71e-3

ρ = m/(V+ΔV)

## MARK: Density of disk with defect (less material) in the surface

D = 40e-3
T = 2.4e-3

A = (pi*0.25)*D^2
V = A*T

#volume reduction due to defects
V = 0.97*V

m = 3.37e-3
ρ = m/(V)


## MARK: HV Power Supply behavior

data_time = [0.7156959526159921,
1.1846001974333658,
4.763079960513325,
9.846989141164855]

data_capacitance = [
        235e-12,
        500e-12,
        2.35e-9,
        5e-9
]

data_capacitance = data_capacitance*1.0e12

plot(data_capacitance, data_time,marker=2)

f(x) = (m,b) -> m*x + b

Loss(x) = sum((data_time - map(X -> f(X)(x[1],x[2]),data_capacitance)).^2)
Loss([data_time[1]/data_capacitance[1],0])
sol = Optim.optimize(Loss,[data_time[1]/data_capacitance[1],0.0],Optim.NelderMead())
sol.minimizer
sol.minimum

f(17.0)(sol.minimizer...)