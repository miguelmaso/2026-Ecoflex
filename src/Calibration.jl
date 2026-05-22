module Calibration

export OneCycleTest
export QuasiStaticTest
export CreepTest

export loss, parallel_loss
export experiment_prediction
export r_squared, stats, uncertainty_models_generator

export plot_experiment_legend!
export plot_experiment!
export plot_experiments
export plot_confidence_bands!
export annotate_r2!

include("ExperimentsData.jl")
include("ConstitutiveModeling.jl")
include("ObjectiveFunctions.jl")
include("ExperimentsPlots.jl")

end
