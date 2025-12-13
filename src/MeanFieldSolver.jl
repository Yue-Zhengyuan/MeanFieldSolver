module MeanFieldSolver

export MeanFieldModel, MFUtils, BrillouinZone
export mf_equations, free_energy, solve_mf!

using Printf
using Optim

"""
Abstract type of mean field models
"""
abstract type MeanFieldModel end

include("utils.jl")
include("bz.jl")
include("solve.jl")

include("models/tJCantedAB.jl")
include("models/BCSSpinlessFermion.jl")

end
