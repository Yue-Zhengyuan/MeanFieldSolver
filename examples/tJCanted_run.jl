using Parameters
using LinearAlgebra
using MeanFieldSolver
using MeanFieldSolver.MFUtils
import MeanFieldSolver as MFS
include("../models/tJCantedBipartite.jl")

# fixargs = Dict(pairs((t = 0.0, J = 1.0, β = 800.0, B = 0.0, D = 0.0, δ = 0.0)))
# varkeys = [:A, :μ, :dλ]
# varvals = [0.6049, -0.4453, 0.0016]

fixargs = Dict(pairs((t = 5.0, J = 1.0, β = 800.0, δ = 0.04)))
varkeys = [:A, :B, :D, :μ, :dλ]
varvals = [0.5838 0.05796 -0.02893 -1.617 0.00135]

varargs = Dict(zip(varkeys, varvals))
model = tJCantedBipartite(; args = merge(varargs, fixargs))

solve_mf!(model, varkeys; verbose = true)
