using OrderedCollections
using MeanFieldSolver
using MeanFieldSolver: tJCanted

# fixargs = OrderedDict(pairs((t=0.0, J=1.0, β=800.0, B=0.0, D=0.0, δ=0.0)))
# varkeys = (:A, :μ, :dλ)
# varvals = [0.6049, -0.4453, 0.0016]

fixargs = OrderedDict(pairs((t=5.0, J=1.0, β=1000.0, δ=0.03)))
varkeys = (:A, :B, :D, :μ, :dλ)
varvals = [0.5838 0.05796 -0.02893 -1.617 0.00135]

varargs = OrderedDict(zip(varkeys, varvals))

solve_mf!(varargs, fixargs, tJCanted; verbose=true)
