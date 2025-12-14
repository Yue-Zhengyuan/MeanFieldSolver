"""
Mean field self consistency equations for `model`.
"""
function mf_equations(model::MeanFieldModel)
    error("A set of mean field equations should be provided for each specific model.")
end

"""
Mean field free energy per lattice unit cell for `model`.
"""
function free_energy(model::MeanFieldModel)
    error("The expression of free energy should be provided for each specific model.")
end


"""
    solve_mf!(
        model::M, varkeys::Vector{Symbol};
        verbose::Bool = true,
        optimizer = (target, x0) -> Optim.optimize(target, x0)
    ) where {M <: MeanFieldModel}

Mean field self-consistency equation solver. 
The optimization algorithm is specified by the function `optimizer`,
which takes two arguments `target` (the self-consistency equation cost)
and `x0` (initial guess of the mean field solution).
"""
function solve_mf!(
        model::M, varkeys::Vector{Symbol};
        verbose::Bool = true,
        optimizer = (target, x0) -> Optim.optimize(target, x0)
    ) where {M <: MeanFieldModel}
    time0 = time()
    if verbose
        println("Solving self-consistency equations...")
        println("Initial parameters:")
        for k in sort(collect(keys(model.args)))
            println("   ", k => model.args[k])
        end
        println("Variables: ", varkeys)
    end
    x0 = collect(model.args[v] for v in varkeys)

    function target(x)
        for (key, val) in zip(varkeys, x)
            model.args[key] = val
        end
        return first(mf_equations(model))
    end
    result = optimizer(target, x0)
    # solution
    x = Optim.minimizer(result)
    # cost function
    cost = Optim.minimum(result)
    # fill solution into model.args
    for (key, val) in zip(varkeys, x)
        model.args[key] = val
    end
    cost_terms = mf_equations(model)[2]
    # free energy
    fe = free_energy(model)
    # free energy at T = 0
    β0 = model.args[:β]
    model.args[:β] = Inf
    fe0 = free_energy(model)
    # restore temperature
    model.args[:β] = β0
    time1 = time()
    if verbose
        println("----------------")
        println("Optimization finished. Solution:")
        for (k, v) in zip(varkeys, x)
            println("   ", k => v)
        end
        println("               Cost = ", cost)
        println("         Cost terms = ", cost_terms)
        println("        Free energy = ", fe)
        println("Free energy (T = 0) = ", fe0)
    end
    println(@sprintf("Time elapsed: %.3f s", time1 - time0))
    return cost, cost_terms, fe, model
end
