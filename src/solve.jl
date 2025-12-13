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

function solve_mf!(
        model::M, varkeys::Vector{Symbol};
        verbose::Bool = true,
    ) where {M <: MeanFieldModel}
    time0 = time()
    if verbose
        println("Solving...")
        println("Variables:\n", varkeys)
        println("Initial parameters:\n", model.args)
    end
    x0 = collect(model.args[v] for v in varkeys)

    function target(x)
        for (key, val) in zip(varkeys, x)
            model.args[key] = val
        end
        return first(mf_equations(model))
    end

    result = Optim.optimize(target, x0, Optim.Options(; g_tol = 1.0e-15))
    # solution
    x = Optim.minimizer(result)
    # cost function
    cost = Optim.minimum(result)
    # fill solution into model.args
    for (key, val) in zip(varkeys, x)
        model.args[key] = val
    end
    cost_terms = first(mf_equations(model))
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
        println("Cost = ", cost)
        println("Cost terms = ", cost_terms)
        println("Free energy         = ", fe)
        println("Free energy (T = 0) = ", fe0)
    end
    println(@sprintf("Time elapsed: %.3f s", time1 - time0))
    return cost, cost_terms, fe, model
end
