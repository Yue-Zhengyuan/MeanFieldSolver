function solve_mf!(
    varargs::OrderedDict{Symbol,Float64},
    fixargs::OrderedDict{Symbol,Float64},
    model;
    verbose::Bool=true,
)
    if verbose
        println("Solving...")
        println("Variables:\n", varargs)
        println("Fixed parameters:\n", fixargs)
    end
    varkeys = keys(varargs)
    x0 = collect(values(varargs))
    allargs = merge(varargs, fixargs)

    function target(x)
        args = allargs
        for (key, val) in zip(varkeys, x)
            args[key] = val
        end
        return model.mf_equations(args)[2]
    end

    result = Optim.optimize(target, x0, Optim.Options(; g_tol=1e-15))
    # solution 
    x = Optim.minimizer(result)
    # cost function
    cost = Optim.minimum(result)
    # cost terms
    for (key, val) in zip(keys(varargs), x)
        varargs[key] = val
        allargs[key] = val
    end
    cost_terms = model.mf_equations(allargs)[1]
    # free energy
    fe = model.free_energy(allargs)
    allargs2 = deepcopy(allargs)
    allargs2[:β] = Inf
    fe0 = model.free_energy(allargs2)
    if verbose
        println("----------------")
        println("Optimization finished")
        println(varargs)
        println("Cost = ", cost)
        println("Cost terms = ", cost_terms)
        println("Free energy         = ", fe)
        println("Free energy (T = 0) = ", fe0)
    end
    # println()
    return cost, cost_terms, fe
end
