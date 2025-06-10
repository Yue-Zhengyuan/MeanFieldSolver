module MFUtils

export theta, log1pexp, log1mexp
export nb, nf, mysgn, myphase

"step function"
theta(x) = (x > 1e-16) ? 1 : 0

"Calculate `log(1 + exp(x))`"
function log1pexp(x)
    return x > 35 ? x : log1p(1 + expm1(x))
end

"Calculate `log(1 - exp(x))`"
function log1mexp(x)
    return log1p(-1 - expm1(x))
end

"Bose distribution function `nb(x) = 1 / (exp(x) - 1)`"
function nb(x)
    if x > 0.0
        return 1 / expm1(x)
    else
        error("Encountered negative/zero beta * E for bosons")
    end
end

"Fermi distribution function `nf(x) = 1 / (exp(x) + 1)`"
function nf(x)
    return 1 / (expm1(x) + 2)
end

"modified sign function (when `x = 0`, return 1 instead of 0)"
function mysgn(x::Real)
    return x >= 0 ? 1 : -1
end

"for any number `x = |x|e^(ia)`, return `e^(ia)`"
function myphase(x::Complex)
    return cis(angle(x))
end

end
