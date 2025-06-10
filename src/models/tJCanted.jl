module tJCanted

using MeanFieldSolver.MFUtils
using Parameters
using LinearAlgebra
using ..MeanFieldSolver: get_BZ_ks
export free_energy, mf_equations

const Ns = [64, 64]
const N = prod(Ns)
# nearest neighbor unit cell vectors (lattice basis) on honeycomb lattice
const del = hcat([0.0, 0.0], [0.0, -1.0], [1.0, -1.0])
const z = size(del, 2)
const ks = get_BZ_ks(Ns)

function gamma(k)
    return sum(cis(2π * k' * e) for e in eachcol(del))
end

function cal_p(args)
    @unpack t, J, B, D = args
    return t * D + J * B / 2.0
end

function cal_q(args)
    @unpack J, A = args
    return J * A
end

function cal_r(args)
    @unpack t, J, B, D = args
    return t * B
end

function xi_b(k, args)
    p = cal_p(args)
    return p * gamma(k)
end

function xi_f(k, args)
    r = cal_r(args)
    return 2 * r * gamma(k)
end

function Delta(k, args)
    q = cal_q(args)
    return q * gamma(k)
end

function chi(k, args)
    @unpack λ = args
    chi_k2 = λ^2 - abs2(Delta(k, args))
    # if round(chi_k2; digits=14) == 0
    #     chi_k2 = 0
    # end
    return sqrt(chi_k2)
end

function Energy_f(k, args)
    ξf = abs(xi_f(k, args))
    @unpack μ = args
    return ξf - μ, -ξf - μ
end

function Energy_b(k, args)
    ξb = abs(xi_b(k, args))
    χ = chi(k, args)
    Eb_A, Eb_B = ξb + χ, -ξb + χ
    # @assert Eb_kB >= 0
    return Eb_A, Eb_B
end

"""
Calculate critical value of the spinon chemical potential λ 
at which the spinons become gapless
"""
function cal_λ0(args)
    λ0 = sqrt(cal_p(args)^2 + cal_q(args)^2) * z
    return λ0
end

"Free energy per site"
function free_energy(args)
    @unpack t, J, β, δ, A, B, D, μ, dλ = args
    κ = 1 - δ
    λ = sqrt(cal_λ0(args)^2 + dλ^2)
    fe = 0.0
    if β < Inf && β > 0
        for k in ks
            Eb_A, Eb_B = Energy_b(k, args)
            Ef_A, Ef_B = Energy_f(k, args)
            fe +=
                (1 / β / N) * (log1mexp(-β * Eb_A) + log1mexp(-β * Eb_B)) +
                (Eb_A + Eb_B - 2 * λ) / (2 * N)
            fe -= (0.5 / β / N) * (log1pexp(-β * Ef_A) + log1pexp(-β * Ef_B))
        end
        fe += z * (-2 * t * B * D + (J / 2) * (2 * A^2 - B^2))
        fe += μ * δ - λ * κ
        fe += -(1 / 8) * z * J * κ^2
    elseif β == Inf
        fe = z * (2 * t * B * D + J * (B^2 / 2 - A^2 - κ^2 / 8))
    end
    return fe
end

"""
Self-consistecy equations for t-J canted state

# Parameters

- t, J: tJ model parameters
- β:    1 / temperature (zero-T corresponds to `β = Inf`)
- δ:    number density of holons
- A:    spinon singlet pairing (AFM)
- B:    spinon hopping (FM)
- D:    holon hopping
- λ:    spinon chemical potential
- μ:    holon chemical potential
"""
function mf_equations(args)
    @unpack J, β, δ, A, B, D, μ, dλ = args
    κ = 1 - δ
    @assert A >= 0
    λ = sqrt(cal_λ0(args)^2 + dλ^2)
    args[:λ] = λ
    p, q, r = cal_p(args), cal_q(args), cal_r(args)
    sp, sr = mysgn(p), mysgn(r)
    eqs = zeros(5)
    for k in ks
        γ = abs(gamma(k))
        χ = chi(k, args)
        Eb_A, Eb_B = Energy_b(k, args)
        nb_A, nb_B = nb(Eb_A * β), nb(Eb_B * β)
        Ef_A, Ef_B = Energy_f(k, args)
        nf_A, nf_B = nf(Ef_A * β), nf(Ef_B * β)
        eqs[1] += (nf_A + nf_B) / (2 * N)
        eqs[2] += (1 / N) * (λ / χ * (1 + nb_A + nb_B) - 1)
        eqs[3] += (A == 0) ? 0 : q * γ^2 / χ * (1 + nb_A + nb_B) / (2 * z * N)
        eqs[4] += (B == 0) ? 0 : (nb_A - nb_B) * γ / (2 * z * N)
        eqs[5] += (D == 0) ? 0 : (nf_A - nf_B) * γ / (2 * z * N)
    end
    eqs -= [δ, κ, A, sp * B, sr * D]
    return eqs, norm(eqs)
end

end
