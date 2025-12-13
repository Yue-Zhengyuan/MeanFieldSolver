"""
Slave fermion mean field t-J model,
on 2D bipartite lattice with canted ansatz.

# Lattice basis vectors (under Cartesian basis)

The nearest neighbor distance is taken as 1.

- Honeycomb lattice: 

# Parameters

- t:    Nearest neighbor hopping coefficient
- J:    Nearest neighbor super-exchange coupling
- β:    1 / temperature (zero-T corresponds to `β = Inf`)
- δ:    number density of holons (hole-doping)
- A:    spinon singlet pairing (AFM)
- B:    spinon hopping (FM)
- D:    holon hopping
- λ:    spinon chemical potential
- dλ:   ±√(λ² - λ₀²), where λ₀ is the minimum (condensation) value of λ
- μ:    holon chemical potential

# Reference

Physical Review B 111, 174518 (arXiv 2301.02274)
"""
@kwdef struct tJCantedAB <: MeanFieldModel
    # Brillouin zone
    bz::BrillouinZone{2} = BrillouinZone((64, 64))
    # nearest unit cell vectors (under lattice basis)
    nbs::Matrix{Float64} = hcat([0.0, 0.0], [0.0, -1.0], [1.0, -1.0])
    # model parameters
    args::Dict{Symbol, Float64}
end

module _tJCantedAB

    using Parameters
    using LinearAlgebra
    using MeanFieldSolver: tJCantedAB
    using MeanFieldSolver.MFUtils

    function gamma(k, model::tJCantedAB)
        return sum(cis(2π * k' * e) for e in eachcol(model.nbs))
    end

    function cal_p(model::tJCantedAB)
        @unpack t, J, B, D = model.args
        return t * D + J * B / 2
    end

    function cal_q(model::tJCantedAB)
        @unpack J, A = model.args
        return J * A
    end

    function cal_r(model::tJCantedAB)
        @unpack t, J, B, D = model.args
        return t * B
    end

    function xi_b(k, model::tJCantedAB)
        p = cal_p(model)
        return p * gamma(k, model)
    end

    function xi_f(k, model::tJCantedAB)
        r = cal_r(model)
        return 2 * r * gamma(k, model)
    end

    function Delta(k, model::tJCantedAB)
        q = cal_q(model)
        return q * gamma(k, model)
    end

    function chi(k, model::tJCantedAB)
        λ = model.args[:λ]
        chi_k2 = λ^2 - abs2(Delta(k, model))
        # if round(chi_k2; digits=14) == 0
        #     chi_k2 = 0
        # end
        return sqrt(chi_k2)
    end

    function Energy_f(k, model::tJCantedAB)
        ξf = abs(xi_f(k, model))
        μ = model.args[:μ]
        return ξf - μ, -ξf - μ
    end

    function Energy_b(k, model::tJCantedAB)
        ξb = abs(xi_b(k, model))
        χ = chi(k, model)
        Eb_A, Eb_B = ξb + χ, -ξb + χ
        # @assert Eb_kB >= 0
        return Eb_A, Eb_B
    end

    """
    Calculate critical value of the spinon chemical potential λ 
    at which the spinons become gapless
    """
    function cal_λ0(model::tJCantedAB)
        z = size(model.nbs, 2)
        λ0 = sqrt(cal_p(model)^2 + cal_q(model)^2) * z
        return λ0
    end

    function free_energy(model::tJCantedAB)
        @unpack t, J, β, δ, A, B, D, μ, dλ = model.args
        z = size(model.nbs, 2)
        N = prod(size(model.bz))
        κ = 1 - δ
        λ = sqrt(cal_λ0(model)^2 + dλ^2)
        fe = 0.0
        if β < Inf && β > 0
            for k in model.bz.ks
                Eb_A, Eb_B = Energy_b(k, model)
                Ef_A, Ef_B = Energy_f(k, model)
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

    function mf_equations(model::tJCantedAB)
        @unpack J, β, δ, A, B, D, μ, dλ = model.args
        z = size(model.nbs, 2)
        N = prod(size(model.bz))
        κ = 1 - δ
        @assert A >= 0
        λ = sqrt(cal_λ0(model)^2 + dλ^2)
        model.args[:λ] = λ
        p, q, r = cal_p(model), cal_q(model), cal_r(model)
        sp, sr = mysgn(p), mysgn(r)
        eqs = zeros(5)
        for k in model.bz.ks
            γ = abs(gamma(k, model))
            χ = chi(k, model)
            Eb_A, Eb_B = Energy_b(k, model)
            nb_A, nb_B = nb(Eb_A * β), nb(Eb_B * β)
            Ef_A, Ef_B = Energy_f(k, model)
            nf_A, nf_B = nf(Ef_A * β), nf(Ef_B * β)
            eqs[1] += (nf_A + nf_B) / (2 * N)
            eqs[2] += (1 / N) * (λ / χ * (1 + nb_A + nb_B) - 1)
            eqs[3] += (A == 0) ? 0 : q * γ^2 / χ * (1 + nb_A + nb_B) / (2 * z * N)
            eqs[4] += (B == 0) ? 0 : (nb_A - nb_B) * γ / (2 * z * N)
            eqs[5] += (D == 0) ? 0 : (nf_A - nf_B) * γ / (2 * z * N)
        end
        eqs -= [δ, κ, A, sp * B, sr * D]
        return norm(eqs), eqs
    end

end

mf_equations(model::tJCantedAB) = _tJCantedAB.mf_equations(model)
free_energy(model::tJCantedAB) = _tJCantedAB.free_energy(model)
