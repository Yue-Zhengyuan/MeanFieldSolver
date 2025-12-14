"""
BCS mean field spinless fermion model on Bravais Lattice

# Parameters in args

- t:    Nearest neighbor hopping coefficient
- g:    Nearest neighbor attraction
- β:    1 / temperature (zero-T corresponds to `β = Inf`)
- x (or x1re, x2re, x2im, ...):
        independent components of nearest neighbor fermion pairing
        (imag(Δ1) is assumed to be 0)
- μ:    chemical potential
- δ:    number density of fermions 
        (should be determined self-consistently from other parameters)

# Reference

Physical Review B 81, 024504 (arXiv 0908.2805)
"""
struct BCSSpinlessFermion{D} <: MeanFieldModel
    # lattice name
    lattice::Symbol
    # pairing name
    pairing::Symbol
    # Brillouin zone
    bz::BrillouinZone{D}
    # nearest neighbor vectors (under lattice basis)
    nbs::Matrix{Float64}
    # model parameters
    args::Dict{Symbol, Float64}
    # force nonzero pairing
    force_pairing::Bool
end

const spinless_fermion_pairings = Dict(
    :triangle => [:none, :f, :pip],
    :square => [:none, :pip, :px, :py, :ppp, :pmp]
)

function BCSSpinlessFermion(
        lattice::Symbol, pairing::Symbol, bz::BrillouinZone, args::Dict{Symbol}; force_pairing::Bool = false
    )
    @assert lattice in (:triangle, :square) "Unrecongnized lattice."
    @assert pairing in spinless_fermion_pairings[lattice] "Unrecongnized pairing."
    # lattice basis vectors
    basis = if lattice == :triangle
        [[1.0; 0.0];; [-0.5; sqrt(3) / 2]]
    else # default: :square
        [[1.0; 0.0];; [0.0; 1.0]]
    end
    nbs = if lattice == :triangle
        [[1.0; 0.0];; [-0.5; sqrt(3) / 2];; [-0.5; -sqrt(3) / 2]]
    else
        [[1.0; 0.0];; [0.0; 1.0]]
    end
    nbs = inv(basis) * nbs
    return BCSSpinlessFermion(lattice, pairing, bz, nbs, args, force_pairing)
end

module _BCSSpinless

    using Parameters
    using LinearAlgebra
    using MeanFieldSolver: BCSSpinlessFermion, BrillouinZone, sum_bz
    using MeanFieldSolver: spinless_fermion_pairings
    using MeanFieldSolver.MFUtils

    """
    For different pairing patterns, 
    construct `Δs` from the independent components
    """
    function get_Δs(model::BCSSpinlessFermion)
        # number of neighbor vectors
        Δs, pairing = nothing, model.pairing
        if model.lattice == :square
            if pairing == :px
                x = model.args[:x]
                Δs = [x, 0.0]
            elseif pairing == :py
                x = model.args[:x]
                Δs = [0.0, x]
            elseif pairing == :ppp
                x = model.args[:x]
                Δs = [x, x]
            elseif pairing == :pmp
                x = model.args[:x]
                Δs = [x, -x]
            elseif pairing == :pip
                x = model.args[:x]
                Δs = [x, im * x]
            else
                @unpack x1re, x2re, x2im = model.args
                Δs = [x1re, x2re + 1.0im * x2im]
            end
        elseif model.lattice == :triangle
            if pairing == "f"
                x = model.args[:x]
                Δs = fill(x, 3)
            elseif pairing == :pip
                x = model.args[:x]
                Δs = x * [1.0, cispi(2 / 3), cispi(4 / 3)]
            else
                @unpack x1re, x2re, x2im, x3re, x3re = model.args
                Δs = [x1re, x2re + im * x2im, x3re + im * x3im]
            end
        else
            error("Unrecognized lattice")
        end
        return Complex.(Δs)
    end

    """
    For different pairing patterns, 
    pick out the independent components `xs` from `Δs`
    """
    function get_xs(Δs::Vector{<:Complex}, lattice::Symbol, pairing::Symbol)
        @assert pairing in spinless_fermion_pairings[lattice]
        # fix overall phase of Δs
        tmp = imag(Δs[1])
        @assert isapprox(tmp, 0.0; atol = 1.0e-16) "Got nonzero imag(Δs[1]) = $(tmp)."
        if lattice == :square
            if pairing == :none
                xs = [real(Δs); imag(Δs)[2:end]]
            elseif pairing == :py
                @assert imag(Δs[2]) == 0.0
                xs = [real(Δs[2])]
            else
                xs = [real(Δs[1])]
            end
        elseif lattice == :triangle
            if pairing == :none
                xs = [real(Δs); imag(Δs)[2:end]]
            else
                xs = [real(Δs[1])]
            end
        else
            error("Unrecognized lattice")
        end
        return xs
    end

    function collect_xs(model::BCSSpinlessFermion)
        args = model.args
        if model.lattice == :square && model.pairing == :none
            return [args[:x1re], args[:x2re], args[:x2im]]
        elseif model.lattice == :triangle && model.pairing == :none
            return [args[:x1re], args[:x2re], args[:x3re], args[:x2im], args[:x3im]]
        else
            return [args[:x]]
        end
    end

    """Pair potential
    ```
    Δ_k = 2 i ∑_η Δ_η sin(kη)
    ```
    where `Δ_η = g_η <c_i c_{i+η}>`
    """
    function Delta(k::Vector{Float64}, model::BCSSpinlessFermion)
        Δs = get_Δs(model)
        return 2.0im * sum(Δ * sin(2 * π * k' * η) for (Δ, η) in zip(Δs, eachcol(model.nbs)))
    end

    """
    Band function
    ```
    ξ_k = -μ + ∑_η (t_η e^{ikη} + h.c.)
    ```
    When all t_η are real, the expression reduces to
    ```
    ξ_k = -μ + ∑_η 2 t_η cos(kη)
    ```
    """
    function xi(k::Vector{Float64}, model::BCSSpinlessFermion)
        μ, t = model.args[:μ], model.args[:t]
        # only consider isotropic nearest neighbor hopping
        return -μ + 2 * t * sum(cos(2 * π * k' * e) for e in eachcol(model.nbs))
    end

    """
    Quasi-particle dispersion
    ```
    E_k = sqrt(ξ_k^2 + |Δ_k|^2)
    ```
    """
    function Energy(k::Vector{Float64}, model::BCSSpinlessFermion)
        return sqrt(xi(k, model)^2 + abs2(Delta(k, model)))
    end

    """
    Calculate density of fermion 
    ```
    δ = <c†_i c_i>
    ```
    """
    function dens_k(k::Vector{Float64}, model::BCSSpinlessFermion)
        β = model.args[:β]
        ξ_k = xi(k, model)
        E_k = Energy(k, model)
        n_k = nf(β * E_k)
        return real(-1 / 2 * ξ_k / E_k * (1 - 2n_k))
    end
    function cal_dens(model)
        return 0.5 + sum_bz(k -> dens_k(k, model), model.bz)
    end

    """
    Self-consistecy equations of Δ_η
    ```
    Δ_η = g_η <c_i c_{i+η}>
    ```
    """
    function cal_Δk(k::Vector{Float64}, i::Int, model::BCSSpinlessFermion)
        g, β = model.args[:g], model.args[:β]
        Delta_k = Delta(k, model)
        E_k = Energy(k, model)
        n_k = nf(β * E_k)
        η = model.nbs[:, i]
        return (-1.0im * g / 2) * Delta_k / E_k * (1 - 2n_k) * sin(2 * π * k' * η)
    end

    function cal_Δs(model::BCSSpinlessFermion)
        return collect(
            sum_bz(k -> cal_Δk(k, i, model), model.bz) for i in 1:size(model.nbs, 2)
        )
    end

    function mf_equations(model::BCSSpinlessFermion)
        xs1 = collect_xs(model)
        Δs2 = cal_Δs(model)
        xs2 = get_xs(Δs2, model.lattice, model.pairing)
        diff = if model.force_pairing
            1 .- (xs2 ./ xs1)
        else
            xs2 - xs1
        end
        return norm(diff)^2, diff
    end

    function freeE_k(k::Vector{Float64}, model::BCSSpinlessFermion)
        β = model.args[:β]
        ξ_k = xi(k, model)
        E_k = Energy(k, model)
        fe_k = 1 / 2 * (ξ_k - E_k)
        if β < Inf && β > 0
            fe_k += (-1 / β) * log1pexp(-β * E_k)
        end
        return fe_k
    end

    function free_energy(model::BCSSpinlessFermion)
        μ, g = model.args[:μ], model.args[:g]
        Δs = get_Δs(model)
        δ = cal_dens(model)
        fe = norm(Δs)^2 / g + μ * δ + sum_bz(k -> freeE_k(k, model), model.bz)
        return fe
    end

end

mf_equations(model::BCSSpinlessFermion) = _BCSSpinless.mf_equations(model)
free_energy(model::BCSSpinlessFermion) = _BCSSpinless.free_energy(model)
