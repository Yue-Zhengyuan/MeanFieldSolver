"""
Momenta in the D-dimensional 1st Brillouin zone (under reciprocal lattice basis)
```
    kᵢ = (mᵢ + δ)/Nᵢ ∈ (−1/2, 1/2] 
```
with the shift δ=0 if `pbc==true` (periodic), and δ=1/2 if `pbc==false` (anti-periodic).
"""
struct BrillouinZone{D}
    # number of sites in each dimension
    Ns::NTuple{D, Int}
    # boundary condition in each dimension (true: PBC; false: anti-PBC)
    pbcs::NTuple{D, Bool}
    # all momenta in 1st Brillouin zone
    ks::Array{Vector{Float64}, D}
end

"""
Create the 1st Brillouin zone, where `Ns` and `pbcs` specify 
the number of sites and boundary condition along each dimension.
"""
function BrillouinZone(Ns::NTuple{D, Int}, pbcs::NTuple{D, Bool}) where {D}
    ranges = Vector{UnitRange{Int}}(undef, D)
    for (i, (N, pbc)) in enumerate(zip(Ns, pbcs))
        if pbc
            if iseven(N)
                mn = -(N ÷ 2) + 1
                mx = N ÷ 2
            else
                mn = -((N - 1) ÷ 2)
                mx = (N - 1) ÷ 2
            end
        else # anti-pbc
            if iseven(N)
                mn = -(N ÷ 2)
                mx = (N ÷ 2) - 1
            else
                # for odd N the range is same as periodic
                mn = -((N - 1) ÷ 2)
                mx = (N - 1) ÷ 2
            end
        end
        ranges[i] = mn:mx
    end
    # take the Cartesian product of all ranges
    shifts = collect(pbc ? 0.0 : 0.5 for pbc in pbcs)
    ks = map(Iterators.product(ranges...)) do ms
        collect(ms .+ shifts) ./ Ns
    end
    return BrillouinZone{D}(Ns, pbcs, ks)
end
function BrillouinZone(Ns::NTuple{D, Int}, pbcs::Bool) where {D}
    return BrillouinZone(Ns, ntuple(_ -> pbcs, D))
end

Base.size(bz::BrillouinZone) = bz.Ns
