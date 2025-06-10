function _get_klim(Ns::Vector{Int})
    M0s = collect((N % 2 == 0) ? -div(N, 2) : -div(N - 1, 2) for N in Ns)
    M1s = collect((N % 2 == 0) ? div(N, 2) - 1 : div(Ns - 1, 2) for N in Ns)
    return M0s, M1s
end

"""
Calculate all discrete values of momenta in the 1st Brillouin zone.
`Ns` is the number of lattice sites along each lattice basis vector.
"""
function get_BZ_ks(Ns::Vector{Int})
    M0s, M1s = _get_klim(Ns)
    ks = collect(Ms ./ Ns for Ms in Iterators.product(map(range, M0s, M1s)...))
    return ks
end
