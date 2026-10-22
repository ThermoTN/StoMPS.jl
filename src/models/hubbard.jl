"""Construct a snake-ordered Hubbard lattice with independent OBC/PBC axes.

`lattice_kind=:square` selects square geometry. `:triangular` uses
primitive vectors (sqrt(3)/2, 1/2), (0, 1), integer sites (x,y), and periods
(L,0), (0,W) in those primitive coordinates. The OBC/PBC cylinder is YC.
The triangular branch uses U1SU2 and returns `sublattice=nothing`.

The SU2SU2 Hamiltonian includes the chemical-potential term at mu=U/2.
U1SU2 returns the physical H and a separate number MPO. Each undirected bond
is counted once, also for periodic axes of length two.
"""
function hubbard_model(c; lattice_kind=:square)
    c.L > 1 && c.W >= 1 || throw(ArgumentError("require L>1 and W>=1"))
    c.boundary_x in (:OBC, :PBC) && c.boundary_y in (:OBC, :PBC) ||
        throw(ArgumentError("boundaries must be :OBC or :PBC"))
    c.W == 1 && c.boundary_y != :OBC &&
        throw(ArgumentError("a chain requires transverse OBC"))
    lattice_kind in (:square, :triangular) || throw(ArgumentError("unknown lattice kind"))
    lattice_kind == :triangular && c.symmetry != :u1su2 &&
        throw(ArgumentError("triangular Hubbard geometry supports U1SU2 only"))
    bc = if c.boundary_x == c.boundary_y == :OBC
        OpenBoundaryCondition()
    elseif c.boundary_x == c.boundary_y == :PBC
        CompositeBoundaryCondition(PeriodicBoundaryCondition((c.L, 0)),
            PeriodicBoundaryCondition((0, c.W)))
    elseif c.boundary_x == :PBC
        PeriodicBoundaryCondition((c.L, 0))
    else
        PeriodicBoundaryCondition((0, c.W))
    end
    sites = [(x,y) for x in 1:c.L for y in 1:c.W]
    lattice = if lattice_kind == :square
        SquareLattice(((1.0, 0.0), (0.0, 1.0)), sites, bc)
    else
        TriangularLattice(((sqrt(3)/2, 0.5), (0.0, 1.0)), sites, bc)
    end |> Snake!
    bonds = collect(neighbor(lattice; d=1.0))
    positions = [coordinate(lattice, i) for i in 1:c.L*c.W]
    sublattice = lattice_kind == :square ?
        [iseven(round(Int, x+y)) ? :A : :B for (x,y) in positions] : nothing
    if c.symmetry == :su2su2
        c.mu == c.U/2 || throw(ArgumentError("SU2SU2 requires mu=U/2"))
        all(sublattice[i] != sublattice[j] for (i,j) in bonds) ||
            throw(ArgumentError("SU2SU2 requires a bipartite lattice; periodic axes must be even"))
        H = SU2SU2Fermion.half_filled_hubbard_hamiltonian(
            (; lattice, bonds, sublattice), (; hopping=c.t, interaction=c.U))
        N = nothing
        physical_space = su2su2_born_space()
    elseif c.symmetry == :u1su2
        tree = InteractionTree(c.L*c.W)
        for (i,j) in bonds, (a,b) in ((i,j),(j,i))
            addIntr!(tree, U1SU2Fermion.FdagF, (a,b), (true,true), -c.t;
                Z=U1SU2Fermion.Z, name=(:Fdag,:F))
        end
        for i in 1:c.L*c.W
            addIntr!(tree, U1SU2Fermion.nd, i, c.U; name=:nd)
        end
        H = AutomataMPO(tree)
        N = hubbard_number_mpo(c.L*c.W)
        physical_space = u1su2_born_space()
    else
        throw(ArgumentError("unsupported Hubbard symmetry"))
    end
    return (; lattice, bonds, positions, sublattice, H, N, physical_space)
end

"""U1SU2 MPO for sum(S_i ⋅ S_j) over distinct, undirected `bonds`.

Uses physical spin-1/2 operators, with singlet/triplet eigenvalues -3/4, 1/4.
The returned MPO represents the extensive spin sum with these eigenvalues.
"""
function hubbard_spin_correlation_mpo(nsites::Integer, bonds)
    pairs = [(min(i,j), max(i,j)) for (i,j) in bonds]
    !isempty(pairs) && all(ij -> 1 <= ij[1] < ij[2] <= nsites, pairs) ||
        throw(ArgumentError("spin correlations require nonempty, distinct-site bonds"))
    allunique(pairs) || throw(ArgumentError("each undirected spin bond must occur once"))
    tree = InteractionTree(nsites)
    for pair in pairs
        addIntr!(tree, U1SU2Fermion.SS, pair, (false,false), 1.0; name=(:S,:S))
    end
    return AutomataMPO(tree)
end

function hubbard_number_mpo(nsites)
    tree = InteractionTree(nsites)
    for i in 1:nsites
        addIntr!(tree, U1SU2Fermion.n, i, 1.0; name=:n)
    end
    return AutomataMPO(tree)
end
