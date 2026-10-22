"""
Local SU(2)charge × SU(2)spin tensors for the half-filled Hubbard model.

The two factors in every sector are ordered as `(charge SU(2), spin SU(2))`.
The local Hilbert space decomposes as

    (C=1/2, S=0) + (C=0, S=1/2),

where the first multiplet contains the empty and doubly occupied states and
the second contains the singly occupied states, as in Sec. II B of
[arXiv:2609.35754v1](https://arxiv.org/html/2609.35754v1).
"""
module SU2SU2Fermion

using ..StoMPS: su2su2_born_space
using FiniteMPS
using FiniteMPS.TensorKit

export half_filled_hubbard_hamiltonian

const pspace = su2su2_born_space()
const charge_sector = Irrep[SU₂ × SU₂](1 // 2, 0)
const spin_sector = Irrep[SU₂ × SU₂](0, 1 // 2)

"""Fermion parity `(-1)^n`: even on the charge doublet and odd on the spin doublet."""
const Z = let
    parity = TensorMap(ones, Float64, pspace, pspace)
    block(parity, spin_sector) .= -1
    parity
end

"""
The half-filled on-site operator

    n_up n_down - (n_up + n_down) / 2.

Multiplying this tensor by `U` implements both the interaction and
`mu = U/2` in the charge-SU(2)-invariant on-site term.
"""
const half_filled_onsite = let
    onsite = TensorMap(zeros, Float64, pspace, pspace)
    block(onsite, spin_sector) .= -1 / 2
    onsite
end

"""Build the SU(2)charge × SU(2)spin half-filled Hubbard MPO on a lattice datum.

`lattice_data` must provide `lattice`, `bonds`, and `sublattice`; `config` must
provide `hopping` and `interaction`. The tensor pair already contains hopping
plus its Hermitian conjugate, so each unique undirected bond is inserted once.
"""
function half_filled_hubbard_hamiltonian(lattice_data, config)
    tree = InteractionTree(size(lattice_data.lattice))
    for bond in lattice_data.bonds
        first_site = first(bond)
        hopping_tensor = lattice_data.sublattice[first_site] == :A ?
            hopping_AB : hopping_BA
        hopping_name = lattice_data.sublattice[first_site] == :A ?
            (:SU2SU2F_AB_left, :SU2SU2F_AB_right) : (:SU2SU2F_BA_left, :SU2SU2F_BA_right)
        addIntr!(
            tree,
            hopping_tensor,
            bond,
            (true, true),
            -config.hopping;
            Z=Z,
            name=hopping_name,
        )
    end
    for site in 1:size(lattice_data.lattice)
        addIntr!(
            tree,
            half_filled_onsite,
            site,
            config.interaction;
            name=:half_filled_onsite,
        )
    end
    return AutomataMPO(tree)
end

"""
Two fermionic rank-3 tensors whose contracted bond is

    sum_sigma (c_i_sigma^dagger c_j_sigma + h.c.)

for opposite sublattices. The auxiliary irrep is `(C=1/2, S=1/2)`.
The `sqrt(2)` factors make the interaction-tree strength equal to the physical
hopping coefficient. FiniteMPS supplies the Jordan-Wigner parity string from
`Z` for the fermionic contraction.
"""
function _hopping_tensor_pair(charge_gauge::Int)
    charge_gauge in (-1, 1) || throw(ArgumentError("charge gauge must be +1 or -1"))
    auxiliary = Rep[SU₂ × SU₂]((1 // 2, 1 // 2) => 1)
    left = sqrt(2) * TensorMap(ones, Float64, pspace, pspace ⊗ auxiliary)
    right = sqrt(2) * TensorMap(ones, Float64, auxiliary ⊗ pspace, pspace)

    block(left, charge_sector) .*= charge_gauge
    block(right, charge_sector) .*= charge_gauge
    # This relative reduced-matrix-element sign makes the scalar contraction
    # Hermitian after FiniteMPS applies the right-end parity convention.
    block(right, spin_sector) .*= -1
    return left, right
end

"""Hopping tensor for an ordered MPS bond whose first site is on sublattice A."""
const hopping_AB = _hopping_tensor_pair(1)

"""Hopping tensor for an ordered MPS bond whose first site is on sublattice B."""
const hopping_BA = _hopping_tensor_pair(-1)

end # module SU2SU2Fermion
