# Tensor and symmetry conventions

These conventions apply to package kernels and to any example that constructs
FiniteMPS tensors. Validate dimensions, spaces, orientations, and values
together.

## TensorMap and MPSTensor order

`TensorMap(codomain, domain)` maps `domain -> codomain`. A prime such as `V'`
denotes the dual space; preserve this duality when matching tensor legs.

For a rank-3 MPS tensor use

```text
A[left_bond, physical ; right_bond]
```

For the rank-4 MPO/tail form use

```text
A[left_bond, physical_out ; physical_in_or_trivial_tail, right_bond]
```

The first two legs are codomain legs and the final virtual leg is the right
bond. Neighboring tensors must satisfy
`domain(A[i])[end] == codomain(A[i+1])[1]`; boundary virtual spaces are
`trivial(d)` unless an API explicitly requests another total sector.

In `@tensor`, list free codomain indices before `;` and free domain indices
after it. Contracted labels occur exactly twice. In `permute(A,
codomain_indices, domain_indices)`, both tuples carry type-level orientation.
Use TensorKit permutations to preserve this orientation and symmetry data.

TensorKit/FiniteMPS operations include `TensorMap`, `AbstractTensorMap`,
`axpy!`, `rmul!`, `norm`, `blocks`, isometries, and `space`.

## Symmetric spaces and fusion

Use the FiniteMPS definitions of symmetric spaces and preserve exact quantum
numbers such as `1//2`. The neutral boundary is `trivial(d)`. Fusion argument
order is semantic, especially in product and braided categories. Changes to
factor order, duality or association require the corresponding categorical maps.

Born marker spaces use dictionary-backed `SortedVectorDict` storage and must
match the storage representation of the input physical space.

`sectors(d)` follows the space's iteration order.
Persist indices such as `r_list` together with the exact sector ordering and
translate them through the same local space when interpreting saved data.
Track multiplicity and carrier degeneracy separately. Each local
``\mathrm{SU(2)_{charge}}\times\mathrm{SU(2)_{spin}}`` multiplet has carrier
dimension two and multiplicity one, as defined in [Born sampling](born_sampling.md).

The canonical U(1)×SU(2) backend constructs its public sector list explicitly as
`q=-1,0,+1`, defining the occupancy indices in that order. The
corresponding marker carrier dimensions are `[1,2,1]`. Its remaining-charge
dynamic program uses public `LocalLeftTensor`/`LocalRightTensor` wrappers and
explicit TensorKit contractions.
