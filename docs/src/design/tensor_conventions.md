# Tensor and symmetry conventions

These conventions apply to package kernels and to any example that constructs
FiniteMPS tensors. Tensor dimensions alone are not sufficient validation:
compare spaces, orientations, and values.

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

Use TensorKit/FiniteMPS operations such as `TensorMap`, `AbstractTensorMap`,
`axpy!`, `rmul!`, `norm`, `blocks`, isometries, and `space` rather than reaching
into dense storage.

## Symmetric spaces and fusion

Use the FiniteMPS definitions of symmetric spaces and preserve exact quantum
numbers such as `1//2`. The neutral boundary is `trivial(d)`. Fusion argument
order is semantic, especially in product and braided categories. Changes to
factor order, duality or association require the corresponding categorical maps.

For the validated Z₂ stack, physical, odd auxiliary, selected-sector, fused,
and trivial spaces use dictionary-backed `SortedVectorDict` storage. Preserve
the storage representation of an input space.

`sectors(d)` iteration order is implementation order, not a physical label.
Persist indices such as `r_list` together with the exact sector ordering and
translate them through the same local space when interpreting saved data.
Track multiplicity and carrier degeneracy separately. ``SU(2)\times SU(2)`` carrier dimension
two is not multiplicity two.

The fixed-N U(1)×SU(2) backend constructs its public sector list explicitly as
`q=-1,0,+1`, defining the occupancy indices in that order. The
corresponding marker carrier dimensions are `[1,2,1]`. Its remaining-charge
dynamic program uses public `LocalLeftTensor`/`LocalRightTensor` wrappers and
explicit TensorKit contractions.

## ``SU(2)\times SU(2)`` hopping gauge

The ``SU(2)\times SU(2)`` fermion tensor carries auxiliary irrep `(1/2,1/2)`. TensorKit's
uncoupled eta basis on B differs from A by a pi pseudospin rotation. An
ordered A-to-B bond therefore uses `hopping_AB`; B-to-A uses `hopping_BA`,
which flips the charge-doublet reduced blocks at both ends. Both orientations
require the parity tensor and the validated `sqrt(2)` endpoint factors.

Any gauge or normalization change requires exact two-site A-to-B and B-to-A
Fock-matrix regressions.

## TDVP interface

FiniteMPS owns TDVP evolution. StoMPS workflows call the public `TDVPSweep1!`
interface with explicit truncation and controlled bond expansion. Dependency
and interface changes require explicit version compatibility and focused
numerical validation.
