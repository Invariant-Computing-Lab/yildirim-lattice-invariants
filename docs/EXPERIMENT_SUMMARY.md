# Final Experiment Summary

Reference run: 12 August 2026  
Master seed: `20260812`  
Random lattice-zonotope atlas: `N=2000`  

## Exact results

1. For every nonempty convex polytope `P`, the stated recursion satisfies
   `Conv(P x {-1,0,+1}) = P x [-1,1]`.
2. Starting from a point, the recursive family is therefore the standard cube family `[-1,1]^n`.
3. Exact integer checks for `n=1,...,50` verified the cubical face formula, recurrence, ternary face-count identity, and boundary Euler identity.
4. Exact adversarial examples show that Ehrhart data alone do not determine structural oriented-matroid class, and oriented-matroid class alone does not determine lattice-polarized identity.
5. The supplied rational one-parameter family crosses the targeted discriminant wall exactly at `t=0`, where the relevant maximal minor is `12t`.

## Finite-atlas results

For the declared operational class label on the 2000-object rank-3, six-generator lattice-zonotope atlas:

| Signature | Unique signatures | Reconstructibility R | Collision rate |
|---|---:|---:|---:|
| rank | 1 | 0.00000 | 1.000000 |
| f-vector | 1 | 0.00000 | 1.000000 |
| oriented matroid | 4 | 0.16114 | 0.321363 |
| volume | 764 | 0.84382 | 0.001480 |
| Ehrhart | 1765 | 0.97754 | 0.000133 |
| OM + Ehrhart | 1919 | 0.99248 | 0.000043 |
| arithmetic multiplicity profile | 2000 | 1.00000 | 0.000000 |
| OM + arithmetic profile | 2000 | 1.00000 | 0.000000 |

`R=1` here is finite-atlas operational sufficiency, not a universal completeness theorem.

## Quotient-preservation controls

- Structural `GL(r,R)` + signed-permutation quotient, `m=6`: 1000/1000 canonical preservation passes.
- Held-out structural `m=7`: 100/100 passes.
- Lattice `GL(r,Z)` + signed-permutation quotient: 300/300 arithmetic-profile preservation passes.
- Raw, noncanonical signatures were unchanged in 0 of the positive-control transformations, confirming that canonicalization removes genuine representation nuisance rather than testing trivial identity maps.

## Controlled polarized-toric families

For axis-aligned lattice boxes modulo coordinate permutation:

- dimension 3, coefficients 1-10: 220 classes; one unpolarized toric type; 120 volume signatures; 220 Ehrhart signatures.
- dimension 4, coefficients 1-7: 210 classes; one unpolarized toric type; 140 volume signatures; 210 Ehrhart signatures.

These are controlled-family discrimination results, not a general classification theorem.
