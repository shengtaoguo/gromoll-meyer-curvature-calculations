# Curvature Calculations for the Gromoll–Meyer Sphere

Companion to *A Positively Curved Metric on the Gromoll–Meyer Sphere*.

## Run

With Python 3.11, run from the repository root:

```sh
python -m pip install -r verification/requirements.txt
python verification/verify_gm.py
```

A successful run prints `PASS` and the A and B curvature bounds.

## Paper correspondence

The [main](verification/verify_gm.py#L554) routine calls [build_model](verification/gm_geometry.py#L782) to reconstruct the source
identities, then verifies the polynomial inequalities. Numbers below refer
to the September 11 manuscript; Section 8 describes the exact computations.

| Paper location | Code entry | What to compare |
| --- | --- | --- |
| Equations (2.6)–(2.7), frame conversion (5.2) | [_H_tilde](verification/gm_geometry.py#L267), [_J_jets](verification/gm_geometry.py#L281) | The tensors `H/sqrt(2)` and `J` in the rational Lie frame. |
| Equations (5.1), (7.2)–(7.5) | [_ASource](verification/gm_geometry.py#L318) | The product metric, A parameters, horizontal commuting pair, projection, and Gram determinant `15`. |
| Equations (5.4)–(5.10), Lemma 7.1 | [_ASource](verification/gm_geometry.py#L318), [_principal_gradient](verification/gm_geometry.py#L488), [_raw_constraints](verification/gm_geometry.py#L586) | The connection and O'Neill terms, first curvature variation, all 54 point-and-vector derivatives, and affine horizontal constraints. |
| Equations (7.6)–(7.7), Appendix B | [_explicit_source_dual](verification/gm_geometry.py#L604), [_twenty_variable_model](verification/gm_geometry.py#L669), [_assemble](verification/gm_geometry.py#L732) | The feasible dual, five homogeneous columns, projection identities, coordinate change, and quadratic model. |
| Proposition 7.2, equation (7.8), Appendix C | [build_model](verification/gm_geometry.py#L782), [_ASource.J_variation](verification/gm_geometry.py#L412) | The full symbolic parameter identities, both row poles, and intrinsic first variation of `J`. |
| Equations (7.10)–(7.12) | [reconstruct_five](verification/verify_gm.py#L323) | The constant inverse, matrix polynomial, complete quadratic gain, and five degree-four Bernstein coefficients. |
| Equation (7.13), final bound in Proposition 7.3 | [verify_irrational_model](verification/verify_gm.py#L299), [main](verification/verify_gm.py#L554) | The irrational-part identity, its square bound, and the rational comparison giving the A Taylor bound `37/560`. |
| Proposition 7.3, Section 8.1 | [primitive_numerator](verification/verify_gm.py#L517), [direct_box_coefficients](verification/verify_gm.py#L400), [verify_cover](verification/verify_gm.py#L414) | Positive denominators, both antipodal identities, the stated bidegrees, and all 38 exact bounds on the complete cover. |
| Rational constants after (6.10), Proposition 6.1 | [check_b_bound](verification/verify_gm.py#L539) | The two norm sums and the final B bound `1/48`; the geometric identity is proved in Section 6. |

## Reading the variables

In [_assemble](verification/gm_geometry.py#L732), `matrix`, `drive`, `baseR`, and `baseI` are the matrix,
linear term, and two base terms in (7.7). In the code, `mass` means the row
parameter `m = |a|^2`. The output of [build_model](verification/gm_geometry.py#L782) contains `A0`, `A1`,
`d0`, `d1`, and the three coefficients of each base term in (7.8).
The base terms retain the paper's `15/2` normalization and enter the feasible lower bound (7.9).

The metadata field `drive` uses the representative `s >= 0`.
`signed_drive` records `s*(d0 + mass*d1)`, with `s^2 = 1-mass`; this signed
identity is checked on the complete symbolic row family.
`A_b_zero_row_pole_drive_exact` refers to the A-family row pole `b = 0`
(`m = 1`).

In [reconstruct_five](verification/verify_gm.py#L323), `small`, `reduced_drive`, and `retained` are the
paper's `S`, `b_*`, and `v_*`. The code's `x` denotes `x_0`; the factor `s`
enters through `s^2 = 1-m` in the gain (7.11). The returned `elevated` list
is the five functions in (7.12).

## Optional result record

```sh
python verification/verify_gm.py --output result.json
```

The result records every comparison's rectangle, function index, bidegree,
exact minimum Bernstein coefficient and its position. Function numbers use
the order in (7.12). These coefficients
belong to the primitive integer numerator after mapping its rectangle to the
unit square. They are not numerical estimates of sectional curvature.

## Scope

The checks use exact rational and integer arithmetic, with no sampled
curvature values or saved polynomial inputs. Both children of every split
are checked, and every function is covered at each terminal rectangle.

Sections 2–5 of the paper establish the quotient geometry, complete zero set,
flat preservation, small-parameter argument, and curvature-variation formula.
The B routine checks the norm sums and final rational constant; its full
geometric argument, including the constant multiplier, is in Section 6.
