# Curvature Calculations for the Gromoll–Meyer Sphere

Computations accompanying *A Positively Curved Metric on the Gromoll–Meyer
Sphere*. The scripts verify the source identities and polynomial inequalities
used in the paper.

## Run

Use Python 3.11 and run from the repository root:

```sh
python3 -m pip install -r verification/requirements.txt
python3 verification/verify_gm.py
```

Success prints `PASS` and the curvature bounds; failures exit with a nonzero
status. The default run creates no output files. Assertions must be enabled:
do not use Python's `-O` option.

The two-script verification was tested with Python 3.11.9 on Windows and
took about three minutes on one CPU. The dependency versions are pinned in
[requirements.txt](verification/requirements.txt).

## Paper correspondence

The [main](verification/verify_gm.py#L538) routine calls
[build_model](verification/gm_geometry.py#L779) for the source identities,
then checks the polynomial inequalities. References below use the paper's
equation and proposition numbers; Section 8 describes the calculation.

| Paper location | Code entry | What to compare |
| --- | --- | --- |
| Equations (2.6)–(2.7), frame conversion (5.2) | [_H_tilde](verification/gm_geometry.py#L267), [_J_jets](verification/gm_geometry.py#L281) | The tensors `H/sqrt(2)` and `J`, expressed in the rational Lie frame. |
| Equations (5.1), (7.2)–(7.5) | [_ASource](verification/gm_geometry.py#L318) | The product metric, A-family parameters, horizontal commuting pair, projection, and Gram determinant `15`. |
| Equations (5.4), (5.7)–(5.10) | [_ASource](verification/gm_geometry.py#L318), [_principal_gradient](verification/gm_geometry.py#L488), [_raw_constraints](verification/gm_geometry.py#L586) | The metric and O'Neill terms, differentiation in all 54 point-and-vector variables, and the affine horizontal constraints. |
| Equations (7.6)–(7.7), Appendix B | [_explicit_source_dual](verification/gm_geometry.py#L604), [_twenty_variable_model](verification/gm_geometry.py#L669), [_assemble](verification/gm_geometry.py#L732) | The feasible dual, five homogeneous columns, projection identities, coordinate change, and quadratic model. |
| Proposition 7.1, equation (7.8), Appendix C | [build_model](verification/gm_geometry.py#L779), [_ASource.J_variation](verification/gm_geometry.py#L412) | The full symbolic mass identities, both row poles, and the intrinsic first variation of `J`. |
| Equations (7.10)–(7.12) | [reconstruct_five](verification/verify_gm.py#L320) | The constant inverse, four-dimensional matrix polynomial, full quadratic gain, and five degree-four Bernstein coefficients. |
| Equation (7.13), final bound in Proposition 7.2 | [verify_irrational_model](verification/verify_gm.py#L296), [main](verification/verify_gm.py#L538) | The irrational-part identity, its square bound, and the rational comparison giving the A Taylor bound `37/560`. |
| Proposition 7.2, Section 8.1 | [primitive_numerator](verification/verify_gm.py#L501), [direct_box_coefficients](verification/verify_gm.py#L397), [verify_cover](verification/verify_gm.py#L411) | Positive denominators, both antipodal identities, and all 38 exact Bernstein bounds over the complete 15-rectangle cover. |
| Rational constants after (6.10), final bound in Proposition 6.1 | [check_b_bound](verification/verify_gm.py#L523) | The two rational norm sums and the final B bound `1/48`; the geometric source identity is proved by hand in Section 6. |

## Reading the variables

In [_assemble](verification/gm_geometry.py#L732), `matrix`, `drive`, `baseR`,
and `baseI` are the matrix, linear term, and two base terms in equation (7.7).
The output of [build_model](verification/gm_geometry.py#L779) contains their
mass coefficients: `A0`, `A1`, `d0`, `d1`, and the three coefficients of each
base term in equation (7.8). The base terms retain the paper's `15/2`
normalization.

In [reconstruct_five](verification/verify_gm.py#L320), `small`,
`reduced_drive`, and `retained` are the paper's `S`, `b_*`, and `v_*`.
The code's `x` denotes `x_0`. The factor `s` enters through `s^2 = 1-m`
when forming the gain in equation (7.11).
The returned `elevated` list is exactly the five functions in equation (7.12).

## Scope

The checks use exact rational and integer arithmetic, with no sampled curvature
values or saved polynomial data. The [subdivision](verification/verify_gm.py#L47)
is checked on both children of every split and for all five functions on each
terminal rectangle.

The B routine evaluates the norm sums in the Cauchy–Schwarz estimate. Its
geometric source identity is proved in Section 6. The quotient geometry,
zero-plane classification, smoothness, and small-parameter argument are in
Sections 2–4.
