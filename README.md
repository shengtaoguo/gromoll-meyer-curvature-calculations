# Curvature Calculations for the Gromoll–Meyer Sphere

This repository contains the exact algebraic verification accompanying
*A Positively Curved Metric on the Gromoll–Meyer Sphere*. The scripts
reconstruct the curvature model from the metric formulas and verify the
polynomial lower bounds over the complete parameter domain.

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

## Files

| File | Purpose |
| --- | --- |
| [verify_gm.py](verification/verify_gm.py) | Reconstructs the polynomial witness and verifies all exact sign bounds. |
| [gm_geometry.py](verification/gm_geometry.py) | Rebuilds the source geometry, dual equations, and full symbolic parameter identities. |

## Verification

The geometry calculation includes differentiated horizontality, the complete
point-and-plane first variation, the feasible dual vectors, and the identities
in the row-mass parameter. The verifier then reconstructs five angular
functions, checks their positive denominators and both antipodal symmetries,
and verifies all 38 Bernstein bounds on a cover of 15 rectangles.

The resulting A-family Taylor coefficient is greater than `37/560`, and its
reduced second derivative is greater than `1/8`. The B-family check verifies
the rational constants in the paper's direct Cauchy–Schwarz estimate, giving
a reduced second derivative greater than `1/48`.

The identity and sign checks use exact rational and integer arithmetic. No
curvature samples or saved polynomial data are inputs. The quotient geometry,
complete zero-plane classification, smoothness, and small-parameter completion
are proved in the paper.
