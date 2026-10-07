# Real pw.x output, Cr2Se2As, two diagonalisation thresholds

Bands-only reruns on the same 152-point circle, reusing the same converged
60 Ry density (no scf rerun), with `diago_thr_init` loosened from the
production 1e-13 to 1e-10 and 1e-6.

1e-6 Ry corresponds to a nominal eigenvalue error bound of 13.6 ueV, which
is the same scale as the 27 ueV noise floor being investigated -- so if the
per-k-point diagonalisation threshold set that floor, this run would show
it. It does not: see PLAN.md P3/P5.
