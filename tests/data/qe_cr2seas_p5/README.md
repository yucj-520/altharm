# Real pw.x output, Cr2SeAs PBE+U at two further cutoffs

Same structure, same 152-point circle, same k_z as `qe_cr2seas_dense/`, with
scf re-converged in each basis (a density converged at 60 Ry is not a clean
starting point for an 80 Ry band structure).

- `ecut66_*` — 66 Ry / ecutrho 264. Predicted by `count_planewaves` to have
  the *smallest* basis-set discontinuity anywhere in 50-100 Ry (14.5 ppm).
- `ecut80_*` — 80 Ry / ecutrho 320. Predicted dirty (975 ppm), and a higher
  cutoff, so it tests "does raising ecutwfc help?" as well.

Both converged in 35 iterations to E_F = 8.1169 / 8.1166 eV with zero net
magnetization. Committed because together with the 60 Ry run they are the
evidence that falsified the basis-set-discontinuity hypothesis -- see
PLAN.md P4/P5.
