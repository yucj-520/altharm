# Real pw.x output, Cr2Se2As PBE+U, dense 152-point circle

Same run as `tests/data/qe_cr2seas_p3/` (same scf, same converged density,
same structure, same k_z = 0.25 b3), but a `bands` step re-run on a
152-point circle (l_max = 36 sizing, `4*l_max+8` points) instead of 56,
reusing the already-converged charge density directly (no scf rerun).

Produced 2026-09-20 specifically to answer P3's open question from the
56-point run: whether the ~0.4-1.2 meV floor seen there at the
odd-non-multiple-of-3 forbidden harmonics was genuine DFT noise or an
artefact of not having enough k-points to resolve the angular structure
this band actually carries. Answer: mostly the latter -- see PLAN.md P3.
