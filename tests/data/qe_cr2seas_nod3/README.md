# Real pw.x output, Cr2Se2As, same run with the vdW correction removed

Identical to `qe_cr2seas_dense/` except `vdw_corr = 'grimme-d3'` is absent.
Every eigenvalue is bit-identical, which is the evidence that Grimme's
dispersion correction -- a classical term depending only on nuclear
positions -- does not enter the Kohn-Sham Hamiltonian and therefore cannot
affect the harmonics at fixed geometry.

Committed because it is what lets the defumat cross-check (P7) proceed
without D3, which defumat does not implement, rather than having to rerun
everything under D2 to match.
