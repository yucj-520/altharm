# Real pw.x output, Cr2Se2As, 304-point circle

Bands-only run reusing the same converged 60 Ry density as
`qe_cr2seas_dense/`, on a circle with twice the points (304, Nyquist l=152).

This is the run that confirmed the aliasing prediction: doubling the k-count
drops the noise floor from 27.5 to 1.0 ueV -- a factor of 27 -- while the
allowed harmonics move by 0.2 ueV. See PLAN.md P3.
