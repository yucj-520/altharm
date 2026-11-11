---
layout: default
title: "How far to trust the numbers"
nav_order: 5
---

# How far to trust the numbers

**In one sentence:** for a given Kohn–Sham Hamiltonian, `altharm` resolves ΔE(θ) to tens
of μeV; it does not tell you whether that Hamiltonian is the right one.

## The floor is the error of the fit, not of the coefficient

Measured on Cr₂Se₂As (median over the eight allowed harmonics):

| Source of error | Change of the coefficients | Seen by the noise floor? |
|---|---|---|
| Hubbard U | about 21 % per eV | no |
| Smearing width chosen too narrow | 12.9 % | no |
| PAW instead of ultrasoft pseudopotentials | 9.7 % | no |
| PBEsol instead of PBE | 7.3 % | no |
| Atomic positions at relaxation precision (10⁻³ Å) | below 4 % | partly |
| Number of points on the circle | 18 % (36 % in the worst case) | no, the floor can even read zero |
| The noise floor itself | 0.04 % | |

Three consequences:

- **Quote a coefficient together with its settings**, in particular the value of U. The
  reference values for Cr₂Se₂As are values at U = 3.0 eV. Linear response gives
  U = 6.39 eV; at that value the coefficients are smaller by a median of 58 %, while the
  ratio of the second to the leading harmonic changes by only 1.5 percentage points.
- **A lower floor does not mean a more accurate result.** Changing the functional or the
  pseudopotential moved the coefficients by 7–10 % and lowered the floor at the same time.
- **A converged density does not mean converged harmonics.** Magnetisation, Hubbard
  occupations and total energy agreed to three digits across a change of smearing that
  moved the coefficients by 27 %.

## Choosing the number of points on the circle

The floor is read from the forbidden harmonics. If the number of points N is a multiple
of the period of the allowed set (6 for l = 3, 9, 15, …), the unresolved high harmonics
fold exactly onto the allowed coefficients: the floor reads zero while the coefficients
are wrong by about 18 %.

**Rule: choose N coprime to that period.** For Cr₂Se₂As, 151, 155 or 157 points are about
five times more accurate than 152 at the same cost.

```python
from altharm.sampling import sampling_report
print(sampling_report(152, [3, 9, 15, 21]).message)
```

```
n = 152 shares the factor 2 with the allowed-set period p = 6: the tail reaches the
allowed coefficients at alias order 3. Usable, but a count coprime to 6 (e.g. 151) is
strictly better at the same cost.
```

## What the floor is good for

It tells you whether `lmax` is high enough, whether the circle is sampled densely enough,
and whether a harmonic is resolved at all. In this project it improved the resolution
from 9 meV to about 1 μeV.

The measurements behind this page are in
[NOTES.md](https://github.com/yucj-520/altharm/blob/master/NOTES.md) and, in full, in
[PLAN.md](https://github.com/yucj-520/altharm/blob/master/PLAN.md) (P3, P12, P13, P14, P19).
