---
layout: default
title: "Results"
nav_order: 7
---

# Results

The calculations behind this page are stored in the repository (`tests/data/`) and are
checked by the test suite. You can reproduce every table with one command.

## Cr₂Se₂As

PBE+U (U = 3.0 eV), band nearest the Fermi level, 152 points on a circle at
k_z = 0.25 b₃.

```bash
altharm tests/data/qe_cr2seas_sm24/data-file-schema.xml --lmax 46
```

| l | 3 | 9 | 15 | 21 | 27 | 33 | 39 | 45 |
|---|---|---|---|---|---|---|---|---|
| Amplitude (meV) | 64.293 | 43.447 | 3.603 | 8.662 | 1.876 | 2.583 | 1.182 | 0.985 |

- Symmetry allows exactly the odd multiples of 3, and nothing in the plane k_z = 0.
- The splitting is carried by two harmonics of comparable size, l = 3 and l = 9.
- Noise floor: about 0.02 meV, more than a thousand times below the two leading
  harmonics.
- At identical DFT settings, the result agrees within 1 % with an independent Wannier
  model for all eight allowed harmonics.
- The values depend on U; see [How far to trust the numbers](reliability.md).

## MnTe

```bash
altharm tests/data/mnte_kz025_u4/data-file-schema.xml
```

Allowed l = 3, 9, 15, …; the leading harmonic is l = 3 (g-wave, see
[Concepts](concepts.md)), and the next one, l = 9, is ten times smaller.

## CoF₂

```bash
altharm tests/data/cof2_u4_kz025/data-file-schema.xml
```

Allowed l = 2, 6, 10, …; the leading harmonic is l = 2 (d-wave) with 6.9 meV, followed by
l = 6 with 1.9 meV.

## Two implementations

`qe_native/harmwave.x` is a second implementation compiled into Quantum ESPRESSO. It uses
the symmetry finder and the linear algebra of Quantum ESPRESSO instead of spglib and
NumPy, and agrees with the Python route within 10⁻⁹ eV on the same data.
