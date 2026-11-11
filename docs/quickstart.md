---
layout: default
title: "Quick start"
nav_order: 3
---

# Quick start

Three steps lead from a structure to the table of harmonics.

## Step 1. Self-consistent calculation

Run a normal spin-polarised `scf` calculation with `pw.x`. Give the two magnetic
sublattices separate species (for example `Cr1` and `Cr2`) with opposite
`starting_magnetization`.

## Step 2. Bands on a circle

Copy the `scf` input and change three things:

- `calculation = 'bands'`,
- `nosym = .true.` and `noinv = .true.` (otherwise the circle is folded into the
  irreducible wedge),
- the `K_POINTS` card, which `altharm` writes for you:

```python
import numpy as np
from altharm.kpath import circle
from altharm.qe.writer import cart_to_tpiba, kpoints_card

alat_bohr = 7.0408                              # celldm(1) of the scf input
b_inplane = 4 * np.pi / (np.sqrt(3) * 3.7258)   # 1/Angstrom, hexagonal cell, a = 3.7258 A
b_z = 2 * np.pi / 16.918                        # 1/Angstrom, c = 16.918 A

points = circle(0.25 * b_inplane, 151, center=(0.0, 0.0, 0.25 * b_z))
print(kpoints_card(cart_to_tpiba(points, alat_bohr)))
```

```
K_POINTS tpiba
151
0.2886775292 0.0000000000 0.0550572010 1.0000000000
0.2884276523 0.0120085499 0.0550572010 1.0000000000
...
```

Two rules for this step:

- **Repeat the `HUBBARD` card.** A `bands` run does not inherit the Hubbard parameters of
  the `scf` run. Without the card, Quantum ESPRESSO silently computes plain PBE bands.
- **Choose the number of points coprime to the period of the allowed harmonics.** Check
  with

  ```python
  from altharm.sampling import sampling_report
  print(sampling_report(151, [3, 9, 15, 21]).message)
  ```

`altharm.qe.writer.bands_input` writes the complete input file if you prefer not to edit
by hand.

## Step 3. Run altharm

```bash
altharm path/to/prefix.save/data-file-schema.xml
```

You can try this at once with a calculation stored in the repository:

```bash
altharm tests/data/qe_cr2seas_sm24/data-file-schema.xml
```

```
source                 : tests/data/qe_cr2seas_sm24/data-file-schema.xml
k-points on the circle : 152
circle radius          : 0.486816 1/A  (spread 1.94e-10, planarity 0.00e+00)
circle axis            : [+0.0000 -0.0000 +1.0000], offset +0.1907 x radius
band                   : 44  (E - E_F = +0.0373 eV, closest approach below 20.1 meV, above 32.1 meV)
symprec                : 5e-05
operations kept        : 24  (12 sublattice-preserving, 12 exchanging)
allowed l (<= 38)     : [3, 9, 15, 21, 27, 33]
|dE| max               : 73.4384 meV

noise floor (RMS over forbidden l) :    0.0215702 meV
worst single forbidden harmonic    :    0.0650803 meV
unfitted residual (RMS)            :      1.34508 meV

   l     status       amplitude     /floor
   1  forbidden       0.0247022         --
   2  forbidden     6.82348e-05         --
   3    allowed         64.2928     2980.6
   ...
   9    allowed         43.4467     2014.2
  15    allowed         3.60328      167.0
  21    allowed         8.66167      401.6
  27    allowed         1.87635       87.0
  33    allowed         2.58318      119.8
```

## How to read the output

| Line | Meaning |
|---|---|
| `circle radius`, `circle axis` | the circle found in the k-points; spread and planarity should be close to zero |
| `band` | the band that was analysed and its energy relative to the Fermi level. Check that this is the band you mean. |
| `symprec`, `operations kept` | the symmetry tolerance that was selected and the number of symmetry operations that keep or exchange the magnetic sublattices |
| `allowed l` | the harmonics that symmetry permits |
| `noise floor` | root mean square of the amplitudes at the forbidden l |
| `unfitted residual` | the part of ΔE(θ) that the fit up to `lmax` does not capture; in this example it is the allowed harmonics above l = 38 |
| table | amplitude of each harmonic in meV; `/floor` is the amplitude divided by the noise floor |

A harmonic is resolved if it is allowed and stands well above the floor.

## Options

| Option | Meaning |
|---|---|
| `--band N` | analyse band N (0-based) instead of the band nearest the Fermi level |
| `--lmax L` | highest harmonic fitted (default: number of k-points divided by 4) |
| `--symprec X` | fix the symmetry tolerance instead of scanning it |
| `--soc` | the calculation included spin–orbit coupling |
| `--scan` | also print the scan of the symmetry tolerance |
