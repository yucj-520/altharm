---
layout: default
title: "Pitfalls"
nav_order: 6
---

# Pitfalls

Each of these produced a wrong but plausible result at some point in this project.
For the first four, `altharm` does the work itself or shows you what to check.

| Pitfall | What goes wrong | What to do |
|---|---|---|
| **Symmetry tolerance** | With the default tolerance of spglib, the mirror that makes Cr₂Se₂As an altermagnet is lost and every l is reported as allowed. | `altharm` scans the tolerance and refuses the unsafe case. Use `--scan` to see the scan. |
| **Band index** | Eigenvalues are sorted by energy. "Band 44" is a different state in a model with a different number of bands, and the n-th band changes character at a crossing. | Check the reported E − E_F. Follow a band by character (`altharm.bandtrack`). |
| **Plane of the circle** | A circle in the y–z plane analysed in the x–y frame gives a plausible but wrong set of allowed l. | The plane is read from the k-points. |
| **Angles from fractional coordinates** | In a hexagonal cell this turns the circle into an ellipse and mixes harmonics. | `altharm` uses Cartesian k only. |
| **Self-consistency on the circle** | Re-converging the density on the band k-points gave 9 meV of noise against a 1 meV signal. | Use a non-self-consistent run (`calculation = 'bands'`). |
| **Hubbard card not inherited** | A `bands` run without the `HUBBARD` card silently gives plain PBE bands. | Repeat the card; check with `altharm.qe.reader.hubbard_is_active`. |
| **Pseudopotentials for DFT+U** | SG15 ONCV files have no atomic wavefunctions, so `pw.x` stops. | Use families with atomic wavefunctions (GBRV, pslibrary PAW). |
| **Valence configuration** | The default pslibrary Se PAW file includes the 3d states, which changes the electron count and the meaning of the band index. | Check the number of electrons when changing pseudopotentials. |
| **Number of points on the circle** | A multiple of the period of the allowed set hides the error. | Choose it coprime; see [How far to trust the numbers](reliability.md). |
| **`celldm(1)` with `CELL_PARAMETERS angstrom`** | `pw.x` stops with "lattice parameter specified twice". | The input writer of `altharm` uses `alat` units. |

More detail on each point is in
[NOTES.md](https://github.com/yucj-520/altharm/blob/master/NOTES.md).
