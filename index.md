---
layout: default
title: "Home"
nav_order: 1
---

# altharm

**Angular harmonics of the altermagnetic spin splitting, directly from DFT**

In an altermagnet the spin splitting ΔE(k) = E↑(k) − E↓(k) changes sign as k moves around
a circle in reciprocal space. `altharm` measures the shape of that oscillation.

| | |
|---|---|
| **Input** | a finished Quantum ESPRESSO band calculation on a circle of k-points |
| **Output** | the amplitude of every angular harmonic l in ΔE(θ) = Σ [a_l cos(lθ) + b_l sin(lθ)] |
| **Quality control** | symmetry decides which l may be non-zero; the amplitude found at the forbidden l is the noise floor of the calculation |

## Quick start

```bash
pip install git+https://github.com/yucj-520/altharm.git
altharm path/to/prefix.save/data-file-schema.xml
```

The second line prints one row per harmonic: allowed or forbidden, its amplitude in meV,
and how far it stands above the noise floor. See [Quick start](quickstart.md) for the
whole workflow.

## Contents

| Page | Content |
|---|---|
| [Installation](install.md) | `pip install`, requirements |
| [Quick start](quickstart.md) | from a Quantum ESPRESSO run to the table of harmonics |
| [Concepts](concepts.md) | harmonics, allowed and forbidden l, the noise floor, d/g/i-wave labels |
| [How far to trust the numbers](reliability.md) | what the noise floor does and does not measure |
| [Pitfalls](pitfalls.md) | numerical traps found in real calculations |
| [Results](results.md) | Cr₂Se₂As, MnTe and CoF₂ |
| [Reference](reference.md) | command line, Python functions, other routes |

## Source code, questions, citation

- Source code: <https://github.com/yucj-520/altharm>
- Questions and bug reports: <https://github.com/yucj-520/altharm/issues>
- Full discussion of the method: [NOTES.md](https://github.com/yucj-520/altharm/blob/master/NOTES.md)
- If you use `altharm`, please cite the repository
  ([`CITATION.cff`](https://github.com/yucj-520/altharm/blob/master/CITATION.cff)).
