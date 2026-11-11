# altharm

**Angular harmonics of the altermagnetic spin splitting, directly from DFT**

[![tests](https://github.com/yucj-520/altharm/actions/workflows/tests.yml/badge.svg)](https://github.com/yucj-520/altharm/actions/workflows/tests.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://github.com/yucj-520/altharm/blob/master/LICENSE)
![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)

In an altermagnet the spin splitting ΔE(k) = E↑(k) − E↓(k) changes sign as k moves around
a circle in reciprocal space. `altharm` measures the shape of that oscillation:

- **Input:** a finished Quantum ESPRESSO band calculation on a circle of k-points.
- **Output:** the amplitude of every angular harmonic l in
  ΔE(θ) = Σ [a_l cos(lθ) + b_l sin(lθ)].
- **Quality control:** the magnetic symmetry decides which l may be non-zero. Whatever
  amplitude appears at a forbidden l is numerical noise, so every result comes with its
  own noise floor.

Everything is read from the calculation itself: the structure, the plane of the circle,
the band nearest the Fermi level and a safe symmetry tolerance.

## Example

```bash
altharm tests/data/qe_cr2seas_sm24/data-file-schema.xml
```

```
k-points on the circle : 152
band                   : 44  (E - E_F = +0.0373 eV, ...)
allowed l (<= 38)     : [3, 9, 15, 21, 27, 33]
noise floor (RMS over forbidden l) :    0.0215702 meV

   l     status       amplitude     /floor
   1  forbidden       0.0247022         --
   2  forbidden     6.82348e-05         --
   3    allowed         64.2928     2980.6
   ...
   9    allowed         43.4467     2014.2
  15    allowed         3.60328      167.0
  21    allowed         8.66167      401.6
```

This is Cr₂Se₂As: the splitting is carried by l = 3 and l = 9 together (64 and 43 meV),
both thousands of times above the noise floor of 0.02 meV.

## Installation

```bash
pip install git+https://github.com/yucj-520/altharm.git
```

Python 3.10 or newer. NumPy, SciPy, spglib and ASE are installed automatically.

## Quick start

1. **Run Quantum ESPRESSO.** A spin-polarised `scf` calculation, then a
   `calculation = 'bands'` run with `nosym` and `noinv` on a circle of k-points.
   `altharm` can write that input for you
   ([how](https://yucj-520.github.io/altharm/quickstart.html)).
2. **Run altharm** on the XML file of the bands run:

   ```bash
   altharm path/to/prefix.save/data-file-schema.xml
   ```

3. **Read the table.** A harmonic is resolved if it is `allowed` and its amplitude is well
   above the floor (column `/floor`).

## Before you quote a number

- **The noise floor is the error of the fit, not of the coefficient.** The Hubbard U, the
  functional, the pseudopotential and the smearing change the coefficients by 7–21 %
  without changing the floor.
- **Choose the number of k-points on the circle coprime to the period of the allowed
  harmonics** (for example 151 instead of 152 when the allowed l are 3, 9, 15, …).
  `altharm.sampling.sampling_report` tells you.
- **Check the band.** A band index is not portable between calculations; the output
  reports E − E_F of the band it used.

Details: [How far to trust the numbers](https://yucj-520.github.io/altharm/reliability.html).

## Validated on

| Material | Result |
|---|---|
| Cr₂Se₂As (PBE+U, U = 3 eV) | allowed l = 3, 9, 15, …; 64.3 and 43.4 meV at l = 3 and 9. At identical settings the result agrees within 1 % with an independent Wannier model for all eight allowed harmonics. |
| MnTe | g-wave: leading harmonic l = 3 on the circle |
| CoF₂ | d-wave: leading harmonic l = 2 |

A second implementation compiled into Quantum ESPRESSO (`qe_native/harmwave.x`) agrees
with the Python route within 10⁻⁹ eV. The test suite has more than 350 tests on stored
calculations.

## Documentation

**Website: <https://yucj-520.github.io/altharm/>**

| Page | Content |
|---|---|
| [Installation](https://yucj-520.github.io/altharm/install.html) | `pip install`, requirements |
| [Quick start](https://yucj-520.github.io/altharm/quickstart.html) | from a Quantum ESPRESSO run to the table of harmonics |
| [Concepts](https://yucj-520.github.io/altharm/concepts.html) | harmonics, allowed and forbidden l, the noise floor, d/g/i-wave labels |
| [How far to trust the numbers](https://yucj-520.github.io/altharm/reliability.html) | what the noise floor does and does not measure |
| [Pitfalls](https://yucj-520.github.io/altharm/pitfalls.html) | numerical traps found in real calculations |
| [Results](https://yucj-520.github.io/altharm/results.html) | Cr₂Se₂As, MnTe and CoF₂ |
| [Reference](https://yucj-520.github.io/altharm/reference.html) | command line, Python functions, other routes |

## Also in this repository

- [`qe_native/`](https://github.com/yucj-520/altharm/tree/master/qe_native): `harmwave.x`, the same analysis as an executable compiled into Quantum ESPRESSO.
- `altharm.wannier`: the same harmonics from a Wannier90 `hr.dat` model.
- [`NOTES.md`](https://github.com/yucj-520/altharm/blob/master/NOTES.md): the full discussion of the method and its limits.
- [`PLAN.md`](https://github.com/yucj-520/altharm/blob/master/PLAN.md): the complete development record, including results that were later withdrawn.

## Citation

If you use `altharm`, please cite this repository (see
[`CITATION.cff`](https://github.com/yucj-520/altharm/blob/master/CITATION.cff)).

## Licence and author

`altharm` is distributed under the GNU General Public License, version 3 or later.

Author: Cuiju Yu, Aalto University. Questions and bug reports are welcome in the
[issue tracker](https://github.com/yucj-520/altharm/issues).

Parts of the code were developed with the assistance of an AI coding tool (Claude). The
results are validated against stored calculations by the test suite.
