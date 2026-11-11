---
layout: default
title: "Reference"
nav_order: 8
---

# Reference

## Command line

```
altharm <prefix>.save/data-file-schema.xml [--band N] [--lmax L] [--symprec X] [--soc] [--scan]
```

| Argument | Meaning | Default |
|---|---|---|
| `xml` | the `data-file-schema.xml` of a spin-polarised `bands` run on a circle of k-points | |
| `--band N` | 0-based band index | the band nearest the Fermi level |
| `--lmax L` | highest harmonic fitted | number of k-points divided by 4 |
| `--symprec X` | symmetry tolerance for spglib | scanned automatically |
| `--soc` | the run included spin–orbit coupling, so odd l are not excluded by time reversal | off |
| `--scan` | also print the scan of the symmetry tolerance | off |

`python -m altharm ...` is equivalent. The exit status is 2 if no safe symmetry tolerance
is found.

## Python

```python
from altharm import analyze_qe

a = analyze_qe("path/to/prefix.save/data-file-schema.xml")
print(a.summary())            # the report printed by the command line
print(a.amplitudes_mev())     # {l: amplitude in meV} for the allowed l
print(a.noise.floor)          # noise floor in eV
```

| Function | Purpose |
|---|---|
| `altharm.analyze_qe` | the complete analysis of a Quantum ESPRESSO run |
| `altharm.kpath.circle`, `altharm.kpath.sphere_grid` | k-points on a circle or a sphere |
| `altharm.qe.writer.kpoints_card`, `bands_input` | `K_POINTS` card and complete `bands` input for `pw.x` |
| `altharm.qe.reader.read_bands`, `read_structure` | eigenvalues and structure from the XML file |
| `altharm.basis.CircularBasis`, `SphericalBasis` | the angular basis functions |
| `altharm.symmetry.spin_operations`, `allowed_l`, `symprec_scan` | symmetry operations with their sublattice character, allowed harmonics, tolerance scan |
| `altharm.project.project` | least-squares fit of the harmonics |
| `altharm.noise.noise_report` | noise floor from the forbidden harmonics |
| `altharm.sampling.sampling_report` | check of the number of points on the circle |
| `altharm.bandtrack.track_by_continuity` | follow a band through crossings |
| `altharm.wannier` | the same harmonics from a Wannier90 `hr.dat` model |

## Three routes to the same coefficients

| Route | Status |
|---|---|
| Quantum ESPRESSO, non-self-consistent run on a circle | implemented, validated on Cr₂Se₂As, MnTe and CoF₂ |
| Wannier90 model (`altharm.wannier`) | implemented, validated on Cr₂Se₂As |
| `harmwave.x` inside Quantum ESPRESSO (`qe_native/`) | implemented, agrees with the Python route within 10⁻⁹ eV |

## Units

Energies are in eV inside the package and in meV in the printed report. k-vectors are
Cartesian, in 1/Å.

## Further reading

- [NOTES.md](https://github.com/yucj-520/altharm/blob/master/NOTES.md): the full discussion of the method and its limits.
- [PLAN.md](https://github.com/yucj-520/altharm/blob/master/PLAN.md): the complete development record.
- [`qe_native/README.md`](https://github.com/yucj-520/altharm/blob/master/qe_native/README.md): building and using `harmwave.x`.
