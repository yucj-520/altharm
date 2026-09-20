# altharm

`altharm` extracts the angular harmonics of an altermagnet's spin splitting

$$\Delta E(\mathbf k) = \varepsilon_\uparrow(\mathbf k) - \varepsilon_\downarrow(\mathbf k)$$

directly from plane-wave DFT, and expands it on a circle in k-space (a monolayer,
or a fixed-$k_z$ cut of a 3D material) as

$$\Delta E(\theta) = \sum_l \big[a_l \cos(l\theta) + b_l \sin(l\theta)\big],$$

or, for a genuinely 3D splitting, on a sphere in real spherical harmonics
$Y_{lm}$. The leading surviving harmonic classifies the altermagnet: $l=2$ is
d-wave, $l=4$ is g-wave, $l=6$ is i-wave. The group has already established,
from a Wannier tight-binding fit, that MON is i-wave — $\Delta E \propto
\cos(6\theta)$ — and that $l=6$ there is not a fit result but a consequence of
the magnetic order breaking the crystal's mirror planes while leaving its
$C_3$ rotation intact.

## Why do this again, from DFT

The group's harmonics so far come from Wannier models fit to DFT — the
Cr2SeAs model reproduces the DFT bands to 0.07 meV, which is a tight fit by
any standard. The reason to also do it directly from DFT, in the supervisor's
words, is that "there were some terms that are not captured by the Wannier …
hence it would be good to have both (especially given that the DFT high
harmonic part would be genuinely new)."

The mechanism is specific, not a vague worry about truncation. A hopping at
radius $R$ contributes to $\Delta E$ on a circle of radius $k_0$ a term that
Jacobi-Anger expansion turns into a sum over **all** $l$, weighted by
$J_l(k_0 R)$ — a Bessel function that falls off sharply once $l$ exceeds
$k_0 R$. A Wannier model built from hoppings out to some $R_\mathrm{max}$
therefore doesn't forbid a harmonic of index $l > k_0 R_\mathrm{max}$, it
suppresses it exponentially — which in finite numerical precision is
indistinguishable from not being there at all. Whatever high-$l$ physics is
real gets damped by the truncation, precisely because it's high-$l$. A
plane-wave DFT calculation has no hopping range to truncate at, so it has no
equivalent damping — which is also why it is the only one of the two routes
that can tell you whether a high harmonic is real.

## The noise floor is the point

What DFT has instead of truncation is noise: self-consistency artefacts, a
plane-wave basis that isn't perfectly smooth in k, and everything else in
["known numerical traps"](#known-numerical-traps) below. So the organising
idea of this package is:

**the symmetry-forbidden harmonics measure the noise floor.**

Fit *every* $l$ up to some $l_\mathrm{max}$, not only the ones the magnetic
space group allows. A coefficient at a forbidden $l$ cannot be physical —
symmetry says it must be exactly zero — so whatever amplitude it carries in a
real calculation is the calculation's own numerical noise, measured on the
same footing as the signal it's being compared against. A coefficient at an
*allowed* $l$ is only believable once it stands well above that floor. This
turns "is the high-harmonic signal real?" from a judgement call into a number
with a number next to it.

## Three routes to the same $a_l$

1. **Quantum ESPRESSO**, `pw.x` run non-self-consistently on an explicit
   circle of k-points.
2. **Wannier90**, interpolating from an `hr.dat` — the group's existing
   baseline, already validated on Cr2SeAs.
3. **[defumat](https://github.com/joselado/defumat)**, the group's
   differentiable JAX plane-wave code. Being differentiable, it can deliver
   $\partial\Delta E/\partial\mathbf k$ exactly rather than by finite
   differences, which sidesteps some of the discretisation issues a sampled
   circle runs into.

These three share no machinery — different code, different approximations,
different failure modes. Agreement between them, harmonic by harmonic, is the
validation; disagreement says which route to distrust for which $l$.

## Known numerical traps

Things this project is built to measure rather than assume away. The numbers
belong in `PLAN.md`; this is the list of what to watch for.

- **Self-consistency on the band k-mesh is itself a noise source.**
  Re-converging the density on a line or circle of k-points (VASP
  `ICHARG=1`) gave this group ≈9 meV of eigenvalue noise against a ≈1 meV
  splitting — noise several times the signal. A frozen-density
  non-self-consistent run (QE `calculation='bands'`, VASP `ICHARG=11`) makes
  every k-point see the same potential and should remove this source
  entirely; how far it actually drops is `PLAN.md` P3.
- **Basis-set discontinuity.** At fixed `ecutwfc`, the number of plane waves
  inside $|\mathbf k+\mathbf G|^2 < E_\mathrm{cut}$ changes discretely as
  $\mathbf k$ moves around the circle. That puts small steps into
  $\varepsilon(\mathbf k)$ that are broadband in $\theta$, and broadband noise
  contaminates the high-$l$ end of the fit specifically — the same end the
  Wannier truncation argument says is where the new physics is.
- **Angles must come from Cartesian $\mathbf k$.** Computing $\theta$ from
  fractional coordinates applies the lattice metric to the angle; for a
  hexagonal cell this stretches a circle into an ellipse and mixes harmonics
  by a factor of $\sqrt3$. This has bitten the group before. `altharm.basis`
  only ever takes Cartesian directions for exactly this reason.
- **Band index is not band identity.** Eigenvalues come back sorted by
  energy, so at a band crossing the "$n$-th band" swaps character and
  $\Delta E_n(\theta)$ jumps discontinuously. A fit has to track a band by
  what it *is*, not by where it sits in the sorted list.

## Usage

What exists today is the angular basis and the symmetry projector — not yet a
one-call "fit my DFT run" entry point. Both pieces below run as shown.

**The angular basis** (`altharm.basis`) is the design matrix a fit will
eventually be built on:

```python
import numpy as np
from altharm.basis import CircularBasis

basis = CircularBasis(lmax=3)
dirs = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])   # Cartesian, unit-norm

print(basis.l_of())          # [0 1 1 2 2 3 3]  -- one column per (l, cos/sin)
print(basis.evaluate(dirs))
# [[ 1.  1.  0.  1.  0.  1.  0.]
#  [ 1.  0.  1. -1.  0. -0. -1.]]
```

**The symmetry projector** (`altharm.symmetry`) derives which $l$ a magnetic
structure allows numerically, from spglib's operations tagged by whether each
one preserves or exchanges the magnetic sublattices — no character table is
looked up by hand. Run on a two-site honeycomb Néel antiferromagnet (opposite
moments on the two sites, related by an inversion centre):

```python
from altharm.symmetry import spin_operations, allowed_l

ops = spin_operations(cell, scaled_positions, numbers, magmoms,
                       include_time_reversal=True)
print(sorted(allowed_l(basis, ops, dirs)))
# []
```

The empty set is the correct answer, not a null result: a centrosymmetric
compensated antiferromagnet has $E_\uparrow(\mathbf k) = E_\downarrow(\mathbf
k)$ everywhere by $\mathcal{PT}$ symmetry, so *no* harmonic survives, and this
case is not an altermagnet at all. Getting that textbook case right, from the
same numerical machinery that will later be pointed at a real altermagnet, is
the sanity check underneath `PLAN.md` P1.

## Install

```bash
pip install -e .
```

Python 3.10 or newer. Dependencies are NumPy, SciPy, spglib and ASE; `pip`
installs them. The `dev` extra adds `pytest`, and `plot` adds matplotlib.

## Status

The harmonic machinery, the QE route, the Wannier route, and a second independent
implementation compiled directly into Quantum ESPRESSO are all written and tested.
99 tests pass.

| | |
|---|---|
| **Written and tested** | `units.py`; `basis.py`; `kpath.py`; `symmetry.py` (the numerical projector, plus the `symprec` guard below); `project.py`; `noise.py`; `bandtrack.py` (band identity by continuity); `planewaves.py` (an independent plane-wave-count sphere counter); `wannier.py`; `qe/` (writer + reader) |
| **Validated against real data** | the Wannier route, on the group's Cr2SeAs v03 model; the QE route, on the group's real Cr2SeAs PBE+U production run — see below |
| **Validated against real `pw.x` output** | bcc iron (interface validation, `tests/data/qe_fe_bands/`) and the group's Cr2SeAs run (production, `tests/data/qe_cr2seas_p3/`) |
| **A second, independent implementation** | `qe_native/harmwave.f90` — a new executable compiled directly into Quantum ESPRESSO (linked against its own `PW/src`/`PP/src`, using its own `find_sym()` and LAPACK, no XML, no spglib, no numpy) that agrees with the Python route to $10^{-9}$ eV on the same real data. See `qe_native/README.md` |
| **Not started** | the defumat route; the 3D spherical-harmonic path beyond its unit tests; a proper noise-floor measurement (see below — the current one is inconclusive, not negative) |

### What has been measured

On the group's Cr2SeAs Wannier model (PBE+U, 38 orbitals per spin, 583
$\mathbf R$-vectors), at $k_z = 0.25\,b_3$ on the band nearest $E_F$:

- The symmetry module, from the structure alone, allows exactly the **odd
  multiples of three**, $l \in \{3, 9, 15, 21, \ldots\}$, and **nothing at all**
  in the $k_z = 0$ plane. Both follow by hand: $C_{3z}$ preserves the magnetic
  sublattices so $l$ must be a multiple of 3, while $M_z$ exchanges them and,
  combined with time reversal, gives $C_{2z}$ with character $-1$, so $l$ must
  be odd.
- The Wannier Hamiltonian, which knows nothing of that prediction, puts
  **59.5 meV at $l=3$** and **39.1 meV at $l=9$**, with every symmetry-forbidden
  harmonic below 1.2 meV. The two largest harmonics are the two that are
  allowed. That agreement is the cross-check.
- In the $k_z=0$ plane, where symmetry forbids everything, the model returns
  **below 0.07 meV**. That number is not physics — it is the Wannier model's own
  angular error, and it is the floor the DFT route has to beat.
- The high harmonics need in-plane hoppings out to **≈15 Å**. Cut to 4 Å and
  $l=21$ falls by a factor of six while $l=3$ moves by a third. The falloff is
  not monotonic, which is what the Bessel picture above predicts and what a
  simple exponential-damping story would miss.

On the group's **real Cr2SeAs PBE+U DFT run** (USPP, Grimme-D3, a 56-point circle
at the same $k_z$), on the band nearest $E_F$ ($E-E_F = +0.04$ eV):

- **DFT reproduces the Wannier fit to under 1%**, from an entirely independent
  eigenvalue source: $l=3$ at 59.13 meV (Wannier: 59.48), $l=9$ at 38.86 meV
  (Wannier: 39.06). This is the strongest single cross-check in the project so
  far, and it very nearly wasn't measured correctly — see the trap below.
- The noise-floor question this route exists to answer is **not yet resolved,
  and is more interesting than a clean number would have been**: every
  forbidden *even* $l$ is clean (sub-μeV), but every forbidden *odd,
  non-multiple-of-three* $l$ (1, 5, 7, 11) sits at 0.4–1.2 meV — not negligible
  against this band's own smaller allowed harmonics at higher $l$ (15, 21,
  measured on the Wannier route at 3.2 and 8.5 meV). The basis-set
  discontinuity from `ecutwfc` is *ruled out* as the cause, precisely because
  it is exactly even in $l$ by construction (measured: its own harmonic
  spectrum peaks at $l=6$, an order of magnitude over its other even
  harmonics, and vanishes exactly at every odd $l$) — a mechanism excluded is
  still a finding. What actually produces the odd-$l$ floor, and whether a
  denser circle would resolve it, is open; see `PLAN.md` P3/P4.

**A trap this measurement caught in itself.** The first attempt used band
index 22 — copied from the Wannier model's own "band nearest $E_F$" index,
without checking that the two models number their bands the same way. They
don't: this DFT run has 80 bands (16 semicore + 6 Se/As-4s + 38 p-d + 20
buffer) against the Wannier model's 38, so band 22 in the DFT run is a
dispersionless semicore state 4.4 eV below $E_F$, not the band this project
cares about — and its near-zero splitting made the noise floor measured on it
look spectacular (0.106 μeV) for a completely uninteresting reason. Caught
only by checking the band's mean energy against $E_F$, which is exactly the
"a band is identified by its character, not its sorted position" convention
this project already had for following one band across k-points, applied
here across two *different* models' band-counting schemes instead.

### What running real QE changed

The QE writer passed every test against hand-built strings and still produced an
input `pw.x` refuses: `celldm(1)` together with `CELL_PARAMETERS angstrom` stops
the code with *"lattice parameter specified twice"*. The cell is now written in
`alat` units, keeping `celldm(1)` — because that is what fixes the meaning of
`tpiba` in the k-point card, and an ambiguous `tpiba` would put the circle at
the wrong radius without saying so. `tests/test_qe_real_output.py` regresses on
the generated input as well as on the parsed output.

The same run confirms the point of `nosym`/`noinv`: all 24 requested k-points
came back, each at exactly the requested radius and in the order given, rather
than folded into the irreducible wedge.

### One trap worth knowing before you use this

On this very structure, spglib reports the non-magnetic parent as **P3m1 with
no sublattice-exchanging operation** for `symprec` $\le 10^{-5}$, and as
**P-6m2 with six of them** for `symprec` $\ge 5\times10^{-5}$. The horizontal
mirror maps every atom onto its partner to within $1.0\times10^{-5}$ Å, which
is spglib's default tolerance — so the default silently describes an
altermagnet as something that cannot be one, and every $l$ including $l=0$
comes out allowed. `spin_operations` now refuses that case by name and points
at `symprec_scan`. **Scan every new structure before believing its harmonics.**

## License

GPL v3 or later — see [LICENSE](LICENSE).
