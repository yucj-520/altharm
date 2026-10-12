# altharm

> **Status (2026-09-21).** The workflow has been validated end to end on
> Cr2Se2As, and the symmetry analysis has additionally been checked against the
> literature classifications of MnTe (g-wave) and CoF2 (d-wave). The test suite
> contains 239 tests. The reference coefficients were obtained with a
> **24x24x8 k-mesh, `degauss` = 0.0150, band 44, a 152-point circle at
> $k_z = 0.25b_3$ and $U = 3.0$ eV**: $l=3$: 64.293, $l=9$: 43.447,
> $l=15$: 3.603, $l=21$: 8.662, $l=27$: 1.876, $l=33$: 2.583, $l=39$: 1.182 and
> $l=45$: 0.985 meV, with a noise floor of 23.6 μeV. Each of these values is a
> coefficient *at the stated value of $U$*; see
> ["What the noise floor does not measure"](#what-the-noise-floor-does-not-measure).
> `PLAN.md` contains the complete development record, including the stages at
> which earlier conclusions were found to be incorrect and were withdrawn.


`altharm` extracts the angular harmonics of the spin splitting of an altermagnet,

$$\Delta E(\mathbf k) = \varepsilon_\uparrow(\mathbf k) - \varepsilon_\downarrow(\mathbf k),$$

directly from plane-wave density functional theory (DFT), and expands it on a
circle in k-space (for a monolayer, or for a fixed-$k_z$ cut through a
three-dimensional material) as

$$\Delta E(\theta) = \sum_l \big[a_l \cos(l\theta) + b_l \sin(l\theta)\big],$$

or, for a genuinely three-dimensional splitting, on a sphere in real spherical
harmonics $Y_{lm}$. The leading non-vanishing harmonic classifies the
altermagnet. Note, however, that the conventional d-, g- and i-wave labels count
**nodal surfaces**, not the index $l$; the two coincide only in the absence of a
horizontal nodal plane. Rutile CoF2 is d-wave with $l=2$: it has two vertical
nodal planes. MnTe is **g-wave with $l=3$**: the $\cos 3\theta$ term contributes
three vertical nodal planes, and the $k_z=0$ plane provides a fourth. Reading
"g-wave" as "$l=4$" would therefore make the correct MnTe result obtained here
appear incorrect. Earlier work based on a Wannier tight-binding fit has
established that MON is i-wave, with $\Delta E \propto \cos(6\theta)$, and that
$l=6$ in that case is not an artefact of the fit but a consequence of the
magnetic order breaking the mirror planes of the crystal while preserving its
$C_3$ rotation symmetry.

## Motivation for a direct DFT approach

The harmonics reported so far have been obtained from Wannier models fitted to
DFT. The Cr2Se2As model reproduces the DFT bands to within 0.07 meV, which is a
close fit by any standard. Nevertheless, a Wannier model can omit contributions
that are present in the underlying DFT calculation. Computing the harmonics
directly from DFT therefore provides an independent check on the Wannier results
and, in addition, gives access to the high-harmonic content of the spin
splitting, which the Wannier approach cannot reliably resolve.

The mechanism is specific and is not merely a general concern about truncation.
A hopping at distance $R$ contributes to $\Delta E$ on a circle of radius $k_0$ a
term that the Jacobi-Anger expansion decomposes into a sum over **all** $l$,
weighted by $J_l(k_0 R)$, a Bessel function that decays rapidly once $l$ exceeds
$k_0 R$. A Wannier model that includes hoppings up to a cutoff $R_\mathrm{max}$
therefore does not forbid harmonics with $l > k_0 R_\mathrm{max}$; it suppresses
them exponentially, which at finite numerical precision is indistinguishable
from their absence. Any physical high-$l$ contribution is thus damped by the
truncation precisely because it is of high order. A plane-wave DFT calculation
has no hopping range to truncate and hence no equivalent damping. For the same
reason, it is the only one of the two approaches that can establish whether a
high harmonic is physical.

## The noise floor as the central quantity

In place of truncation, DFT is subject to numerical noise: self-consistency
artefacts, a plane-wave basis that is not perfectly smooth in k, and the other
effects listed under ["Known numerical traps"](#known-numerical-traps) below.
The organising principle of this package is therefore:

**the symmetry-forbidden harmonics measure the noise floor.**

All harmonics up to some $l_\mathrm{max}$ are fitted, not only those allowed by
the magnetic space group. A coefficient at a forbidden $l$ cannot be physical,
since symmetry requires it to vanish exactly; any amplitude it carries in an
actual calculation is therefore the numerical noise of that calculation,
measured on the same footing as the signal with which it is compared. A
coefficient at an *allowed* $l$ is considered reliable only when it lies well
above this floor. This turns the question of whether a high-harmonic signal is
real from a matter of judgement into a quantitative comparison.

## Note on the chemical formula and directory names

The material is **Cr2Se2As**: the unit cell contains 4 Cr, 4 Se and 2 As atoms,
which reduces to Cr2Se2As rather than Cr2SeAs. Earlier versions of this
repository used Cr2SeAs throughout; the name was inherited from a directory named
`Cr2SeAs_wannier_v03` and had not been checked against the atom count. The text
has since been corrected.

**The test fixture directories are still named `qe_cr2seas_*`**, and one path on
disk is still `Cr2SeAs_wannier_v03`. These are identifiers only. Renaming the 22
fixture directories and the 24 test files that reference them would risk
breaking paths without any gain in correctness, so they have been left unchanged
and should be read as labels.

## What the noise floor does not measure

The noise floor described above is an error bar on the **fit**, not on the
coefficient itself, and the two were found to differ by three orders of
magnitude. The following values were measured on Cr2Se2As at $l=3$:

| source of error | size (median over the eight allowed $l$) | is it reflected in the floor? |
|---|---|---|
| Hubbard $U$, per eV | **~21%** | no |
| **the value of $U$ itself** | **production calculations use 3.0 eV; linear response gives 6.39 eV** | no |
| smearing width set too narrow | 12.9% | no |
| **PAW vs. ultrasoft pseudopotential** | **9.7%** (17.8% at $l$=45) | no |
| **PBEsol vs. PBE** | **7.3%** (33.2% at $l$=39) | no |
| atomic positions at relaxation precision (1e-3 Å) | <4% on the coefficients | partially; it degrades the floor by a factor of 7 |
| **the number of points on the circle** | **18% (36% worst case)** | **no; the floor reports _zero_** |
| **the floor itself** | **0.04%** (24 μeV / 64 meV) | — |

The pseudopotential and functional entries are the clearest cases, because they
shift the floor in the **wrong direction**: changing the functional or the
pseudopotential moves the coefficients by 7-10% while *improving* the reported
floor (from 24.3 to 8.2 and 14.3 μeV, respectively). The floor tracks the
unresolved tail of the signal above $l_\mathrm{max}$, not the error. A
calculation can therefore appear cleaner while being less accurate.

These two effects also act specifically on the tail rather than rescaling the
coefficients uniformly: $l=3$ changes by 0.4-1.3%, whereas $l=39$ changes by
12-33%.

**The number of sample points is a different and more serious kind of error.**
The other entries are blind spots; this one is an inversion. The floor is read
from harmonics that are forbidden by symmetry, but *which* operation forbids a
given harmonic determines what its amplitude means. On the 152-point circle, the
67 coefficients forbidden by an operation that the sample grid can represent lie
at 0.0040 μeV, whereas the 15 coefficients forbidden only by C3 (the one
operation that 152 points cannot represent) lie at 24.70 μeV and **constitute**
the published floor. If $N$ is chosen as a multiple of the period of the allowed
set, every operation becomes representable: the floor falls to zero while the
unresolved tail aliases directly onto the reported coefficients, giving a median
error of approximately 18%. Across $N=120\ldots199$, the floor is *monotonically
anti-correlated* with the true error. $N$ should therefore be chosen coprime to
that period; `altharm.sampling.sampling_report` identifies suitable values. The
152-point circle used in this work shares a common factor with the period and is
approximately five times less accurate than a coprime choice at the same
computational cost. See `PLAN.md`, P14.

**The second row is not a sensitivity but a discrepancy.** P12 quantified the
effect of $U$ around 3.0 eV; P19 examined whether 3.0 eV is the appropriate
value and found that the linear-response calculation with QE's `hp.x` gives
**6.386/6.388 eV**, a factor of 2.13 larger. The scan in P12 is monotonic (every
coefficient decreases as $U$ increases; from $U$ = 2 to 5 eV, $l$=3 decreases by
36% and $l$=45 by 61%), so the direction of the error is certain: the published
coefficients are systematically too large. Its magnitude is not known, because
6.39 eV lies outside the 2-5 eV range that was scanned. Two caveats apply to the
linear-response value: the q-mesh was 2x2x1, and `conv_thr_chi` had to be
relaxed to 1e-4. See `PLAN.md`, P19.

Each of these quantities is an *input* to the Kohn-Sham Hamiltonian, whereas the
floor measures what happens once the Hamiltonian is fixed. An accurate summary
of what this package does is therefore: **for a given Hamiltonian, it resolves
$\Delta E(\theta)$ to tens of μeV; it does not establish whether that
Hamiltonian is the correct one.** All five sources have now been measured
(`PLAN.md` P3, P12, P13 and the structure-precision entry).

The practical consequence is that **convergence of the density does not imply
convergence of the harmonics.** The magnetization, the Hubbard occupations and
the total energy all agreed to three significant digits across a change in
smearing that shifted the coefficients by 27%.

The floor nevertheless remains the appropriate tool for its intended purpose:
determining whether $l_\mathrm{max}$ is sufficiently high, whether the circle is
sampled sufficiently densely, and whether a given harmonic is resolved at all.
It improved the resolution of this work by a factor of $10^4$ (from 9 meV to
1 μeV). It should not, however, be quoted as the uncertainty of a coefficient.

## Three independent routes to the same $a_l$

1. **Quantum ESPRESSO**: `pw.x` run non-self-consistently on an explicit circle
   of k-points.
2. **Wannier90**: interpolation from an `hr.dat` file. This is the existing
   baseline and has already been validated on Cr2Se2As.
3. **[defumat](https://github.com/joselado/defumat)**: a differentiable
   plane-wave code written in JAX. Because it is differentiable, it can provide
   $\partial\Delta E/\partial\mathbf k$ exactly rather than by finite
   differences, which avoids some of the discretisation issues that arise on a
   sampled circle.

The three routes share no common machinery: they use different codes and
different approximations, and they have different failure modes. Agreement
between them, harmonic by harmonic, constitutes the validation; disagreement
indicates which route should be distrusted for which $l$.

## Known numerical traps

The following effects are measured by this package rather than assumed to be
negligible. The corresponding numbers are recorded in `PLAN.md`; this section
lists what to watch for.

- **Hubbard calculations restrict the choice of pseudopotentials.** SG15 ONCV
  norm-conserving pseudopotentials have an **empty** `<PP_PSWFC>` block, whereas
  DFT+U requires atomic wavefunctions to construct the Hubbard projectors, so
  `pw.x` stops with
  `determine_hubbard_occ (1): no atomic wavefunctions in pseudopotential file`.
  Any test of pseudopotential sensitivity on a DFT+U system must therefore remain
  within families that include atomic wavefunctions (GBRV ultrasoft, pslibrary
  PAW). See `PLAN.md`, P13.
- **The valence configuration must be checked when changing pseudopotential
  families.** The default pslibrary Se PAW potential (`-dn`) includes the 3d
  states in the valence (16 electrons rather than 6), which would increase the
  electron count of Cr2Se2As from 90 to 130 and silently change which state band
  44 refers to. The `-n` variant has the matching configuration. The CoF2 entry
  in P8 illustrates the consequences when this is not checked.
- **Self-consistency on the band k-mesh is itself a source of noise.**
  Re-converging the density on a line or circle of k-points (VASP `ICHARG=1`)
  produced approximately 9 meV of eigenvalue noise against a splitting of
  approximately 1 meV, i.e. noise several times larger than the signal. A
  non-self-consistent calculation with a frozen density (QE
  `calculation='bands'`, VASP `ICHARG=11`) ensures that every k-point
  experiences the same potential and is expected to remove this source entirely;
  the extent to which it does so is documented in `PLAN.md`, P3.
- **Basis-set discontinuity.** At fixed `ecutwfc`, the number of plane waves
  satisfying $|\mathbf k+\mathbf G|^2 < E_\mathrm{cut}$ changes discretely as
  $\mathbf k$ moves around the circle. This introduces small steps into
  $\varepsilon(\mathbf k)$ that are broadband in $\theta$, and broadband noise
  contaminates specifically the high-$l$ end of the fit, which is the same range
  in which the Wannier truncation argument places the new physics.
- **Angles must be computed from Cartesian $\mathbf k$.** Computing $\theta$ from
  fractional coordinates applies the lattice metric to the angle; for a
  hexagonal cell this distorts a circle into an ellipse and mixes harmonics by a
  factor of $\sqrt3$. This error has occurred in previous work, and for this
  reason `altharm.basis` accepts only Cartesian directions.
- **Band index is not band identity.** Eigenvalues are returned sorted by
  energy, so at a band crossing the "$n$-th band" changes character and
  $\Delta E_n(\theta)$ jumps discontinuously. A fit must therefore follow a band
  by its character rather than by its position in the sorted list.

## Usage

The package currently provides the angular basis and the symmetry projector; a
single-call entry point for fitting a DFT calculation is not yet available. Both
components below run as shown.

**The angular basis** (`altharm.basis`) provides the design matrix on which the
fit is built:

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

**The symmetry projector** (`altharm.symmetry`) determines numerically which
values of $l$ a magnetic structure allows, using the symmetry operations from
spglib, each labelled according to whether it preserves or exchanges the
magnetic sublattices; no character tables are consulted by hand. Applied to a
two-site honeycomb Néel antiferromagnet (opposite moments on the two sites,
related by an inversion centre):

```python
from altharm.symmetry import spin_operations, allowed_l

ops = spin_operations(cell, scaled_positions, numbers, magmoms,
                       include_time_reversal=True)
print(sorted(allowed_l(basis, ops, dirs)))
# []
```

The empty set is the correct result, not a null result: a centrosymmetric
compensated antiferromagnet satisfies $E_\uparrow(\mathbf k) =
E_\downarrow(\mathbf k)$ everywhere by $\mathcal{PT}$ symmetry, so *no* harmonic
survives and the system is not an altermagnet. Reproducing this textbook case
with the same numerical machinery that is subsequently applied to real
altermagnets is the consistency check underlying `PLAN.md`, P1.

## Installation

```bash
pip install -e .
```

Python 3.10 or newer is required. The dependencies (NumPy, SciPy, spglib and
ASE) are installed automatically by `pip`. The `dev` extra adds `pytest`, and
the `plot` extra adds matplotlib.

## Status

The harmonic machinery, the QE route, the Wannier route, and a second,
independent implementation compiled directly into Quantum ESPRESSO have all been
implemented and tested. 99 tests pass.

| | |
|---|---|
| **Implemented and tested** | `units.py`; `basis.py`; `kpath.py`; `symmetry.py` (the numerical projector, together with the `symprec` guard described below); `project.py`; `noise.py`; `bandtrack.py` (band identity by continuity); `planewaves.py` (an independent counter of plane waves inside the cutoff sphere); `wannier.py`; `qe/` (input writer and output reader) |
| **Validated against real data** | the Wannier route, on the Cr2Se2As v03 Wannier model; the QE route, on a production Cr2Se2As PBE+U calculation (see below) |
| **Validated against real `pw.x` output** | bcc iron (interface validation, `tests/data/qe_fe_bands/`) and the production Cr2Se2As calculation (`tests/data/qe_cr2seas_p3/`) |
| **A second, independent implementation** | `qe_native/harmwave.f90`: a new executable compiled directly into Quantum ESPRESSO (linked against its own `PW/src` and `PP/src`, and using its own `find_sym()` and LAPACK, without XML, spglib or NumPy). It agrees with the Python route to within $10^{-9}$ eV on the same data. See `qe_native/README.md`. |
| **Not yet started** | the defumat route; the 3D spherical-harmonic path beyond its unit tests; a complete noise-floor measurement (see below; the current result is inconclusive rather than negative) |

### Measured results

For the Cr2Se2As Wannier model (PBE+U, 38 orbitals per spin, 583
$\mathbf R$-vectors), at $k_z = 0.25\,b_3$ and for the band nearest $E_F$:

- From the structure alone, the symmetry module allows exactly the **odd
  multiples of three**, $l \in \{3, 9, 15, 21, \ldots\}$, and **no harmonics at
  all** in the $k_z = 0$ plane. Both results can be derived by hand: $C_{3z}$
  preserves the magnetic sublattices, so $l$ must be a multiple of 3, while $M_z$
  exchanges them and, combined with time reversal, yields $C_{2z}$ with
  character $-1$, so $l$ must be odd.
- The Wannier Hamiltonian, which contains no information about this prediction,
  gives **59.5 meV at $l=3$** and **39.1 meV at $l=9$**, with every
  symmetry-forbidden harmonic below 1.2 meV. The two largest harmonics are both
  symmetry-allowed. This agreement constitutes the cross-check.
- In the $k_z=0$ plane, where symmetry forbids all harmonics, the model gives
  values **below 0.07 meV**. This value is not physical; it is the intrinsic
  angular error of the Wannier model, and it is the floor that the DFT route must
  improve upon.
- The high harmonics require in-plane hoppings out to **approximately 15 Å**.
  With a 4 Å cutoff, $l=21$ decreases by a factor of six while $l=3$ changes by
  one third. The decay is not monotonic, as predicted by the Bessel-function
  argument above and contrary to what a simple exponential-damping model would
  suggest.

For the **Cr2Se2As PBE+U DFT calculation** (USPP, Grimme-D3, same $k_z$), for
the band nearest $E_F$ ($E-E_F = +0.04$ eV), first on a 56-point circle and
subsequently, to resolve the questions this left open, on a 152-point circle
using the same converged density:

- **DFT reproduces the Wannier fit to within 1% for all eight allowed
  harmonics**, from $l=3$ to $l=45$, using an entirely independent source of
  eigenvalues, and not only for the two leading harmonics. This is the strongest
  cross-check in the project.
- **Noise floor, measured once a sufficient number of k-points was available:**
  27.5 μeV RMS over the forbidden harmonics, with a maximum of 86 μeV for any
  single harmonic. This is more than two orders of magnitude below the
  approximately 9 meV floor obtained with self-consistency on the k-mesh. The
  odd forbidden harmonics on the 56-point circle (0.4-1.2 meV, close to the
  weaker allowed signal of this band) decreased by roughly an order of magnitude
  when the number of k-points was tripled, whereas the allowed harmonics and the
  even forbidden harmonics, which were already small, did not change. This is
  the signature of finite-$k$ aliasing decreasing with sampling density, rather
  than of a fixed noise floor. The basis-set discontinuity associated with
  `ecutwfc` remains *excluded* as the cause of that floor (it is exactly even in
  $l$ by construction) and is now understood to be a minor contribution compared
  with undersampling. See `PLAN.md`, P3, P4 and P6, for the values at each
  intermediate step.

**An error identified by this measurement.** The first attempt used band index
22, copied from the "band nearest $E_F$" index of the Wannier model without
verifying that the two models number their bands in the same way. They do not:
the DFT calculation has 80 bands (16 semicore, 6 Se/As 4s, 38 p-d and 20 buffer
bands), compared with 38 in the Wannier model. Band 22 in the DFT calculation is
therefore a dispersionless semicore state 4.4 eV below $E_F$, not the band of
interest, and its near-zero splitting made the noise floor measured on it appear
exceptionally low (0.106 μeV) for an irrelevant reason. The error was detected
only by comparing the mean energy of the band with $E_F$. This is the same
principle that the package already applies when following a band across
k-points (a band is identified by its character, not by its position in the
sorted list), applied here across the band-numbering schemes of two *different*
models.

### Lessons from running Quantum ESPRESSO

The QE input writer passed every test against hand-written reference strings,
yet still produced an input that `pw.x` rejects: specifying `celldm(1)` together
with `CELL_PARAMETERS angstrom` causes the code to stop with *"lattice parameter
specified twice"*. The cell is now written in `alat` units while retaining
`celldm(1)`, because this parameter defines the meaning of `tpiba` in the
k-point card, and an ambiguous `tpiba` would place the circle at the wrong
radius without any warning. `tests/test_qe_real_output.py` tests both the
generated input and the parsed output.

The same calculation confirms the purpose of `nosym`/`noinv`: all 24 requested
k-points were returned, each at exactly the requested radius and in the
requested order, rather than being folded into the irreducible wedge.

### Important: check the symmetry tolerance

For this structure, spglib identifies the non-magnetic parent structure as
**P3m1 with no sublattice-exchanging operation** for `symprec` $\le 10^{-5}$,
and as **P-6m2 with six such operations** for `symprec` $\ge 5\times10^{-5}$.
The horizontal mirror maps every atom onto its partner to within
$1.0\times10^{-5}$ Å, which equals the default tolerance of spglib. With the
default setting, an altermagnet is therefore silently described as a structure
that cannot be one, and every $l$, including $l=0$, is reported as allowed.
`spin_operations` now explicitly rejects this case and refers the user to
`symprec_scan`. **The symmetry tolerance should be scanned for every new
structure before its harmonics are interpreted.**

## License

GNU General Public License v3 or later; see [LICENSE](LICENSE).
