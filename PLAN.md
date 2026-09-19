# PLAN.md

This is the record. `CLAUDE.md` is the rules — conventions, refusals, what a finished phase
has to leave behind. This file is what each phase actually found: what it establishes, what
independent check settles it, and the number that check produced. A phase that has code but no
entry here is not done, whatever the code does.

Every phase's "measured" line is either a real number with the check that produced it, or the
literal words `not yet` — never a description standing in for a number that could have been
run.

## P0 — Basis orthonormality

**Establishes:** that `CircularBasis` and `SphericalBasis` (`altharm/basis.py`) are genuinely
orthonormal under the measure a least-squares projection will assume ($\int_0^{2\pi}
d\theta$ on the circle, $d\Omega = \sin\theta\, d\theta\, d\phi$ on the sphere), so that fitting
against them recovers $a_l, b_l$ without a separate, easy-to-forget normalisation step.

**Independent check:** the Gram matrix $\langle\phi_i,\phi_j\rangle$ computed by numerical
quadrature that does not go through `basis.py`'s own machinery for anything but evaluating the
functions themselves — a fine trapezoidal grid in $\theta$ for the circle, a Gauss-Legendre
grid in $\cos\theta$ times a uniform $\phi$ grid for the sphere. An integration routine is
independent of the object it integrates.

**Measured:** `SphericalBasis` is orthonormal — Gram matrix diagonal 1.0000 to $1.4\times
10^{-14}$, all off-diagonal entries below that, on a $40\times40$ Gauss-Legendre$\times$uniform
grid up to $l=3$. **`CircularBasis` is not** — its Gram matrix under $\int_0^{2\pi}$ came out
$\mathrm{diag}(2\pi,\ \pi,\pi,\ \pi,\pi,\ \pi,\pi)$ for $l=(0,1,1,2,2,3,3)$ on a
200001-point trapezoidal grid (off-diagonal below $3\times10^{-5}$). It is orthogonal, not
orthonormal: the module docstring's claim holds for the spherical basis only. **Open:** either
rescale `CircularBasis.evaluate` by $1/\sqrt{2\pi}$ at $l=0$ and $1/\sqrt\pi$ at $l\geq1$, or
narrow the docstring's claim to "orthogonal" — and whichever is chosen, P2's projection step
must not assume unit norm from this basis until it's settled.

## P1 — Symmetry projector

**Establishes:** which angular harmonics a given magnetic structure allows, from the numerical
projector in `altharm/symmetry.py` (`spin_operations` tags each operation of the non-magnetic
parent by whether it preserves or exchanges the magnetic sublattices; `allowed_harmonics`
diagonalises the resulting projector block by block in $l$). The point of deriving this
numerically rather than from a looked-up character table is that it holds for whatever cell is
actually handed in — strained, distorted, lower-symmetry — not only the ideal one a table would
cover.

**Independent check:** a selection rule derived by hand from a point group's character table,
for a case simple enough to do that reliably (e.g. a $C_3$-symmetric structure should forbid
every $l$ not a multiple of 3), compared against `allowed_harmonics`'s output for the same
structure. `S^2 = S` (the projector is idempotent) is a necessary but not sufficient check and
belongs alongside the character-table comparison, not instead of it.

**Measured: closed.** Three structures whose answers are unambiguous, each checked against a
rule worked out by hand in the test's own docstring (`tests/test_symmetry.py`):

| structure | expected | `allowed_harmonics` returned |
|---|---|---|
| rutile, the two metal sites antialigned | $l=2$, and the function $\sin 2\theta$ i.e. $k_xk_y$ | $\{2,6\}$ up to $l_{max}=8$; the $l=2$ subspace is one-dimensional and is $\sin 2\theta$ to $10^{-8}$ |
| sublattices related by a lattice translation | nothing — conventional AFM | empty set |
| sublattices related by inversion | nothing — conventional AFM | empty set |

The rutile case is the one the phase asked for: a genuine altermagnet, a non-trivial allowed
set, and the $d$-wave form recovered rather than assumed. $S^2=S$ holds to $2\times10^{-16}$,
and every Cartesian rotation is orthogonal to $10^{-10}$ on a hexagonal cell, which is the
guard against the metric convention being wrong.

On the group's real Cr2SeAs the allowed set at $k_z \ne 0$ came out
$\{3, 9, 15, 21, 27, 33, 39, 45\}$ — the odd multiples of three — and **nothing at all** in
the $k_z = 0$ plane. Both follow by hand: $C_{3z}$ preserves the sublattices so $l$ must be a
multiple of 3; $M_z$ exchanges them and, combined with time reversal, gives $C_{2z}$ with
character $-1$, so $l$ must be odd; and $M_z$ acting alone as $k_z \to -k_z$ forces
$\Delta E = -\Delta E$ at $k_z=0$.

## P2 — Synthetic end-to-end recovery of known $a_l$

**Establishes:** that the pipeline from a basis's design matrix through a least-squares fit
recovers coefficients that were put in by hand, before any DFT data is involved. This is the
phase that would catch a sign error, a missing normalisation (see P0), or an off-by-one in
which column of the design matrix belongs to which $l$.

**Independent check:** build $\Delta E(\theta)$ (or $\Delta E$ on a sphere) from an explicit,
hand-chosen set of $a_l, b_l$ — values the fitting code never sees — sample it on a grid of
directions, add controlled synthetic noise, and confirm the fit recovers the input coefficients
to the accuracy least-squares theory predicts for that noise level. The check is only
independent if the same code that generates the synthetic data isn't reused to do the fit, or
vice versa.

**Measured:** not yet. Blocked on the projection step itself, which doesn't exist yet —
`altharm/basis.py` has the design matrix and `altharm/symmetry.py` has the projector, but
nothing currently turns a design matrix and a set of $(\text{direction}, \Delta E)$ pairs into
fitted coefficients.

## P3 — QE noise floor on a real altermagnet (the big feasibility question)

**Establishes:** how far a frozen-density, non-self-consistent QE run
(`calculation='bands'` from an already-converged charge density) reduces the eigenvalue noise
the group measured from re-converging self-consistently on a band/circle k-mesh — a VASP
`ICHARG=1` run gave ≈9 meV of noise against a ≈1 meV altermagnetic splitting there, noise
several times the signal. This phase decides whether the DFT route is usable at all, before any
of the later phases are worth running.

**Independent check:** fit every $l$ up to $l_\mathrm{max}$ on a real material's circle of
k-points; P1's projector says which $l$ are symmetry-forbidden, and their fitted amplitude *is*
the noise floor, measured on exactly the same footing as the allowed-$l$ signal it will be
compared against. A second, independent estimate of the same floor: run the frozen-density
calculation twice from two different converged starting densities (or two SCF thresholds) and
see whether the forbidden-$l$ amplitude moves.

**Measured:** not yet.

## P4 — Basis-set discontinuity sizing

**Establishes:** the size of the steps that a fixed `ecutwfc`'s discretely-changing plane-wave
count puts into $\varepsilon(\mathbf k)$ as $\mathbf k$ moves around the circle, and which $l$
they land on. The expectation, stated in `README.md`, is that this is broadband in $\theta$ and
so should contaminate the high-$l$ end specifically — the same end the Wannier-truncation
argument says carries the new physics, which is why this has to be measured rather than assumed
small.

**Independent check:** repeat the same circle at two different `ecutwfc` values (or two
different circle radii at one `ecutwfc`) and see whether the resulting artefact in $l$-space
moves the way a plane-wave-counting argument predicts — a genuinely physical high-$l$ signal
would not move with `ecutwfc` at all.

**Measured:** not yet.

## P5 — Convergence in `ecutwfc` and k-density

**Establishes:** that the allowed-$l$ coefficients (for $l \leq l_\mathrm{max}$) have actually
converged in the plane-wave cutoff and in how finely the circle is sampled, as a claim separate
from P3's noise floor and P4's discontinuity — a coefficient can sit above the noise floor and
still not be converged.

**Independent check:** increase `ecutwfc` and the number of circle points independently and
confirm the allowed-$l$ coefficients stop moving while the forbidden-$l$ floor from P3 does not
grow back in.

**Measured:** not yet.

## P6 — Wannier cross-check and the truncation study

**Establishes:** two things at once. First, that the QE route and the group's existing
Wannier90 route (`hr.dat` interpolation, already validated to 0.07 meV against DFT on Cr2SeAs)
agree on the allowed-$l$ coefficients, each within the other's noise floor. Second — the part
that's actually new — exactly how a finite hopping range $R_\mathrm{max}$ damps a chosen
high-$l$ harmonic, turning the Bessel-function suppression argument in `README.md` into a
measured curve rather than an assertion.

**Independent check:** the QE/Wannier agreement is itself two codes sharing no implementation.
The truncation study's own check: build a *synthetic* `hr.dat` containing a known harmonic at a
controlled radius, truncate it at a series of $R_\mathrm{max}$ values, and confirm the measured
damping matches the $J_l(k_0 R_\mathrm{max})$ prediction rather than some other falloff.

**Measured: the truncation half is closed; the QE half is not started.**
`altharm/wannier.py` reads `hr.dat` and the study has been run on the group's v03 Cr2SeAs
model (38 orbitals per spin, 583 R-vectors, reach 42.6 Å), at $k_z = 0.25\,b_3$ on the band
nearest $E_F$. Cutting **in-plane** $|R|$ while leaving the c axis intact:

| $R_\mathrm{in}$ (Å) | $n_R$ | $l=3$ | $l=9$ | $l=21$ |
|---|---|---|---|---|
| 4 | 29 | 79.80 meV (+34 %) | 21.35 meV (−45 %) | 1.33 meV (−84 %) |
| 7 | 53 | 91.32 meV (+54 %) | 29.12 meV (−25 %) | 2.70 meV (−68 %) |
| 11 | 125 | 62.93 meV | 37.73 meV | 8.15 meV |
| 15 | 245 | 58.91 meV (−0.9 %) | 38.98 meV (−0.2 %) | 8.50 meV (+0.7 %) |
| 42.6 (full) | 583 | 59.48 meV | 39.06 meV | 8.44 meV |

So the high harmonics need in-plane hoppings out to **≈15 Å** to converge to 1 %, and a
short-ranged model of the kind one would write down on physical grounds (4–7 Å) gets $l=21$
wrong by a factor of six and $l=3$ wrong by a third. The falloff is **not monotonic** — it
overshoots before settling — which is consistent with the Bessel-function picture in
`README.md` and rules out reading it as simple exponential damping.

**A trap this phase uncovered.** Cutting on the isotropic $|R|$ instead is meaningless for this
material: $c = 16.9$ Å, so any isotropic radius below that deletes every $R_z \ne 0$ hopping at
once, the model loses its $k_z$ dependence, and a splitting that lives only at $k_z \ne 0$
collapses to zero. It looks exactly like truncation killing the harmonics. `truncate` therefore
takes a `plane_normal`, and the circle's plane should always be passed to it.

**Still open:** the QE↔Wannier agreement itself (needs P3), and the synthetic-`hr.dat` check
that the damping really follows $J_l(k_0 R_\mathrm{max})$. A weaker form of the latter is in
`tests/test_wannier.py`: a hand-built model carrying a known $l=6$ component loses over 98 % of
it under a cut that leaves $l=2$ above half.

## P7 — defumat cross-check

**Establishes:** agreement between the QE/Wannier pair and
[defumat](https://github.com/joselado/defumat)'s differentiable route, which — being able to
give $\partial\Delta E/\partial\mathbf k$ exactly rather than by finite differences on a sampled
circle — offers a genuinely independent way to get at the same harmonics, and a way to check
the *fit itself* (a Taylor expansion in $\mathbf k$ against the Fourier coefficients recovered
from a sampled circle, at a point where both should agree).

**Independent check:** three codes that share no machinery (Quantum ESPRESSO's plane-wave
SCF/NSCF, Wannier90's tight-binding interpolation, defumat's JAX plane-wave implementation)
converging on the same $a_l$ is the strongest evidence available in this project. The
derivative-based check above is a second, internal-to-defumat cross-check on top of that.

**Measured:** not yet. The defumat route is not started.

## P8 — Production on the group's materials

**Establishes:** the final $a_l$ (and $b_l$) for the group's actual targets — MON, and whatever
follows Cr2SeAs — each reported together with its noise floor (P3, P4) and its cross-route
agreement (P6, P7), so that a published coefficient is a number with an error bar rather than a
number alone.

**Independent check:** the three-route agreement carries over from P6/P7, plus, where it
exists, comparison against the group's own earlier Wannier-only literature values — MON's
established $\cos(6\theta)$ rule is the obvious first target to reproduce before trusting
anything this route finds that Wannier didn't.

**Measured:** not yet.

## Numerical traps, indexed by the phase that has to deal with each

These are described with their physical reasoning in `README.md`; this is where the number
belongs once a phase measures it.

- **Self-consistency noise on the band k-mesh** (≈9 meV against a ≈1 meV splitting, measured on
  a related VASP `ICHARG=1` calculation, not yet on this project's own QE route) — P3.
- **Basis-set discontinuity** from `ecutwfc`'s discretely-changing plane-wave count — P4.
- **`symprec` deciding the physics.** *Measured, and the most dangerous one found so far.* On
  the group's own Cr2SeAs structure, spglib reports the non-magnetic parent as **P3m1 (#156)
  with no sublattice-exchanging operation at all** for `symprec` ≤ 1e-5, and as **P-6m2 (#187)
  with six of them** for `symprec` ≥ 5e-5. The horizontal mirror maps every atom onto its
  partner to within **1.0e-5 Å**, which is spglib's default tolerance — so *the default returns
  the wrong physics silently*. Without an exchanging operation every $l$ including $l=0$ comes
  out allowed, i.e. a compensated altermagnet would be reported as though it carried a net
  splitting. `spin_operations` now refuses that case by name and points at
  `symprec_scan`; the cliff is regression-tested in `tests/test_cr2seas_regression.py`.
  **Every new structure must be scanned before its harmonics are believed.**
- **Angle convention.** Fractional-coordinate angles apply the lattice metric and mix harmonics
  by $\sqrt3$ on a hexagonal cell. Structurally prevented rather than merely tested:
  `altharm.basis` takes only Cartesian directions. Still worth a regression test the day
  anything reads k-points out of a QE, Wannier90 or defumat output file, since that's the
  boundary where a fractional coordinate could sneak back in.
- **Band index vs. band identity.** A crossing swaps which sorted index carries which band's
  character, and $\Delta E_n(\theta)$ jumps if a fit follows the index instead of the band. No
  phase above has real band-structure data yet (P3 is the first), so this is currently a risk
  to design around rather than a measured number — the QE/Wannier interface should decide its
  band-tracking method before P3 is run, not after a jump is found in the output.
