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

**Measured: the interface is validated, the noise floor is not.**
`bands_input` and `read_bands` have been round-tripped through real
Quantum ESPRESSO 7.2 on bcc iron (LSDA, ONCV, 80 Ry, a 24-point circle of radius
0.40 $2\pi/a$). Against the eigenvalues `pw.x` printed — a different code path
inside QE from the XML, parsed in the test by its own regex — the reader agrees
to **5e-5 eV**, which is the four-decimal rounding of the text output and so the
best this comparison can resolve. $E_F$ and `alat` match the inputs. All 24
k-points came back unfolded, at exactly the requested radius, in order, which is
what `nosym`/`noinv` are there for. The run is committed under
`tests/data/qe_fe_bands/`.

**A bug this phase found that no unit test could.** The writer generated an
input `pw.x` rejects outright: `celldm(1)` and `CELL_PARAMETERS angstrom`
together give *"lattice parameter specified twice"* from `cell_base_init`. Every
synthetic test passed on that input. Fixed by writing the cell in `alat` units
and keeping `celldm(1)`, since `celldm(1)` is what defines `tpiba`.

**Measured, twice -- the first measurement was on the wrong band, and the correction is
itself part of the record.**

The first pass used band index 22 (0-based) because that was the "band nearest $E_F$" index
found earlier for the group's *Wannier* model, which has 38 bands. This DFT run has 80
(16 semicore + 6 Se/As-4s + 38 p-d + 20 buffer, exactly matching the group's own comment in
`nscf.in`), so band 22 there is a deep, nearly dispersionless semicore state 4.4 eV below
$E_F$ -- not the band this project cares about, and its splitting is trivially near zero. That
made the noise floor measured on it look spectacular (0.106 μeV) for the wrong reason: there
was almost no signal to be noisy about. Caught only by checking the band's mean energy against
$E_F$, which the first pass never did -- exactly the "band identity, not sorted position"
convention `CLAUDE.md` already names, applied across two *different* band-counting schemes
rather than across k-points on one path.

**Measured correctly, on band 44** (0-based; $E - E_F = +0.0415$ eV, matching the Wannier
route's own reported band to the same energy within 2 meV), same run, same 56-point circle at
$k_z = 0.25\,b_3$:

| | value |
|---|---|
| $l=3$ (allowed) | 59.13 meV — agrees with the independent Wannier fit (59.48 meV) to 0.6% |
| $l=9$ (allowed) | 38.86 meV — agrees with Wannier (39.06 meV) to 0.5% |
| forbidden **even** $l$ (0,2,4,6,8,10,12) | all $<$ 1 μeV — as clean as the wrong-band measurement suggested |
| forbidden **odd, non-multiple-of-3** $l$ (1,5,7,11) | **0.4 – 1.2 meV** — three orders of magnitude worse |

The DFT-vs-Wannier agreement at $l=3,9$ is a real result this project did not have before, and
a strong one: two routes sharing no code, on a genuine altermagnet, agreeing to under 1%. But
the noise-floor picture the phase set out to measure is now openly unresolved rather than
falsely closed: the forbidden harmonics are not uniformly small, and the ones that are large
(0.4-1.2 meV) are not negligible against this same band's own higher allowed harmonics measured
on the Wannier route ($l=15$: 3.2 meV, $l=21$: 8.5 meV) -- a real DFT measurement at those $l$
would need a floor well under a meV to say anything Wannier does not already show.

**Open question this leaves, sharper than before:** is the odd-non-3 floor genuine DFT
noise, or an artefact of resolving up to $l_\mathrm{max}=12$ with only 56 k-points (a fit that
uses 25 of 56 available points as coefficients is not generously over-determined)? P4 below
rules out one candidate mechanism (the plane-wave count's own discontinuity, which is exactly
even in $l$ and so cannot produce odd-$l$ noise). The next step is a denser circle -- enough
k-points to push $l_\mathrm{max}$ past 21 with real degrees of freedom to spare, so the
forbidden floor can be measured with the same rigor the allowed harmonics just were.

**Resolved with a denser circle.** The 56-point circle's odd-non-multiple-of-3 floor
(0.4-1.2 meV) was suspected to be a finite-$k$ artefact rather than fixed DFT noise; a
152-point circle at the same $k_z$, reusing the already-converged density (no scf rerun),
confirms it. At $l_\mathrm{max}=45$:

| | 56 points | 152 points |
|---|---|---|
| noise floor (RMS, forbidden $l$) | 0.53 meV | **27.5 μeV** |
| worst forbidden harmonic | 1.22 meV ($l=11$) | 86.3 μeV ($l=41$) |
| weakest allowed/floor ratio | ~30$\times$ | **36$\times$** |
| even forbidden $l$ | already clean | still clean, unchanged |

Tripling the k-count dropped the odd-forbidden floor by roughly an order of magnitude while
the allowed harmonics did not move (still agreeing with Wannier to under 1%, at *every* one of
$l=3,9,15,21,27,33,39,45$, not only the leading two) and the even-forbidden floor stayed exactly
as clean as before. That is the signature of finite-$k$ aliasing shrinking with more points, not
of a fixed numerical noise floor that would not care how many k-points were used. **P3's
headline number, properly measured this time:** noise floor 27.5 μeV RMS, worst single
harmonic 86 μeV, against the group's ~9 meV self-consistent-on-mesh floor — over two orders of
magnitude better, on the physically relevant band, with the signal it is being compared against
independently confirmed rather than assumed. The residual beyond $l=45$ (0.80 meV) is not part
of the floor — it is real unresolved structure, consistent with the Wannier route needing
$l$ up to ~48 to fully converge this band.

**What remains from the original plan:** the second independent check (re-run from a different
converged density or SCF threshold) has still not been done. Whether the residual odd-forbidden
floor (now 20-90 μeV, scaling down with k-count) extrapolates to zero as $k\to\infty$ or
saturates at a genuine noise level is also not yet measured — it would take a third k-density
to distinguish a $1/N_k$-type falloff from a plateau.

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

**Measured, on the group's real Cr2SeAs circle** (56 k-points, $k_z = 0.25\,b_3$,
`ecutwfc = 60` Ry): the plane-wave count steps between 10770 and 10801 around the circle -- a
peak-to-trough swing of 31 out of ~10800, about 0.3%. Its own angular harmonic spectrum
(`altharm.planewaves.discontinuity_spectrum`), fit on the exact same basis as the splitting:

| $l$ | 2 | 4 | 6 | 8 | 10 | 12 |
|---|---|---|---|---|---|---|
| amplitude | 1.28 | 0.06 | **12.98** | 0.48 | 1.24 | 0.62 |

and **exactly zero at every odd $l$** -- a symmetry of the sphere-counting itself, not a
coincidence of this circle.

**This phase originally claimed the $l=6$ peak here explained P3's largest forbidden
harmonic.** That was correct arithmetic on the wrong band (P3's band-22 mistake, see above):
band 22's one substantial forbidden harmonic happened to sit at $l=6$, matching this spectrum's
own peak, but it was a coincidence between two nearly meaningless numbers, not a causal link.
On the corrected band 44, the real forbidden-harmonic problem is the *opposite* parity --
$l=1,5,7,11$, all odd -- and the plane-wave discontinuity is exactly even in $l$ by
construction. **It is therefore ruled out, not confirmed, as the source of P3's actual noise.**
A mechanism checked and excluded is still a result; it narrows what P3's open question is
asking about.

**Found: mostly finite-$k$ aliasing, not `ecutwfc`.** Tripling the circle's point count
(56$\to$152, same `ecutwfc`, same converged density, P3) dropped the odd-forbidden floor by
roughly an order of magnitude while the even-forbidden floor and every allowed harmonic stayed
put. A discretisation effect tied to `ecutwfc` would not respond to k-count alone; a
finite-sampling aliasing effect would, and did. The plane-wave-count mechanism this phase
originally measured remains real (it is a genuine, `ecutwfc`-dependent, exactly-even-in-$l$
effect) but is now understood to be a minor contributor next to plain undersampling for this
circle's noise budget.

**Not yet done:** the original plan's own ask, repeating at several `ecutwfc` values (at fixed,
adequate k-density) to size the plane-wave-count mechanism's own contribution on its own, now
that it is no longer confounded with undersampling.

## P5 — Convergence in `ecutwfc` and k-density

**Establishes:** that the allowed-$l$ coefficients (for $l \leq l_\mathrm{max}$) have actually
converged in the plane-wave cutoff and in how finely the circle is sampled, as a claim separate
from P3's noise floor and P4's discontinuity — a coefficient can sit above the noise floor and
still not be converged.

**Independent check:** increase `ecutwfc` and the number of circle points independently and
confirm the allowed-$l$ coefficients stop moving while the forbidden-$l$ floor from P3 does not
grow back in.

**Measured, the k-density half:** already answered by P3's 56$\to$152-point comparison. The
allowed harmonics did not move (each still within 1% of the Wannier route) while the forbidden
floor dropped an order of magnitude, which is the required behaviour: converged signal, shrinking
noise.

**Measured, the `ecutwfc` half — the cheap part first, and it changed what the expensive part
should be.** `count_planewaves` reproduces QE's own plane-wave count exactly (P4), so the size
of the basis-set discontinuity at any candidate cutoff can be measured *without running DFT at
all*. Scanning 50-100 Ry in 2 Ry steps on the real Cr2SeAs circle, reporting the $l=6$ amplitude
of the plane-wave count per million plane waves:

| ecutwfc (Ry) | 54 | 60 | 66 | 72 | 80 | 96 | 100 |
|---|---|---|---|---|---|---|---|
| $l=6$ (ppm) | 255 | **1168** | **14.5** | 48 | 975 | 571 | 29 |

**It does not converge smoothly, and it is not monotone.** Neighbouring cutoffs differ by more
than an order of magnitude, because which reciprocal-lattice points happen to land in a thin
shell at the cutoff radius is an arithmetic property of that particular radius, not a quantity
that anneals as the basis grows. Raising `ecutwfc` is therefore the *wrong instinct* for this
noise source: 80 Ry is 67$\times$ dirtier than 66 Ry while costing 33% more plane waves.

**And the production cutoff is the worst one in the range.** 60 Ry — what the group's setup
uses — measures 1168 ppm, the maximum anywhere in 50-100 Ry; 66 Ry measures 14.5 ppm for 15%
more plane waves. That turns P5 from a convergence scan into a **falsifiable prediction**: if
the plane-wave discontinuity is what sets the even-$l$ floor, a 66 Ry run should show that floor
drop by up to ~80$\times$ while the allowed harmonics stay put. If it does not, the mechanism is
not what sets that floor, and P4's "ruled out for odd $l$" should be extended to even $l$ too.

The scan is pinned in `tests/test_p5_cutoff_selection.py`, including the $N_{pw} \propto
E_{cut}^{1.5}$ sphere-volume law as a sanity check on the counter itself (exponent recovered to
within 0.02).

**Running:** scf + 152-point bands at 66 Ry and at 80 Ry, to test that prediction against real
eigenvalues.

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

**The QE↔Wannier agreement itself: measured, and now exhaustive.** On band 44 of the
group's real Cr2SeAs run, a 152-point circle (P3) puts the DFT-fitted allowed harmonics within
1% of this Wannier model's own fit at **every one of the eight allowed harmonics from $l=3$ to
$l=45$**: 59.17/59.48 ($l=3$), 38.82/39.06 ($l=9$), 3.16/3.18 ($l=15$), 8.42/8.44 ($l=21$),
1.49/1.49 ($l=27$), 2.84/2.84 ($l=33$), 0.99/1.00 ($l=39$), 1.15/1.14 meV ($l=45$; DFT/Wannier).
Two routes sharing no code -- a plane-wave calculation and a tight-binding fit to a different
plane-wave calculation, read by two different parsers, projected by the same `altharm.project`
but starting from entirely independent eigenvalue sources -- agreeing to under 1% across eight
harmonics spanning nearly two orders of magnitude in amplitude is the strongest result this
project has produced.

**Still open:** the synthetic-`hr.dat` check that the truncation damping really follows
$J_l(k_0 R_\mathrm{max})$. A weaker form is in `tests/test_wannier.py`: a hand-built model
carrying a known $l=6$ component loses over 98% of it under a cut that leaves $l=2$ above half.

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

## P9 — A QE-native tool: harmonics computed inside a compiled QE executable

**Establishes:** that the harmonic fit and the magnetic-symmetry classification can be computed
directly by a new compiled Quantum ESPRESSO executable, `harmwave.x` (`qe_native/harmwave.f90`),
built the way `bands.x`/`dos.x` are -- linked against QE's own `PP/src`/`PW/src` libraries,
reading a finished `bands` run through QE's own `read_file()`. It uses QE's own `find_sym()`
(with the magnetic sublattices collapsed to their true chemical species via the pseudopotential's
`psd` field, not `ityp`) for symmetry, and LAPACK's `DGELSS` directly on the in-memory Rydberg
eigenvalues for the fit -- no XML parsing, no spglib, no numpy.

**Independent check:** exactly the point of building it. Two implementations sharing no
language, no symmetry engine, no linear-algebra library and no unit-conversion path, run on the
same real Cr2SeAs data, are compared coefficient by coefficient.

**Measured:** built and validated once, on the group's real Cr2SeAs run (PBE+U, USPP, 56-point
circle at $k_z = 0.25\,b_3$, band nearest $E_F$):

- eigenvalue differences $\Delta E(\theta)$: Fortran (QE's `et`, Ry$\to$eV) vs Python (XML,
  Ha$\to$eV) agree to $10^{-9}$ eV, run in run out of the same underlying k-point list.
- the fitted $(a_l, b_l)$: LAPACK `DGELSS` vs `numpy.linalg.lstsq`, agree to $10^{-9}$ eV
  coefficient by coefficient, for every $l = 0..12$.
- the operation list QE's own `find_sym()` derived, fed into `altharm.symmetry.allowed_harmonics`,
  reproduces the $\{3, 9\}$ selection rule exactly where the independently-fitted real signal
  is (allowed amplitudes over 100$\times$ the forbidden floor) -- see
  `tests/test_harmwave_cross_check.py` and `tests/data/harmwave_cr2seas/`.

**A bug this phase found that a clean compile could not.** A first version's time-reversal
closure only doubled operations that had *already* survived the plane-preservation filter for
the given circle. Cr2SeAs's sublattice-exchanging mirror flips $k_z$ on its own and so fails
that filter directly; its time-reversed partner (which restores $k_z$) needs the mirror in a
closure taken over the *full* symmetry group, before filtering -- not after. The bug compiled
cleanly, ran to completion, and reported the allowed set as $\{0,3,6,9,12\}$ instead of
$\{3,9,15,\ldots\}$: over-inclusive, and not obviously wrong from the fit alone, since $l=0,6,12$
still came out numerically tiny. Only feeding the tool's own operator dump back into
`altharm.symmetry.allowed_harmonics()` and comparing the predicted rule against where the real
signal landed caught it.

**Scope of v1, honestly not yet general:** collinear LSDA only; the circle must lie in the
$k_x k_y$ plane at fixed $k_z$ (hardcoded); band identity is given explicitly rather than
tracked (`altharm.bandtrack`'s continuity tracking is not yet ported to Fortran); symmetrization
and the noise-floor split stay in the already-tested Python code rather than being duplicated in
Fortran. See `qe_native/README.md`.

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
