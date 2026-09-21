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

**Measured: closed.** `altharm/project.py` exists and `tests/test_project.py` does exactly
what this phase asked for: $\Delta E(\theta)$ built from $a_2=-0.4$, $b_2=1.3$, $b_6=0.75$
chosen by hand, fitted with $l_\mathrm{max}=8$, and every coefficient recovered to
$10^{-12}$ with no spurious weight anywhere else. The amplitude is also checked to be
invariant under rotating the angular origin (a physical quantity must not depend on where
$\theta$ is measured from), the least-squares route is checked against an independent
quadrature sum on a sphere, and the weighted fit against normal equations solved separately
in the test. Underdetermined and mismatched-length inputs are refused by name rather than
returning garbage.

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

**Confirmed at 304 points, and the floor is now 1 μeV.** The prediction that follows from
the aliasing diagnosis was tested directly: a 304-point circle (Nyquist $l=152$), bands-only on
the same converged density, gives

| | 152 points | 304 points |
|---|---|---|
| noise floor | 27.49 μeV | **1.01 μeV** |
| worst forbidden harmonic | 86.26 μeV | **2.81 μeV** |
| $l=3$ / $l=45$ | 59.1712 / 1.1468 meV | 59.1714 / 1.1466 meV |

A factor of **27 for twice the k-points**, with every allowed harmonic moving by less than
0.2 μeV — the control that makes it aliasing rather than anything else. The falloff is steeper
than a naive expectation because what aliases in is the signal's own high-$l$ tail: moving
Nyquist from 76 to 152 replaces content from $l\sim30$–76 with content from $l>152$, where the
spectrum is far smaller.

**The practical statement this phase ends on:** the floor is not an intrinsic limit of the DFT
calculation. Basis size, cutoff, diagonalisation threshold and density convergence have each
been measured to do nothing for it; k-point count sets it, and buys it down steeply and
cheaply — the 304-point bands step was about fifteen minutes on 16 cores, reusing the density.
The weakest allowed harmonic ($l=39$, 0.99 meV) goes from ~36$\times$ the floor to
~1000$\times$, which is the headroom needed to go after $l>45$.

**A side result that unblocked P7.** Removing `vdw_corr = 'grimme-d3'` leaves every eigenvalue
**bit-identical**: Grimme's correction depends only on nuclear positions and never enters the
Kohn-Sham Hamiltonian at fixed geometry. defumat does not implement D3, so this is what lets
the cross-check run without it rather than rebuilding every result under D2 — which would have
made two independent D2 implementations part of the comparison.

**The second independent check, now done -- and it is a caveat, not a reassurance.**
The same 152-point circle was re-run from a *different* converged density: a deliberately
coarse 6x6x2 SCF mesh (14 irreducible k-points) against the production 12x12x4 (69). Both
converged to 1e-10 Ry; both then frozen for an identical bands step.

| $l$ | 3 | 9 | 15 | 21 | 27 | 33 | 39 | 45 |
|---|---|---|---|---|---|---|---|---|
| 12x12x4 (meV) | 59.171 | 38.817 | 3.162 | 8.421 | 1.486 | 2.838 | 0.986 | 1.147 |
| 6x6x2 (meV) | 54.246 | 33.135 | 2.847 | 7.885 | 1.222 | 3.154 | 0.821 | 1.478 |
| difference | -8.3% | -14.6% | -10.0% | -6.4% | -17.8% | +11.1% | -16.8% | **+28.9%** |

The noise floor also moved, 27.49 -> 47.41 μeV. These are not two samples of one quantity:
the coarse mesh is a genuinely different self-consistent state ($E_F$ 8.1169 -> 8.1379 eV,
band 44 at +0.0242 eV instead of +0.0415 eV above it). The bands moved, so the harmonics moved.

**This is now the largest known uncertainty in the pipeline, and the forbidden-harmonic floor
does not see it.** The floor measures the *fit*; it says nothing about whether the density the
fit was performed on is converged. At $l=3$ the density effect (4.9 meV) is over 100x both
P5's per-coefficient cutoff spread (~20 μeV) and P3's own noise floor (27 μeV).

**What this opens:** the production 12x12x4 mesh has never itself been shown converged. Going
6 -> 12 changed the answer by up to 29%, so 12 -> 24 may too. Until that is run, every
coefficient in P3/P6 carries an unquantified SCF-mesh uncertainty, and the <1% DFT-Wannier
agreement of P6 should be read as "these two routes agree *on this density*", not as an
absolute accuracy claim. Pinned in `tests/test_p3_scf_density.py`.

**The 24x24x8 density, and the answer is worse than the question.** The production mesh was
not merely unverified -- the coefficients do not converge at all over the range tested.

| $l$ | 3 | 9 | 15 | 21 | 27 | 33 | 39 | 45 |
|---|---|---|---|---|---|---|---|---|
| 6x6x2 (14 k) | 54.246 | 33.135 | 2.847 | 7.885 | 1.222 | 3.154 | 0.821 | 1.478 |
| 12x12x4 (69 k) | 59.171 | 38.817 | 3.162 | 8.421 | 1.486 | 2.838 | 0.986 | 1.147 |
| 24x24x8 (449 k) | 64.364 | 43.356 | 3.566 | 8.681 | 1.880 | 2.607 | 1.196 | 1.002 |

$l=3$ runs 54.2 -> 59.2 -> 64.4 meV, monotone, with the increment *growing* (+4.93, +5.19).
On $\Delta E(\theta)$ itself the change is 5.36 meV (6->12) then 4.90 meV (12->24): a ratio of
**1.09**, where a converging sequence would halve. Summed over the eight allowed harmonics,
12.54 meV then 11.37 meV. **There is no convergence in this range.**

**What has been excluded.** The band is the same one throughout (band 44 sits at +0.0242 /
+0.0415 / +0.0372 eV above $E_F$; the 43/44 gap is 11.6 / 18.0 / 19.3 meV and never closes, so
no reordering). The inputs differ only in the k-mesh and in `vdw_corr`, which was separately
measured to leave eigenvalues bit-identical. The calculation is numerically stable: deep
semicore band 20's k-dependent difference between the 12 and 24 runs is **0.003 meV**. And the
density is converged by every usual measure -- absolute magnetization 15.95 $\mu_B$ in all
three, Hubbard occupations agreeing to $10^{-3}$ (4.49067 / 4.49213 / 4.49170), total energies
within 3.7 meV.

So a density that looks converged by moment and energy still moves the valence eigenvalues by
1.5-2.9 meV in their k-dependent part, and that is 4-10% of the splitting and 9-27% of an
individual weak harmonic.

**Resolved: it was Fermi-surface sampling, and the fix works.** Bands 43 and 44 straddle $E_F$
(-0.13 and +0.04 eV) while the production Gaussian smearing is only 0.0037 Ry = 50 meV, so
their occupations track how finely the Fermi surface is sampled. Repeating both meshes at
`degauss = 0.015` Ry:

| | change in $\Delta E(\theta)$, 12 -> 24 | worst harmonic |
|---|---|---|
| degauss 0.0037 (production) | 4.903 meV | 26.5% |
| degauss 0.0150 | **0.541 meV** | **2.2%** |

A factor of **9** on the splitting and **12** on the worst coefficient. Two independent knobs
then land on the same answer, which is what makes it the converged one:

| degauss / mesh | 0.0037 / 12x12x4 | 0.0150 / 12x12x4 | 0.0037 / 24x24x8 | 0.0150 / 24x24x8 |
|---|---|---|---|---|
| $l=3$ (meV) | **59.17** | 63.75 | 64.36 | 64.29 |
| $l=9$ (meV) | **38.82** | 42.92 | 43.36 | 43.45 |

The two 24x24x8 runs differ 4x in smearing yet agree to 0.1-1.7%, so the smearing's own bias is
negligible there: **the converged $l=3$ is 64.3 meV**. The noise floor is unaffected throughout
(27.5 / 24.3 / 25.1 / 23.6 μeV).

**What this costs the group's current setup.** 12x12x4 at degauss 0.0037 underestimates $l=3$
by **8%** and the weak harmonics by over 20%. The cheap fix is `degauss = 0.015` at the same
mesh, which lands within 1-2% of converged for the price of a normal run; 24x24x8 at either
smearing is the reference. Pinned in `tests/test_p3_smearing.py`.

**A caveat that survives.** This is a *metal* with two bands within 130 meV of $E_F$; any
altermagnet with that band topology will have the same problem, and the usual convergence
indicators will not warn about it. Every new material needs this two-knob check, not a
single-mesh run.

**What this does and does not undermine.** It does *not* touch the code-to-code results: QE,
Wannier and defumat are each compared at fixed density, and P7's 0.12% agreement and the
identical 47.4 μeV floor stand. It does *not* touch the noise floor, which measures the fit.
What it undermines is any claim about the **absolute value** of a coefficient: until this is
resolved, every number in P3/P6/P7 should be read as "on this density", and the weak harmonics
($l \geq 27$, 1-2 meV) carry an SCF-mesh uncertainty of order 20%, not 27 μeV. Pinned in
`tests/test_p3_scf_density.py`. Whether the residual odd-forbidden
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

**Done, and the answer is that it contributes nothing measurable (P5).** Repeating at three
`ecutwfc` values at fixed, adequate k-density — chosen so the *predicted* discontinuity spans
80$\times$ (1167.5 / 14.5 / 975.2 ppm at 60 / 66 / 80 Ry) — produced noise floors of
27.5 / 27.0 / 28.7 μeV, a 6% range, with the cutoff predicted to be cleanest measuring
marginally *dirtier* than the production one. The plane-wave count's angular structure is real
and `ecutwfc`-dependent, but it does not reach the eigenvalues at this level. **The mechanism is
therefore ruled out for the whole floor, not merely its odd-$l$ part** — this phase's
contribution to the noise budget is, as far as this measurement can tell, zero.

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

## P6b — The corrected reference values, and what the old settings cost

P3's smearing result invalidates the *absolute* coefficients everywhere else in this project.
Every number in P3, P6 and P7 above was computed at `degauss = 0.0037`. The reference set is
now **24x24x8 at `degauss = 0.0150`**, band 44, 152-point circle:

| $l$ | 3 | 9 | 15 | 21 | 27 | 33 | 39 | 45 |
|---|---|---|---|---|---|---|---|---|
| **reference (meV)** | **64.293** | **43.447** | **3.603** | **8.662** | **1.876** | **2.583** | **1.182** | **0.985** |
| 12x12x4 @ 0.0150 | 63.746 | 42.921 | 3.533 | 8.644 | 1.837 | 2.615 | 1.168 | 1.003 |
| old (12x12x4 @ 0.0037) | 59.171 | 38.817 | 3.162 | 8.421 | 1.486 | 2.838 | 0.986 | 1.147 |
| Wannier (fitted to old) | 59.481 | 39.064 | 3.175 | 8.443 | 1.495 | 2.838 | 0.996 | 1.142 |

Three statements follow, and the third is the one that matters for the group's published work.

1. **The correction is large and not a rescaling:** median 12.9%, worst 26.3% ($l=27$), with
   $l=33$ and $l=45$ moving the *other* way. No single factor fixes an old number.
2. **The affordable setting is good enough:** 12x12x4 at 0.0150 is within 1.2% median / 2.1%
   worst of the reference, for the cost of an ordinary run.
3. **The Wannier model inherits the error.** It was fitted to DFT at the old settings, so it
   still agrees with the *old* DFT to ~1% and sits 12-26% from the corrected reference. P6's
   "<1% at all eight harmonics" was therefore agreement on a **shared error**, not evidence of
   absolute accuracy. Correcting the Wannier side means refitting the wannierisation, which is
   the group's workflow rather than this repo's.

**And the noise floor paid for none of it** -- 27.49 / 24.30 / 23.61 μeV across the three, a
range of under 20% while the signal moved 26%. That is precisely why it never warned us, and it
is the cleanest statement of the method's blind spot: **the forbidden-harmonic floor measures
the fit, not whether the band structure handed to the fit is converged.** Pinned in
`tests/test_p6_corrected_settings.py`.

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

**Measured:** not yet, but the route is no longer blocked — four obstacles were found and
cleared, and each is worth recording because none was visible from the README.

1. **JAX was already here.** Not in the main environment, but in `pyqula_env` (0.6.2, CPU
   device working, `jax.grad` verified). Only `equinox` was missing; installed.
2. **The pseudopotentials are UPF v1, defumat reads only v2.** Its error message names the
   fix, and the converter is one we already built: `q-e/upflib/upfconv.x -u`. Cr and Se
   converted; As was already v2. Converted copies live in
   `altharm_runs/defumat_p7/pseudo_v2/`.
3. **defumat does not implement Grimme-D3**, which the group's production input uses (it has
   D2 only). Rather than rebuild everything under D2 — which would have made two independent
   D2 implementations part of the comparison — this was settled by measurement: a QE bands run
   with `vdw_corr` removed gives **bit-identical eigenvalues**. Grimme's correction depends
   only on nuclear positions and never enters the Kohn-Sham Hamiltonian at fixed geometry, so
   the cross-check can simply omit it. Recorded in P3 and `tests/test_p3_k304.py`.
4. **`/tmp` is node-local on this cluster.** The defumat source sat in a scratchpad under
   `/tmp`; every probe on the login node worked and every SLURM job failed with
   `ModuleNotFoundError`. Anything a batch job touches has to live on `/scratch` or `$HOME`.

**Cost, measured rather than guessed.** On the real Cr2SeAs cell, defumat's SCF costs
$37.8\,\mathrm{s} + 6.8\,\mathrm{s} \times N_k$ per iteration (fit to 1, 2 and 4 k-points;
the 80% fixed part is density, FFT, ultrasoft augmentation and Hubbard, none of which scales
with $k$). For QE's 69 irreducible k-points that is ~8.4 minutes an iteration, ~5 hours for the
35 iterations QE needed — feasible overnight. The band step on the 152-point circle turned out
to be the larger cost and is being timed separately.

**Measured, and P7 is closed for the harmonic comparison.** defumat ran the full pipeline on
the real Cr2SeAs cell: SCF to convergence (45 iterations, 43.8 min) then the 152-point circle
(135.3 min). Compared against a QE run whose input file is **byte-identical** (6x6x2 SCF mesh,
no vdW) -- which matters more here than it normally would, since P3's second check showed the
SCF mesh alone moves these coefficients by 6-29%, so comparing defumat-on-6x6x2 against
QE-on-12x12x4 would have hidden a code difference inside a density difference.

| $l$ | 3 | 9 | 15 | 21 | 27 | 33 | 39 | 45 |
|---|---|---|---|---|---|---|---|---|
| QE 6x6x2 (meV) | 54.2458 | 33.1347 | 2.8469 | 7.8846 | 1.2220 | 3.1537 | 0.8207 | 1.4780 |
| defumat (meV) | 54.2676 | 33.1525 | 2.8439 | 7.8879 | 1.2222 | 3.1536 | 0.8218 | 1.4774 |
| deviation | +0.04% | +0.05% | -0.10% | +0.04% | +0.02% | -0.00% | +0.12% | -0.04% |

**Worst case 0.12%, median 0.04%**, across amplitudes spanning two orders of magnitude. The
SCF free energies agree to all six digits defumat prints (-912.52758604 vs -912.527586 Ry,
0.001 meV), and band 44 sits 0.01 meV from where QE puts it, so this is the same band and not
a lucky index.

**The two codes also report the same noise floor: 47.41 vs 47.43 μeV.** That is the more
interesting half. A floor reproduced to 0.04% by an implementation sharing no linear algebra,
no FFT, no ultrasoft augmentation and no Hubbard code is a property of the *k-sampling and the
physics*, not an artefact of either code -- which is exactly what P3's aliasing diagnosis
predicted and what no single-code measurement could establish.

**The three routes, all on Cr2SeAs band 44:** QE (Fortran plane-wave), Wannier90 (tight-binding
interpolation of a *different* plane-wave run) and defumat (Python plane-wave) agree at every
one of the eight allowed harmonics. QE-vs-defumat is ~10x tighter than QE-vs-Wannier (0.04% vs
0.5%), as it should be: the first pair solve the same equations, the second pair includes an
interpolation. Pinned in `tests/test_p7_defumat.py`.

**Confirmed on the production density too.** The 12x12x4 defumat run finished (SCF 51
iterations / 220 min, bands 125 min). Its SCF free energy reproduces QE's 12x12x4 value minus
QE's own printed D3 term (-912.74971137 + 0.22170600 = -912.528005 Ry) to all six digits, and
band 44 sits 0.021 meV from QE's.

| $l$ | 3 | 9 | 15 | 21 | 27 | 33 | 39 | 45 |
|---|---|---|---|---|---|---|---|---|
| QE (meV) | 59.1712 | 38.8172 | 3.1618 | 8.4211 | 1.4860 | 2.8381 | 0.9864 | 1.1468 |
| defumat (meV) | 59.1296 | 38.7834 | 3.1613 | 8.4180 | 1.4833 | 2.8395 | 0.9843 | 1.1480 |
| deviation | -0.07% | -0.09% | -0.02% | -0.04% | -0.19% | +0.05% | -0.21% | +0.10% |
| Wannier (meV) | 59.4812 | 39.0640 | 3.1750 | 8.4434 | 1.4946 | 2.8382 | 0.9964 | 1.1423 |

Worst 0.21%, median 0.08%. **The noise floors coincide to 0.02%: 27.49 vs 27.48 μeV** --
tighter than the 0.04% seen at the coarse mesh, and the strongest form of the result, since a
floor reproduced by an implementation sharing no linear algebra, no FFT, no ultrasoft
augmentation and no Hubbard code cannot be an artefact of either one.

**P7 is closed on both densities.** Two independent statements come out of it: the three routes
agree on the coefficients, and two of them independently agree on the *error bar*.

**Still open:** the derivative-based check (a Taylor expansion in $\mathbf k$ against the
Fourier coefficients from the sampled circle), and re-running the three-route comparison at
`degauss = 0.015`, since the coefficients above are all on the under-converged smearing.

## P11 — Is the differentiable route real? — **yes, and it buys something other than speed**

**Establishes:** whether `jax.grad` genuinely reaches through defumat's band structure, since the
entire case for developing on defumat rather than on QE rests on it. If it does not, defumat is
only a slower second opinion and the QE route (P9) is the whole method.

**Independent check:** the derivative against a finite difference of the band structure it is
supposed to be the derivative of -- two routes through the same code that share no arithmetic.

**Measured: it works.** defumat already ships `defumat.response.velocity`, which builds
`v = dH/dk` by a `jvp` at a frozen plane-wave sphere, and `band_velocities(calculation, result,
kpoints=, nbnd=)` returns `d(eps_n)/dk` at arbitrary k-points. Against central differences at
three random general-position k-points, all three Cartesian directions, on a 1-atom As cell:

| bands (well separated) | near-degenerate |
|---|---|
| **1e-6** relative | up to 1e-3 relative |

The degenerate case is expected rather than a defect: a band velocity is a diagonal expectation
value, which the eigensolver's arbitrary mixing rotates inside a degenerate manifold, and the
module says so and offers `matrix_elements` (the whole block) for a degeneracy-safe consumer.

It also satisfies every hard requirement this project has: arbitrary k-points (so the circle can
be fed to it directly), the ultrasoft `dS/dk` term (the generalised Hellmann-Feynman derivative
`<psi|dH/dk - eps dS/dk|psi>`, from the same `jvp`), Hubbard `ns` forwarded, and collinear spin.

**Cost: derivatives are ~1.5x the eigenvalues alone**, measured at nk = 24 with 16 bands.

**And that is what changes the conclusion.** Values-plus-derivatives at N points carry the
information of ~2N samples, so the effective Nyquist in $l$ doubles for 1.5x rather than 2x.
That is a real saving and a small one -- and it is *not* the argument for defumat, because QE's
bands step is ~20x faster per k-point to begin with (6 min vs 125 min for the same 152-point
circle). Brute-forcing k-points in QE remains about 12x cheaper than the derivative route in
defumat for the same effective resolution. **The derivative route does not win on the noise
floor.**

**Parameter derivatives are NOT the easy case, and that was measured rather than assumed.**
The hope was that `jax.grad` would give $\partial a_l/\partial U$ or $\partial a_l/\partial R$
in one pass, since those are the two largest uncertainties here and both are invisible to the
forbidden-harmonic floor. It does not, and the reason is structural rather than a software gap.

A frozen-density derivative is correct for $k$ because moving $k$ does not perturb the
self-consistent potential. Moving an atom, or changing $U$, does: the density responds, and the
eigenvalue moves again. Total-energy forces escape this by stationarity -- which is exactly why
`defumat.forces.autodiff` differentiates `frozen_energy` and is right to -- but **an eigenvalue
is not a stationary quantity and has no such protection.**

Measured on a 2-atom As cell, $d\varepsilon_n/dx_2$ at a general k-point, comparing a finite
difference with full SCF re-convergence (the true derivative) against one at frozen density
(what a naive autodiff gives):

| | median | worst |
|---|---|---|
| self-consistent response term, as a fraction of the true derivative | **66%** | **285%** |

Not a correction: the dominant part, and on band 7 it flips the answer's magnitude entirely.
Any parameter derivative of an eigenvalue therefore needs the full linear response --
**Sternheimer / DFPT**, not `jax.grad`.

**Which leaves a narrower but still real case for defumat.** It *has* a Sternheimer solver
(`defumat.response.sternheimer`, transcribing QE's `LR_Modules/cgsolve_all.f90`,
`ch_psi_all.f90`, `orthogonalize.f90`, `h_prec.f90`, `setup_alpha_pv.f90`), and its response
chain is autodiff-generated where QE hand-codes a routine apiece -- its own docstring tabulates
`incdrhoscf`/`addusdbec`/`lr_addusddens` becoming one `jvp` of the density w.r.t. the states,
`newdq`'s `int3` becoming one `jvp` of `newd`. So **adding a new perturbation type** (a Hubbard
$U$ derivative, a strain derivative aimed at inverse design) is genuinely cheaper there than in
PHonon. That is the honest argument: not that autodiff makes parameter derivatives free, but
that it makes a *new* DFPT perturbation cheaper to implement.

**And the pragmatic path for the immediate question is neither.** The $U$ sensitivity does not
need DFPT at all -- three or four QE SCF runs at different $U$ answer it in a few hours, at the
production mesh, with no new code. That should come first, and it decides whether $\partial
a_l/\partial U$ is worth building machinery for.

Scripts in `examples/differentiable/`.

## P12 — The Hubbard $U$: the largest uncertainty in the pipeline

**Establishes:** what $U = 3.0$ eV costs. It is an inherited choice, nothing here had tested
it, and the forbidden-harmonic floor cannot see it -- $U$ is an input to the Hamiltonian and the
floor measures the fit. Third and worst member of the family P3's smearing result opened.

**Independent check:** four full SCF + 152-point bands runs at $U = 2, 3, 4, 5$ eV, identical in
every other respect (12x12x4, `degauss` 0.0150), scored against the hand-derived allowed set.

**Measured:** band 44, amplitudes in meV.

| $l$ | $U$=2 | $U$=3 | $U$=4 | $U$=5 | 2 -> 5 |
|---|---|---|---|---|---|
| 3 | 69.152 | 63.746 | 54.660 | 44.260 | -36% |
| 9 | 45.094 | 42.921 | 36.080 | 28.892 | -36% |
| 15 | **0.000** | 3.533 | 5.358 | 5.940 | switched on |
| 21 | 9.550 | 8.644 | 6.889 | 5.141 | -46% |
| 27 | 1.098 | 1.837 | 1.984 | 1.762 | +61% |
| 33 | 2.962 | 2.614 | 2.069 | 1.457 | -51% |
| 39 | 1.244 | 1.168 | 0.889 | 0.602 | -52% |
| 45 | 1.235 | 1.003 | 0.761 | 0.482 | -61% |
| **noise floor** | 81.8 μeV | 24.3 μeV | 6.5 μeV | 1.3 μeV | /65 |

**21% per eV at the median, 24% at worst** (excluding $l=15$; it alone moves 52% from $U$=3 to
4, which is about its zero crossing rather than a typical sensitivity). Against 12.9% for the
smearing error, under 4% for a 1e-3 A geometry perturbation, 2.4% of $l=45$ for the floor:
**$U$ dominates every other uncertainty this project has measured.**

**It changes the shape, not the scale.** Six harmonics shrink monotonically with $U$ while
$l=27$ rises then falls and $l=15$ grows from nothing, so no rescaling maps one $U$'s spectrum
onto another's. $l=15$ is the sharpest case: at $U=2$ it is 0.3 μeV, **250x below that run's own
noise floor** -- the harmonic is absent, not small -- and its coefficient is positive at
$U$=3,4,5 approaching zero from above at 2, so it crosses zero near $U \approx 2$. A selection
rule says which harmonics *may* be non-zero; which ones *are* is a question about the
interaction strength. Every allowed harmonic stays pure $\sin(l\theta)$ at every $U$, which is
the control that keeps this a physics result rather than a fitting artefact.

**The floor is a function of $U$ too**, 81.8 -> 1.3 μeV. **A quoted noise floor is a property
of (settings, $U$), not of the material.**

**What the floor tracks -- a correction that plain PBE forced.** On $U = 2..5$ alone the floor
falls monotonically as band 44 moves away from $E_F$, and this phase first read that as P3's
Fermi-surface-sampling mechanism seen from the other side. A $U = 0$ run breaks it: band 44 sits
at +0.0332 eV, *between* $U$=2's +0.0166 and $U$=3's +0.0357, yet the floor is **0.92 μeV**, the
smallest of all five and 90x below $U$=2's.

| $U$ | band 44 $-E_F$ (eV) | floor (μeV) | residual, the $l>45$ tail (meV) |
|---|---|---|---|
| 0 | +0.0332 | **0.92** | **0.053** |
| 2 | +0.0166 | 81.76 | 1.146 |
| 3 | +0.0357 | 24.30 | 0.795 |
| 4 | +0.0526 | 6.51 | 0.466 |
| 5 | +0.0656 | 1.26 | 0.232 |

What it does track is the **fit residual**, i.e. the signal living above $l_{max}$ -- exactly
what P3 identified as the thing that aliases in. Five-point correlation **+0.91 against the
residual**, -0.77 against the gap (the gap's apparent -0.94 was an artefact of the $U \geq 2$
subset). The corrected statement is the stronger one: it confirms P3's aliasing diagnosis
mechanistically rather than resting on a coincidence.

**And plain PBE is not a neutral baseline, it is a different physics.** $l=9$ drops 2.8x,
$l=21$ drops 4.5x, and $l=15$ comes out with the **opposite sign** (-5.540 meV against +3.533 at
$U$=3), confirming the zero crossing near $U \approx 2$ from below. $l=3$ is not the $U=2$ trend
continued either -- 64.2 meV at $U$=0 against 69.2 at $U$=2, so 0 to 2 is not monotonic.
Removing $U$ therefore answers nothing: it substitutes one unargued choice for another, and for
Cr 3d it is the physically worse one.

**Consequence, and the next step.** Every coefficient in this project is a coefficient *at
$U = 3.0$ eV* and must be reported as such. The principled fix is not a scan but a
first-principles $U$: QE's `hp.x` computes it by linear response (DFPT on the Hubbard manifold),
and it is now built in this tree. Pinned in `tests/test_p12_hubbard_u.py`.

## P12b — Which conclusions survive the choice of $U$

Standard practice is to pick a reasonable $U$ and check that it does not change the conclusion.
P12's scan is what turns the second half of that sentence into a measurement. Separating the
two groups explicitly, over $U = 2..5$ eV, band 44:

**Robust**

| claim | evidence |
|---|---|
| the allowed set is the odd multiples of 3 | pure symmetry, no DFT input |
| every allowed harmonic is pure $\sin(l\theta)$ | worst cos/sin ratio 5e-4; cos amplitude below that run's own floor at every $U$ |
| the leading two are $l=3$ then $l=9$ | identical ordering at all four $U$ |
| $l=45$ is measurable | 15x (at $U$=2) to 383x (at $U$=5) the floor |
| $l_{max}=45$ suffices and the fit is well posed | residual 0.6-1.9% of signal, condition number 1.4, at every $U$ |

**Not robust**

| claim | why |
|---|---|
| any individual coefficient | $l=3$ spans 43% over the range, ~21% per eV |
| the full ordering | $l=15$ runs from last ($U$=2) to third ($U$=5) |
| $l=15$ existing at all | 0.3 μeV at $U$=2, **below** that run's floor |
| "every allowed harmonic is above the floor" | true at $U$ = 3,4,5, false at $U$ = 2 |

**Not tested, and claimed for neither group:** whether the QE/Wannier/defumat agreement holds
away from $U = 3$. All three were only ever compared there, and the Wannier model is a fit to a
DFT run at that $U$, so closing this means refitting the wannierisation.

So the defensible sentence is: *"We take $U = 3.0$ eV and have verified that the conclusions of
this work are robust for $U \in [2,5]$ eV; the absolute coefficients vary by about 21% per eV
and are reported at $U = 3.0$."* `hp.x` would upgrade that to a first-principles $U$, but it is
**not a prerequisite** -- the methodological results, the three-route validation, the noise
floor and the existence of high harmonics all sit in the robust column.

**A process note that belongs in the record.** The first version of this separation was written
in prose from a script whose pure-sine check had an empty loop body: it printed the conclusion
without computing it. The conclusion happened to be correct. Three neighbouring claims in the
same prose were not -- the floor called "μeV-level" is 1.3-82 μeV, "every allowed harmonic is
measurable" fails at $U$=2, and the three-route agreement was listed as robust having never
been measured away from $U$=3. Pinned, including the weakened forms, in
`tests/test_p12_robustness.py`.

## P13 — The pseudopotential and the functional: the last two untested inputs

**Establishes:** the size of the two remaining input sensitivities the noise floor cannot see.
P12 closed the Hubbard $U$; these were the entries `README.md` listed as "not measured".

**Independent check:** one knob at a time against the corrected reference (12x12x4, `degauss`
0.0150, $U$ = 3.0, band 44), with the band's identity verified by its distance above $E_F$
rather than by its index -- the check P8's CoF2 result showed is not optional.

**Measured**, and it completes the blind-spot table:

| input | median | worst | floor sees it? |
|---|---|---|---|
| Hubbard $U$, per eV | **21%** | -- | no |
| **PAW vs ultrasoft** | **9.7%** | 17.8% ($l$=45) | no |
| **PBEsol vs PBE** | **7.3%** | 33.2% ($l$=39) | no |
| smearing set too narrow | 12.9% | -- | no |
| geometry at 1e-3 A | <4% | -- | partly |
| *the floor's own number* | *0.04%* | -- | -- |

**Both eat the high harmonics specifically.** $l=3$ moves -1.3% (PBEsol) and +0.4% (PAW), while
$l=39$ moves -33.2% and -12.1% and $l=45$ moves -21.7% and -17.8%. Neither is a rescaling: the
leading harmonic is essentially untouched and the tail -- the part this project exists to
measure -- takes the damage.

**And the floor goes *down* while the answer goes wrong:** 24.30 -> 8.15 (PBEsol) and 14.28 μeV
(PAW), because it tracks the residual above $l_{max}$ (0.795 -> 0.495, 0.647 meV) and not the
error. A run whose coefficients are 20% off reports a **better** noise floor than the reference.
That is P12's mechanism correction at its sharpest, and the cleanest single demonstration that
the floor is not an error bar on a coefficient.

**A general constraint found on the way.** The first attempt used SG15 ONCV norm-conserving
pseudopotentials, a bigger contrast with ultrasoft than PAW is. It cannot be done: their
`<PP_PSWFC>` block is **empty**, and DFT+U needs atomic wavefunctions to build the Hubbard
projectors, so `pw.x` stops with `determine_hubbard_occ (1): no atomic wavefunctions in
pseudopotential file`. **A Hubbard calculation restricts which pseudopotential families exist
for it at all** -- a real limitation on this kind of test, not a misconfiguration.

A second trap on the same path: pslibrary's default Se PAW (`-dn`) carries 3d in valence, 16
electrons rather than 6, which takes the cell from 90 electrons to 130 and stops band 44 being
the same band. The `-n` variant restores 6; all three runs are verified at 90 electrons and 54
KS states. Pinned in `tests/test_p13_pseudo_and_functional.py`.

## P8 — Production: MnTe and CoF2, and where the high harmonics actually are

**Establishes:** the method on materials with independent literature classifications, rather
than only on the group's own Cr2SeAs where the answer comes from this project.

**Independent check:** two systems chosen so their symmetry requires **opposite** answers on the
same plane. MnTe's g-wave form $\Delta E \propto k_z k_\perp^3\cos 3\theta$ vanishes identically at
$k_z=0$; CoF2's rutile d-wave $\Delta E \propto k_xk_y$ carries no $k_z$ at all and must not. One
projector, one pipeline, two required answers that cannot both be satisfied by a systematic error.

**Measured, and both hold.** Band nearest $E_F$, 152-point circle of radius $0.25|b_1|$:

| | rms splitting at $k_z=0$ | required |
|---|---|---|
| CoF2 (U=0) | **86.6 meV** | non-zero |
| MnTe (U=0) | **19.25 μeV** | nodal |
| MnTe (U=4) | **0.20 μeV** | nodal |

A ratio of 4500 at $U=0$ and $4\times10^5$ at $U=4$, from symmetry alone.

**And the function is the right one of the two.** $l=2$ has two members, $k_xk_y$ and
$k_x^2-k_y^2$. Rutile must give the first: F sits at $(u,u,0)$ so the $x=y$ diagonal is a mirror
through the Co-F bonds that maps each Co to itself ($\chi=+1$), under which $k_xk_y$ is even and
$k_x^2-k_y^2$ odd. The symmetrizer predicted it before any DFT (keeps 100% of $\sin2\theta$,
0.00% of $\cos2\theta$); the DFT gives cos/sin = **9.5e-15**.

**Neither material has a symprec cliff** -- P6_3/mmc (12+12 ops) and P4_2/mnm (8+8) stable from
1e-6 to 1e-2, where Cr2SeAs flips at 1e-5. So that trap is not universal, but it still has to be
scanned per structure to know.

**MnTe resolves the full ladder to $l=45$.** Band 20, $U=4$: 949.38 / 53.45 / 12.18 / 3.96 /
1.466 / 0.582 / 0.243 / 0.105 meV for $l = 3..45$, the weakest at **4493x the noise floor**.
Three separate reasons, each confirming an earlier phase from the other side:

1. a **semiconductor**, so P3's Fermi-surface sampling problem does not arise -- floor
   0.02-0.03 μeV against Cr2SeAs's 23.6 μeV;
2. **exact Wyckoff positions**, so no relaxation noise -- cos/sin at $l=3$ is 1e-13 to 4e-12
   against 8.6e-10 for relaxed Cr2SeAs at the same $l$;
3. a splitting an order of magnitude larger.

Band 20's ladder is **exponential**, $A \propto e^{-l/6.35}$, $R^2 = 0.997$ over four decades,
beating a power law -- and it passes a check that does not go through the fit: the residual left
above $l_{max}=45$ is 0.0374 meV where the same exponential summed over $l=51..75$ gives 0.0425.
**Band 21 is not** a clean exponential (step ratios 5.4, 10.1, 3.1, 1.25), so this is a band
property and one band is not a law.

**The headline, and it did not go as expected.** Second allowed harmonic over the leading one,
on the band nearest $E_F$:

| | leading | ratio |
|---|---|---|
| CoF2 | $l=2$, 46.7 meV | **0.021%** |
| MnTe | $l=3$, 949.4 meV | **5.6%** |
| **Cr2SeAs** | $l=3$, 64.3 meV | **67.6%** |

**The two canonical altermagnets are essentially single-wave; the group's own material is the
multi-harmonic one** -- its $l=9$ is two thirds of its $l=3$. That is where "the DFT high
harmonic part would be genuinely new" actually lives, and now with two contrasting references.
The comparison survives a caveat pointing the other way: CoF2's allowed set steps by 4 rather
than 6, and a shorter step should decay *less*.

**Two failure modes this phase produced, both recorded as tests.**

*One overwritten circle.* Both MnTe circles are bands runs over one density, and the first job
wrote both into `tmp/MT.save/` with the fixture copied afterwards -- so the $k_z=0$ run silently
replaced the $k_z\ne0$ eigenvalues and the committed fixture was the wrong circle under the
right name. Three tests failed on it, which is the only reason it surfaced. The CoF2 job had it
right: copy out between runs.

*Band indices do not survive a change of $U$.* At $U=0$ CoF2's band 30 is the valence maximum
at -0.256 eV with a 540 meV gap to band 31; at $U=4$ both sit within 21 meV of $E_F$ and are
**4.07 meV apart**. Comparing by index across that shows up as $l=2$ apparently collapsing from
46.7 to 6.9 meV with the floor jumping 150x. Not a failed calculation -- both converged, moments
agree (5.65 vs 5.72 $\mu_B$), compensation holds, Hubbard occupations symmetric. This is
`CLAUDE.md`'s band-identity rule in a form the project had only guarded across *models with
different band counts*, not across *values of a parameter in one code*. What does survive is the
symmetry: $l=2$ stays pure $k_xk_y$ to 1e-11 on every band at every $U$, because a selection
rule is a statement about the crystal and not about a band.

Pinned in `tests/test_p8_mnte.py` (14) and `tests/test_p8_cof2.py` (10).

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

**Re-verified on the corrected settings.** After P3's smearing result moved every coefficient,
`harmwave.x` was re-run on the 24x24x8 / degauss 0.0150 reference (the wavefunctions had to be
regenerated first -- see the trap below). It reproduces the Python route to **4.8e-15 eV** at
all eight allowed harmonics, and independently computes the same noise floor, **23.61 μeV**.
The Fortran tool therefore tracks the corrected reference with no change to its own code.

**A trap in this project's own disk hygiene.** `harmwave.x` reads a finished run through QE's
`read_file()`, which needs the **wavefunctions**, not just the charge density and XML. Deleting
`tmp/AB.wfc*` to save quota makes it fail with `Error in routine read_wfc (2)`. The density is
enough to redo the bands step and regenerate them (~6 min), but wavefunctions are not a
disposable intermediate the way they are for the Python route.

**v2: the plane is read from the k-points, not assumed.** The $k_x k_y$ restriction was
hardcoded in two places -- a check that every k-point shared one $k_z$, and a filter testing
whether row 3 of $R$ was $(0,0,1)$. v2 fits the plane to the k-list (normal from the first three
non-collinear points, every point then verified), measures $\theta$ in an in-plane frame
$(e_1,e_2)$ written into the report, and generalises the filter to its coordinate-free form:
$k' = Rk$ preserves $n\cdot k$ for every $k$ exactly when $R^{T} n = n$. The plane's *offset*
needs no separate test, since $R^{T}n = n$ gives $n\cdot(Rk) = n\cdot k$ automatically.

For $n = +\hat z$ this reduces to v1's row-3 test and the frame becomes $(x,y)$ exactly, so a
$k_x k_y$ circle must reproduce v1 bit for bit. On the corrected Cr2SeAs reference it does:
**all eight allowed harmonics identical** ($l=9$ differs by 2.3e-13, the rest by 0), largest
absolute difference over all 91 coefficients **1e-14 eV**, normal detected as $(0,0,1)$, offset
0.0550574 = $a/4c$ as it should be, and exactly four new lines in the report. Pinned in
`tests/test_p9_arbitrary_plane.py`.

**What that does not show:** v2 reproduces v1 *where v1 was valid*. No circle with a non-$+\hat
z$ normal has been run, so v2's correctness for a general plane is untested. Closing it needs a
bands run on e.g. a $k_y k_z$ circle with the allowed set compared against
`altharm.symmetry`'s prediction for that plane -- the same cross-check that caught v1's
time-reversal bug.

**Scope of v2, honestly still not general:** collinear LSDA only; band identity is given
explicitly rather than tracked (`altharm.bandtrack`'s continuity tracking is not yet ported to
Fortran); symmetrization and the noise-floor split stay in the already-tested Python code rather
than being duplicated in Fortran. See `qe_native/README.md`.

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
- **Structure precision.** *Measured.* The relaxation that produced this structure converges
  forces, not coordinates, and leaves atoms good to ~1e-3 A; the sublattice-exchanging mirror
  $M_z$ holds to **1.6e-5 A**. Two orders of magnitude separate them. Random Cartesian
  displacements of per-component RMS `amp`, fixed seed, rigid translation removed:

  | amp (A) | worst allowed-$l$ shift | noise floor | $l=0$ |
  |---|---|---|---|
  | 0 | -- | 24.3 μeV | 10.7 μeV |
  | 1e-4 | 0.36% | 29.4 μeV | 76.3 μeV |
  | **1e-3** (relaxation) | **3.70%** | **166.4 μeV** | **676.0 μeV** |
  | 1e-2 | 34.8% | 1695.3 μeV | 6708.4 μeV |

  **Relaxation precision is good enough for the coefficients and not for the error bar**: at
  1e-3 A the strong harmonics move under 0.3% while the floor degrades 7x. Since the floor is
  what says how much to trust a coefficient, such a structure yields numbers one can mostly
  believe and no honest way to say so. Recovering the intrinsic ~24 μeV needs ~1e-4 A, so the
  structure must be **symmetrised before fitting**, not taken raw from the relaxation.
  $l=0$ is the sharpest single indicator -- the angular average must be exactly zero for a
  compensated altermagnet, and it reads 676 μeV at 1e-3 A.

  On the symmetry side the collapse is a cliff at **1.99e-5 A**: the allowed set goes from 8
  straight to all 49, never losing two or three. The symprec that survives an amplitude `amp`
  is ~5.4x `amp`, so a 1e-3 A structure needs symprec >= 5.4e-3, i.e. **500x spglib's default**.

  **The forbidden spectrum identifies which generator broke**, which is more than an error bar:
  z-displacements break $M_z$ (even multiples of 3 leak in, worst $l=12$, and $l=0$ rises to
  536 μeV), in-plane displacements break $C_3$ (non-multiples leak in, worst $l=2$, while
  $l=0$ stays 17x smaller because the sublattice compensation survives). This reproduces the
  spglib classification exactly (z-noise -> P3m1, 6 preserving, 0 exchanging; in-plane -> P1).
  Corollary: in-plane noise leaves every coefficient within 0.18% of baseline while multiplying
  the floor by 6.7 -- **coefficient stability is not evidence that a structure is symmetric
  enough**. Pinned in `tests/test_geometry_precision.py` (15 tests).

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

---

## P10 — Why a short-range Wannier model still carries l = 45 — **resolved**

This answers the question that motivated the DFT route in the first place: *"there were some
terms that are not captured by the Wannier"*. The worry is exactly right for a linear
response, and exactly wrong for an eigenvalue. Both halves are now measured.

**The linear argument.** On a circle of radius $k_0$,
$e^{i\mathbf{k}\cdot\mathbf{R}} = \sum_l i^l J_l(k_0 R_\parallel) e^{il(\theta-\phi_R)}$, so any
quantity *linear* in $H(\mathbf{R})$ has its $l$-th harmonic weighted by $J_l(k_0 R)$, which is
exponentially small until $k_0 R \sim l$. For Cr2SeAs $k_0 = 0.4868\ \text{Å}^{-1}$:

| $l$ | 3 | 9 | 15 | 21 | 27 | 33 | 39 | 45 |
|---|---|---|---|---|---|---|---|---|
| range needed $l/k_0$ (Å) | 6.2 | 18.5 | 30.8 | 43.1 | 55.5 | 67.8 | 80.1 | 92.4 |
| in unit cells | 1.7 | 5.0 | 8.3 | 11.6 | 14.9 | 18.2 | 21.5 | 24.8 |

The group's model reaches 25.8 Å in plane. By this argument it should produce **nothing** above
$l\approx 13$.

**The control** (`tests/test_p10_nonlinearity.py`). On one and the same Hamiltonian truncated to
8 Å, three quantities, two of them manifestly linear in $H(\mathbf{R})$:

| $l$ | matrix element $H_{22}$ (linear) | trace $\mathrm{Tr}\,H$ (linear, gauge-invariant) | eigenvalue $\Delta E$ (non-linear) |
|---|---|---|---|
| 3 | 2.1e-01 | 9.6e+00 | **59.05** |
| 9 | 8.1e-05 | 3.6e-03 | **42.71** |
| 15 | 3.7e-10 | 1.6e-08 | **2.43** |
| 21 | 2.0e-13 | 1.3e-11 | **8.66** |
| 33 | 5.0e-13 | 7.1e-12 | **2.54** |
| 45 | 1.7e-13 | 2.6e-12 | **0.87** |

(meV.) The linear channels die on the Bessel schedule and hit the numerical floor by $l=21$.
The eigenvalue splitting does not: at $l=45$ it is **ten orders of magnitude** larger. With the
full 583-vector model the linear cutoff simply moves out as predicted ($l=27$ at 3.6e-9,
$l\ge33$ at the floor) while $\Delta E$ is unchanged.

**Conclusion.** The high harmonics are *not* hopping-range physics. They are generated by the
**non-linearity of the eigenvalue problem** — level repulsion at the six avoided crossings
(P4), which is a sharp feature in $\theta$ and therefore has high harmonics however short-ranged
$H(\mathbf{R})$ is. Twenty-nine hoppings inside 6 Å already give 14% of the converged $l=45$.

**But truncation is not free either.** $l=45$ runs 0.155 → 0.869 → 0.947 → 1.110 → 1.143 meV
over cutoffs 6 / 8 / 11.5 / 15 / ∞ Å. So a Wannier model needs ~15 Å for the weakest harmonics
to be within 5% — a factor of ~7, not the 40 orders of magnitude the linear picture predicts.

**What this means for the two routes.** It explains P6 (DFT ↔ Wannier agreeing to <1% at all
eight harmonics) rather than leaving it a coincidence, and it reframes the DFT route's value:
DFT is *not* needed to recover harmonics a Wannier model structurally cannot represent — it can
represent them. DFT is needed because it removes the fit's *quantitative* truncation error and
because it needs no Wannierisation at all, which matters for screening new materials.
