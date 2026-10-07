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
orthonormal: the module docstring's claim holds for the spherical basis only.

**Closed 2026-10-07, by the second option.** `basis.py`'s module docstring now states that the
circular basis is "orthogonal on a uniform circle but **not normalised** -- the Gram matrix over
N sample points is diag(N, N/2, N/2, ...), not the identity", and warns that anything comparing a
circular coefficient with a spherical one must not assume otherwise. `CircularBasis.evaluate`
was *not* rescaled, which keeps the coefficients in the units people expect (a_l in eV, directly
the amplitude of cos(l θ)).

The projection step is safe by construction rather than by convention: `project` uses
`np.linalg.lstsq`, which is invariant under any scaling of the basis. **Verified** by scaling
the whole basis by 7 and refitting -- the reconstructed function agrees to 7.1e-15 while the
coefficients differ by exactly 7.000. So no part of the pipeline assumes unit norm.

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

On the group's real Cr2Se2As the allowed set at $k_z \ne 0$ came out
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

**Measured, on the group's real Cr2Se2As circle** (56 k-points, $k_z = 0.25\,b_3$,
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
all*. Scanning 50-100 Ry in 2 Ry steps on the real Cr2Se2As circle, reporting the $l=6$ amplitude
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
Wannier90 route (`hr.dat` interpolation, already validated to 0.07 meV against DFT on Cr2Se2As)
agree on the allowed-$l$ coefficients, each within the other's noise floor. Second — the part
that's actually new — exactly how a finite hopping range $R_\mathrm{max}$ damps a chosen
high-$l$ harmonic, turning the Bessel-function suppression argument in `README.md` into a
measured curve rather than an assertion.

**Independent check:** the QE/Wannier agreement is itself two codes sharing no implementation.
The truncation study's own check: build a *synthetic* `hr.dat` containing a known harmonic at a
controlled radius, truncate it at a series of $R_\mathrm{max}$ values, and confirm the measured
damping matches the $J_l(k_0 R_\mathrm{max})$ prediction rather than some other falloff.

**Measured: the truncation half is closed; the QE half is not started.**
`altharm/wannier.py` reads `hr.dat` and the study has been run on the group's v03 Cr2Se2As
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
group's real Cr2Se2As run, a 152-point circle (P3) puts the DFT-fitted allowed harmonics within
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

**Cost, measured rather than guessed.** On the real Cr2Se2As cell, defumat's SCF costs
$37.8\,\mathrm{s} + 6.8\,\mathrm{s} \times N_k$ per iteration (fit to 1, 2 and 4 k-points;
the 80% fixed part is density, FFT, ultrasoft augmentation and Hubbard, none of which scales
with $k$). For QE's 69 irreducible k-points that is ~8.4 minutes an iteration, ~5 hours for the
35 iterations QE needed — feasible overnight. The band step on the 152-point circle turned out
to be the larger cost and is being timed separately.

**Measured, and P7 is closed for the harmonic comparison.** defumat ran the full pipeline on
the real Cr2Se2As cell: SCF to convergence (45 iterations, 43.8 min) then the 152-point circle
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

**The three routes, all on Cr2Se2As band 44:** QE (Fortran plane-wave), Wannier90 (tight-binding
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
than only on the group's own Cr2Se2As where the answer comes from this project.

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
1e-6 to 1e-2, where Cr2Se2As flips at 1e-5. So that trap is not universal, but it still has to be
scanned per structure to know.

**MnTe resolves the full ladder to $l=45$.** Band 20, $U=4$: 949.38 / 53.45 / 12.18 / 3.96 /
1.466 / 0.582 / 0.243 / 0.105 meV for $l = 3..45$, the weakest at **4493x the noise floor**.
Three separate reasons, each confirming an earlier phase from the other side:

1. a **semiconductor**, so P3's Fermi-surface sampling problem does not arise -- floor
   0.02-0.03 μeV against Cr2Se2As's 23.6 μeV;
2. **exact Wyckoff positions**, so no relaxation noise -- cos/sin at $l=3$ is 1e-13 to 4e-12
   against 8.6e-10 for relaxed Cr2Se2As at the same $l$;
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
| **Cr2Se2As** | $l=3$, 64.3 meV | **67.6%** |

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
same real Cr2Se2As data, are compared coefficient by coefficient.

**Measured:** built and validated once, on the group's real Cr2Se2As run (PBE+U, USPP, 56-point
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
the given circle. Cr2Se2As's sublattice-exchanging mirror flips $k_z$ on its own and so fails
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
$k_x k_y$ circle must reproduce v1 bit for bit. On the corrected Cr2Se2As reference it does:
**all eight allowed harmonics identical** ($l=9$ differs by 2.3e-13, the rest by 0), largest
absolute difference over all 91 coefficients **1e-14 eV**, normal detected as $(0,0,1)$, offset
0.0550574 = $a/4c$ as it should be, and exactly four new lines in the report. Pinned in
`tests/test_p9_arbitrary_plane.py`.

**What that did not show, and how it was closed.** v2 reproduces v1 *where v1 was valid*, which
left its correctness on a general plane untested. **That gap is now closed** (job 20383535, run
2026-09-22, analysed 2026-10-01): a great circle in the $k_y k_z$ plane, normal $(1,0,0)$,
offset 0, frame $(e_1,e_2) = (\hat y, \hat z)$, 4 of 12 operations kept. Its allowed set is
$l = 2,4,\dots,24$, matching `altharm.symmetry`'s independent prediction for that plane, and
the coefficients are sharper still -- even $l$, sine only, with odd $l$ at 1e-14 eV. See P14.1,
which also records what that run revealed about the noise floor itself.

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
  the group's own Cr2Se2As structure, spglib reports the non-magnetic parent as **P3m1 (#156)
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

- **The circle's plane and the basis frame must be the same plane.** *Measured, and it was
  silent.* `CircularBasis` takes `e1`/`e2` and measures the angle there; handing it directions
  from another plane returns angles, a design matrix and an allowed set, all plausible and all
  wrong. On a 152-point Cr2Se2As circle in the $k_yk_z$ plane:

  | basis frame | $d\cdot n$ spread | off-block max in $S$ | allowed $l$ |
  |---|---|---|---|
  | default $(x,y)$ -- wrong | 1.96 | **165** | {3,6,8,9,10,11,12,13,14,15,16,19,21} |
  | $(y,z)$ -- correct | 5.6e-17 | 1.4e-07 | **{2,4,6,...,24}** |

  The correct frame gives the even harmonics, which is what the four surviving operations
  require (one is a two-fold axis sending $\theta \to \theta+\pi$). The wrong one gives thirteen
  scattered values with no group-theoretic reading, because the symmetrizer stops being
  block-diagonal in $l$ and `allowed_harmonics`'s per-block diagonalisation loses its premise.
  `angles()` now refuses directions that are not on a cone about its own normal; found while
  setting up P9's non-$+\hat z$ verification, and pinned in
  `tests/test_circular_basis_frame.py`.

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
exponentially small until $k_0 R \sim l$. For Cr2Se2As $k_0 = 0.4868\ \text{Å}^{-1}$:

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

---

## P14 — the sample count: the floor is anti-correlated with the error

**Status: mechanism established two independent ways; the DFT confirmation is prepared and
has not been run.**

P3/P12/P13 established five inputs the forbidden-harmonic floor is blind to. Blindness is
survivable. This phase found a sixth where it is worse than blind.

### P14.1 The y-z circle (job 20383535, run 2026-09-22, read 2026-10-01)

harmwave.x v2 on a great circle in the $k_y k_z$ plane, same density as `sm24`. It reads the
plane from the k-points correctly: normal $(1,0,0)$, offset 0, frame $(e_1,e_2)=(\hat y,\hat z)$,
4 of 12 operations kept. The allowed set matches `altharm.symmetry`'s independent prediction,
$l = 2,4,\dots,24$ — and the coefficients are sharper than that prediction: **even $l$, sine
only**, odd $l$ at 1e-14 eV, even-$l$ cosine at 5.7e-7 eV, even-$l$ sine up to 20.4 meV. Pure
sine means $\Delta E$ is odd in $k_z$; consistently, $\Delta E(k_z{=}0) = -3.4$ μeV.

This closes P9's open item ("no circle with a non-+ẑ normal has been run").

### P14.2 The floor has channels, and they differ by $10^7$

The exact zeros prompted a check that harmwave.x does not symmetrise before fitting — it does
not (`design` is raw cos/sin columns, then `DGELSS`, harmwave.f90:285-295). So the zeros are a
property of the data. Asking the data directly how well it obeys $\Delta E(R\mathbf k) =
\chi\,\Delta E(\mathbf k)$ for each operation that permutes the sample points:

| circle | $\chi=+1$ | $\chi=-1$ |
|---|---|---|
| $k_yk_z$ | 2.8e-14 eV | 1.09 μeV (+ a separate 4.74 μeV rigid $l=0$ offset) |
| $k_xk_y$ | 0.072 μeV | 1.60 μeV |

$\chi=+1$ relates $\varepsilon_\uparrow(\mathbf k)$ to $\varepsilon_\uparrow(R\mathbf k)$ — same
spin block, exact by construction, so what is left is k-point placement roundoff. $\chi=-1$
relates $\varepsilon_\uparrow$ to $\varepsilon_\downarrow$ and holds only as well as the
self-consistent density respects sublattice exchange.

### P14.3 Where the published 23.6 μeV floor actually comes from

No operation that permutes the samples has a residual above 1.6 μeV, so the floor cannot come
from one. Splitting all 90 coefficients of the $k_xk_y$ run by what forbids each:

| group | n | median | max |
|---|---|---|---|
| allowed | 8 | 3093 μeV | 64293 μeV |
| forbidden by an operation the grid **can** represent | 67 | **0.0040 μeV** | 1.11 μeV |
| forbidden **only** by C3, which it **cannot** | 15 | **24.70 μeV** | 75.0 μeV |

The published floor is the third row. 120° is not a multiple of 360°/152, so C3 is the one
operation this grid cannot represent, and the floor measures the aliasing of that single
constraint. This gives P3's empirical result (9 meV → 1.01 μeV over 56 → 152 → 304 points) a
mechanism rather than a correlation.

### P14.4 The rule, and the inversion

A material's splitting contains only allowed harmonics at every order; on Cr2Se2As the allowed
set is $l \equiv 3 \pmod 6$, a lattice of period $p=6$. The part above $l_\mathrm{max}$ aliases
to $l' = |l - mN|$, which is allowed exactly when $mN \equiv 0 \pmod p$. So the tail first
reaches the signal at alias order

$$m = p / \gcd(N, p).$$

Measured on a synthetic signal whose truth is known exactly (only allowed harmonics, no noise,
tail fitted to the real coefficients), sweeping every $N$ from 120 to 199:

| $\gcd(N,6)$ | $m$ | reported floor | true error (median / max) |
|---|---|---|---|
| 6 | 1 | **0 μeV** | **18.1% / 36.2%** |
| 3 | 2 | 160 μeV | 3.7% / 6.3% |
| 2 | 3 | 221 μeV | 1.4% / 2.6% |
| 1 | 6 | 239 μeV | **0.29% / 1.04%** |

Monotone, no exceptions in that range, and stable across five tail decay exponents. When $p \mid
N$ the strongest alias of every unresolved harmonic lands directly on a reported coefficient and
**no forbidden channel receives anything at all**: the floor reads zero and the answer is ~20%
wrong. **The floor is monotonically anti-correlated with the error it is supposed to report.**

Practical rule: **take $N$ coprime to the period of the allowed set; never a multiple of it.**
This project's own 152 is $\gcd=2$ — not the catastrophic case, but ~5x worse than a coprime
count at identical cost, and chosen by accident.

### P14.5 The y-z run is itself the catastrophic case

Its allowed set is every even $l$, so $p=2$, and $N=152$ is even: alias order 1. The signature
is exactly as predicted — it has no C3-style channel at all, and every forbidden amplitude sits
at the machine/SCF level. Its own coefficients therefore carry the tail as error, and its
forbidden channels cannot report it.

### P14.6 The mechanism, measured on real eigenvalues at two sample counts

The synthetic tests control the truth but have no SCF noise; the classification is real data at
one sample count. The 304-point circle closes the gap between them with no new DFT run, by
being used twice: fitted to its own Nyquist ($l=151$) it measures the high-$l$ tail directly,
and subsampled to every second point it *is* a 152-point circle on identical eigenvalues.

On a 152-point grid $\sin(l\theta_j)$ with $l = 152-l'$ folds to $-\sin(l'\theta_j)$, so the
prediction is signed: $b_{l'} = -b_{152-l'}$.

| $l'$ | from | $b_{l'}$ at N=152 | $-b_{152-l'}$ at N=304 | residual |
|---|---|---|---|---|
| 5 | 147 | −33.755 μeV | −32.617 μeV | −3.49% |
| 11 | 141 | +40.339 | +39.324 | +2.58% |
| 17 | 135 | −44.636 | −43.724 | −2.08% |
| 23 | 129 | +54.666 | +53.835 | +1.54% |
| 29 | 123 | −60.772 | −60.059 | −1.19% |
| 35 | 117 | +74.249 | +73.591 | +0.89% |
| 41 | 111 | −86.259 | −85.661 | −0.70% |

**All seven signs correct**, magnitudes within 0.7–3.5%, and the residual shrinks monotonically
with $l'$ — it is what the second fold partner ($l = 152+l'$, above 304's Nyquist) would supply,
implying $b_{157} \approx -1.14$ μeV and $b_{193} \approx -0.60$ μeV.

**So the forbidden-harmonic "noise floor" is not noise.** It is the signal's own unresolved
tail, folded down, and it can be predicted from an independent measurement of that tail to a few
percent. Control: the eight allowed coefficients agree between the two counts to <0.02%, as
required since both are $\gcd=2$ and neither is corrupted at first order.

### What is NOT done

P14.6 measures the *mechanism* on real data, but not the *catastrophic case*: no subsample of a
304-point grid is divisible by 3, so the $m=1$ regime — floor reading zero while the
coefficients absorb the tail — remains shown only synthetically. The DFT confirmation is
prepared in `altharm_runs/nsweep` (N=150/151/153 on $k_xk_y$ and N=151 on $k_yk_z$, all reusing
the `sm24` density, each with its own outdir) and **has not been submitted**. Prediction: N=150
reports a floor near 0.004 μeV and shifts the coefficients by order 10%; YZ151 shifts relative
to the existing 152-point y-z result. Until then P14.3 is an observation about one run and
P14.4 is a statement about a model.

Code: `altharm/sampling.py`. Tests: `tests/test_p14_sampling_aliasing.py` (14).
Notebook: `notebooks/P14_sampling_aliasing.ipynb` (executed).

### P14.7 The rule across three materials, and what bounds the damage

The allowed-set period is a property of the magnetic point group, not of this project's
conventions, so the rule has to be checked where that period differs. On the existing P8
fixtures:

| material | N | p | gcd | m | floor | next/leading |
|---|---|---|---|---|---|---|
| MnTe U=4 | 152 | 6 | 2 | 3 | 0.0234 μeV | 5.630% |
| **CoF2 U=0** | 152 | 4 | **4** | **1** | 0.0369 μeV | **0.021%** |
| Cr2Se2As U=3 | 304 | 6 | 2 | 3 | 1.0082 μeV | **65.601%** |

**CoF2's allowed set is $l = 2,6,10,14$ — period 4 — and 152 is divisible by 4, so the
published P8 CoF2 runs have been in the $m=1$ regime all along.** That went unnoticed.

It did them no harm, and the reason is the refinement that matters: **the damage at $m=1$ is
proportional to the unresolved tail**, and CoF2 barely has one. Verified synthetically by
scanning the tail amplitude at fixed $N=150$ — the error tracks the tail and vanishes with it.

Two consequences:

1. **A low floor does not diagnose $m=1$.** MnTe is at $m=3$ and its floor (0.0234 μeV) is
   *lower* than CoF2's at $m=1$ (0.0369 μeV). Both are low because both materials are nearly
   single-wave, not because of their sampling.
2. **The failure mode is not uniform across materials — it bites hardest on the interesting
   ones.** A steep harmonic ladder is nearly immune to any $N$; a slow one is at risk. Cr2Se2As,
   whose second harmonic is 66% of its first, is precisely the kind of material this project
   exists to study, and precisely the kind $m=1$ would destroy.

**Corollary for P8, now measured** (job 20813368, 2026-10-07). CoF2 rerun at $N = 151$, which is
coprime to 4 and so sits at $m = 4$, the safest sampling available at that period -- same
density, same circle, only the count changed:

| $l$ | $N$=152 ($m$=1) | $N$=151 ($m$=4) | difference |
|---|---|---|---|
| 2 (leading) | 46.67131 meV | 46.67130 meV | **+0.000%** |
| 6 | 0.00999 | 0.00999 | +0.034% |
| 10 | 0.00469 | 0.00468 | +0.241% |
| ~~14~~ | ~~0.00000~~ | ~~0.00001~~ | ~~−78%~~ |

**$l=14$ must not be quoted**: the two values are 1.9 neV and 8.8 neV against a floor of
35.3 neV -- 18× and 4× *below* the floor, so that 78% compares two numbers neither run can
resolve. Excluding
it, the meaningful coefficients agree to **0.034% median, 0.241% worst**.

Against Cr2Se2As at the same contrast (2.18% median, 10.00% worst, P14.8), CoF2 is roughly
**40× less affected**. So "the damage is proportional to the unresolved tail" holds on real
data, from the protected side as well as the exposed one: **P8's CoF2 numbers stand, and the
reason they stand is that CoF2 has no tail, not that the sample count was chosen well.**

---

## P15 — the k_z = 0 residual is a rigid offset, not an angular error

**Open since P8.** MnTe's $k_z = 0$ plane is nodal: symmetry forbids *every* $l$ there, so the
whole measured splitting is residual. It was 19 μeV at $U=0$ and 0.20 μeV at $U=4$, and the note
recorded only that this is "not a fixed floor". What it actually is, measured 2026-10-06 on the
existing `mnte_kz0` / `mnte_kz0_u4` fixtures (band 20, 152-point circle, allowed set empty in
both, as `altharm.symmetry` confirms):

| $U$ | total rms | $l=0$ rigid part | angular part | gap (band 20→21) |
|---|---|---|---|---|
| 0 | 19.252 μeV | **19.245 μeV** | 0.500 μeV | 759 meV |
| 4 | 0.197 μeV | **0.197 μeV** | 0.009 μeV | 1813 meV |

**99.96% of it is a rigid, angle-independent shift between the two spin channels.** The angular
content -- the part that would actually corrupt a harmonic -- is 0.5 μeV at $U=0$ and 9 neV at
$U=4$. So the number that looked alarming was never a harmonic error at all: a constant in
$\theta$ is $l = 0$, and $l = 0$ is not among the coefficients anyone reports.

Decomposed by P14's channels, it sits **entirely in $\chi = -1$**:

| $U$ | $\chi=+1$ residual | $\chi=-1$ residual | ratio |
|---|---|---|---|
| 0 | 1.5e-11 eV | 1.0e-06 eV | 65000 |
| 4 | 5.8e-12 eV | 1.7e-08 eV | 3000 |

(The $\chi=+1$ figures are small but *not* machine zero -- quoting them as 0 would overstate
the case by five orders.)

which is the same split P14.2 found on Cr2Se2As: constraints relating $\varepsilon_\uparrow$ to
itself are exact, constraints relating it to $\varepsilon_\downarrow$ hold only as well as the
self-consistent density respects sublattice exchange. That is why it is $U$-dependent and not a
floor -- opening the gap by 2.4x improves it by 98x.

**Consequence.** A nodal plane is a *good* diagnostic of SCF quality, but its residual must be
reported with $l=0$ removed, or it overstates the angular error by a factor of 40 (U=0) to 22
(U=4). The same applies to the $l=0$ offset P14.2 found on the y-z circle (4.74 μeV of a
4.86 μeV total).

Tests: `tests/test_p15_kz0_residual.py`.

---

## P16 — what the floor does as N grows: neither $1/N_k$ nor a generic noise plateau

**Open since P3** ("whether it saturates at a genuine noise level is not yet measured -- it
would take a third k-density"). The third density was never needed: P14 says the floor *is* the
unresolved tail folded down, so it can be predicted from an independent measurement of that
tail. The 56-, 152- and 304-point circles share one density (60 Ry, degauss 0.0037), and the
304-point fit resolves $l$ up to 151, which is the tail the other two are folding.

Reading the floor only from the channels aliasing actually feeds -- $l'$ forbidden, $N-l'$
allowed -- as P14.3 requires:

| $N$ | channels | measured floor | predicted from the 304-point tail | ratio |
|---|---|---|---|---|
| 56 | 3 | 1050.4 μeV | 986.2 μeV | **1.065** |
| 152 | 7 | 54.67 μeV | 53.84 μeV | **1.015** |

Two independent confirmations across a 19x range in the floor, on real eigenvalues, with no
fitted parameter.

**So the floor does not follow a power of $N$ at all.** It follows the *tail's own decay*, which
is a property of the material's harmonic ladder -- fast for a nearly single-wave altermagnet,
slow for a multi-harmonic one. Quoting a $1/N_k$ law would be fitting the wrong thing.

**Where it stops.** The tail cannot be followed below the point where the $\chi = -1$
constraints stop being exact, because those are what make a forbidden harmonic forbidden in the
first place. On this density that error is **0.766 μeV** (and $\chi = +1$ is 0.072 μeV). The
304-point floor is 1.865 μeV, i.e. **2.4x above it** -- so Cr2Se2As at 304 points is approaching
the crossover but has not reached it, and a further increase in $N$ can buy at most about a
factor of two before the floor stops responding.

**The saturation value is therefore not a universal numerical noise level.** It is the
self-consistency quality of that particular density for that particular material, which is why
P15 finds the same quantity controlling MnTe's nodal-plane residual and why opening a gap
improves it.

Tests: `tests/test_p16_floor_scaling.py`.

---

## P17 — the Bessel suppression, on a synthetic `hr.dat` where the answer is known

**Open since P10**, which asked for it in so many words: "build a *synthetic* `hr.dat`
containing a known harmonic at a controlled radius, truncate it at a series of
$R_\mathrm{max}$ values, and confirm the measured damping matches the $J_l(k_0 R_\mathrm{max})$
prediction rather than some other falloff."

The argument is Jacobi-Anger. A hopping at in-plane displacement $R$, angle $\phi$, contributes
$e^{i k_0 R\cos(\theta-\phi)}$ to $H(\mathbf k)$ on a circle of radius $k_0$, whose $l$-th
Fourier coefficient is $i^l J_l(k_0 R)e^{-il\phi}$. So for a **linear** quantity -- a matrix
element, a trace -- the $l$-th harmonic from a shell at radius $R$ is exactly proportional to
$J_l(k_0 R)$.

**Measured, single shell** (square lattice $a = 3$ Å, one orbital, hoppings on one shell with
angular pattern $\cos l\phi$, $k_0 = 0.9$ Å$^{-1}$): the ratio of the measured harmonic to
$J_l(k_0 R)$ is constant across ten shells to **4e-15**, for each of $l = 2, 4, 6, 8$. That is
not a falloff resembling a Bessel function; it *is* the Bessel function, zeros and sign changes
included.

**A lattice artifact had to be excluded first, and it nearly read as a 5% failure of the law.**
A shell whose radius is a Pythagorean hypotenuse is not one angular orbit: $|R| = 5$ holds 12
sites in 3 distinct orbits ($(5,0)$, $(3,4)$, $(4,3)$), so $\cos l\phi$ on it does not excite a
single $l$. Those shells deviate; every other shell is exact.

**Measured, truncation.** For a model with many shells the harmonic is a sum over everything
retained, so it does **not** track $J_l(k_0 R_\mathrm{max})$ as a smooth envelope -- the ratio
wanders over a factor of 70. What survives is the **threshold**: $J_l(x)$ is negligible for
$x \ll l$, so a model truncated below $R \sim l/k_0$ cannot represent harmonic $l$ at all. At
half that radius the harmonic is suppressed by $10^3$ to $10^5$ relative to the converged value,
and the suppression is of the *order* of $J_l(k_0 R_\mathrm{max})$, within about a decade.

So the honest form of the claim is a threshold, not an envelope -- consistent with P10's own
note that the falloff "is not monotonic, it overshoots before settling", and with why a 4-7 Å
model gets $l = 21$ wrong by a factor of six while a 15 Å one does not.

Tests: `tests/test_p17_bessel_synthetic.py`.

---

## P18 — the Taylor-vs-Fourier check cannot be done at the circle's centre, and why

P7 proposed a second, internal cross-check on the fit: "a Taylor expansion in $\mathbf k$
against the Fourier coefficients recovered from a sampled circle, at a point where both should
agree." It was never run, and the reason it cannot be run as stated is worth recording.

The check's content does not need defumat. If $\Delta E(\mathbf k)$ were analytic at the
circle's centre, the lowest-degree polynomial in $(k_x, k_y)$ producing $\cos l\theta$ would be
degree $l$, so

$$a_l \propto k_0^{\,l} \quad\text{as } k_0 \to 0,$$

a parameter-free prediction testable at any radius with the Wannier model (`v03`, 38 Wannier
functions, 583 $R$-vectors). **Measured 2026-10-06**, circles of radius 0.003–0.034 Å$^{-1}$ at
$k_z = 0.25 b_3$:

| band | $l$ | fitted exponent | Taylor predicts | $R^2$ |
|---|---|---|---|---|
| 20 | 3 | 1.42 | 3 | 0.995 |
| 21 | 3 | **1.33** | 3 | 0.996 |
| 22 | 3 | 1.11 | 3 | 0.993 |
| 21 | 9 | −0.22 | 9 | **0.27** |

$l=3$ scales roughly **linearly**, not cubically; $l=9$ does not follow a power law at all.

**The reason is a degeneracy at the centre, and it is a symmetry fact rather than an accident.**
At $k_x = k_y = 0$ the point sits on the $\Gamma$--A line, whose little group has two-dimensional
irreps, so bands come in degenerate pairs: measured gaps of **0.005, 0.009 and 0.010 meV** for
bands 19/20, 21/22 and 24/25. $\Delta E$ is therefore non-analytic there -- a square-root
singularity, not a Taylor series -- and an exponent near 1 is that singularity's fingerprint.

**This corroborates P10 from the opposite direction.** P10 concluded that the high harmonics are
generated by the non-linearity of the eigenvalue problem (level repulsion at avoided crossings)
rather than by hopping range. P18 shows the same non-analyticity defeats the Taylor expansion:
the two statements are the same fact seen from either side. A Taylor-based check would have to
be centred somewhere the bands are non-degenerate, which is not where this project's circles
are, and where it would no longer test the harmonics of interest.

Tests: `tests/test_p18_taylor_vs_fourier.py` (skipped when the v03 `hr.dat` is absent; it lives
outside the repo at `HIGH_THROUGHT_WANNIERcr2X3/1.cr2se2as/Cr2SeAs_wannier_v03/`).

---

## P14.8 — the $m=1$ regime, measured on real DFT

P14 established the rule synthetically and measured the *mechanism* on real eigenvalues
(P14.6), but the catastrophic case itself -- floor reading zero while the coefficients absorb
the tail -- had only a model behind it. Job 20794763: the same `sm24` density (24×24×8,
degauss 0.0150, $U = 3$), **only the sample count changed**.

### $k_x k_y$ plane, allowed set $l \equiv 3 \pmod 6$, $p = 6$

| $N$ | $\gcd$ | $m$ | reported floor | coefficient shift vs $N=152$ (median / max) |
|---|---|---|---|---|
| **150** | 6 | **1** | **0.0032 μeV** | **2.18% / 10.00%** |
| 153 | 3 | 2 | 32.47 μeV | 0.06% / 0.22% |
| **151** | 1 | **6** | **45.66 μeV** | **0.01% / 0.02%** |
| 152 | 2 | 3 | 24.70 μeV | (reference) |

**Anti-correlated and monotone, three points, on real data.** As $m$ rises the error falls
(2.18 → 0.06 → 0.01%) while the reported floor rises (0.003 → 32.5 → 45.7 μeV).

$N = 153$ is what makes this an ordering rather than two extremes: it lands squarely between on
both axes.

**$N = 150$ has no aliasing-fed forbidden channel at all.** $150 \equiv 0 \pmod 6$, so
"$l'$ forbidden and $150 - l'$ allowed" has no solution -- the tail lands *entirely* on allowed
harmonics. The floor there had to be read from the unfed channels. That is the definition of
$m = 1$, not a statistical accident.

**Control:** $N = 151$ and $N = 152$ are two different counts, neither at $m=1$, and they agree
to **0.02%**. So the $N=150$ shift is caused by $m=1$, not by changing $N$.

### $k_y k_z$ plane, allowed set every even $l$, $p = 2$ -- a second period

The existing 152-point y-z circle is even, so $\gcd(152,2) = 2$ and $m = 1$; 151 is odd, giving
$m = 2$, the maximum available at $p = 2$.

| | forbidden odd-$l$ channels | $l = 8$ coefficient |
|---|---|---|
| 152 points ($m=1$) | **0.0000 μeV** | 0.0529 meV |
| 151 points ($m=2$) | **60.1 μeV** | **0.0930 meV** |

The 152-point run reported a **perfect zero** in every odd channel while its weakest harmonic
was **43% off** -- nearly a factor of two. Median shift over the twelve allowed $l$: 1.44%.

### The honest correction to P14.4

The synthetic model predicted a median error of **18%** at $m=1$; the measurement is **2.18%**
on the $k_xk_y$ circle and 1.44% on the y-z one. Mechanism and direction confirmed, **magnitude
an order of magnitude smaller**, because Cr2Se2As's real tail above $l = 45$ decays faster than
the synthetic tail assumed. Quote the measured numbers, not the model's.

### What this cost, and the two guards it bought

The first attempt (20794209) failed twice over. `harmwave.x` was launched bare inside the
allocation and hit an Open MPI PMI error -- the module stack changed after 2026-10-02, so a
script that worked in September no longer does. And the input generator took everything *before*
`K_POINTS`, but `sm24/bands.in` carries its HUBBARD card at line 207, **after** the k-point
list, so all four tasks silently ran plain PBE. `grep -ci hubbard` gave 0 against sm24's 8, and
nothing errored.

The tell was that all three $N$ agreed with each other to 1e-6 while all disagreeing with the
reference: three different sample counts cannot do that unless the difference is systematic
rather than an $N$ effect. The job script now checks the HUBBARD card is present, that Hubbard
is active in `bands.out`, and launches every MPI binary through `srun`.

Tests: `tests/test_p14_sampling_aliasing.py` (the nsweep block).

---

## P19 — the linear-response $U$ is twice the one in production

**Open since P12**, which measured what $U = 3.0$ eV costs but could not say whether 3.0 is the
right number: it is an inherited choice and nothing in this project had tested it. The
linear-response route (QE's `hp.x`, DFPT) answers that directly. Job 20796619, finished
2026-10-06.

| | |
|---|---|
| **linear-response $U$** | **6.3858 eV** (CrU) / **6.3880 eV** (CrD) |
| **$U$ used in every production run** | **3.0 eV** |
| ratio | **2.13x**, a difference of 3.39 eV |

The two magnetic sublattices agree to 0.002 eV, which symmetry requires and is the check that
the run is self-consistent rather than drifting.

**Verified before believing it:** the Hubbard projectors are `ortho-atomic`, matching production;
the SCF underneath ran at $U = 0.0001$ eV, which is the correct linear-response procedure (the
response is taken about $U \approx 0$), not a typo.

### What it implies, and what it does not

P12's scan is monotone -- every coefficient *falls* as $U$ rises:

| $l$ | $U$=2 | $U$=3 | $U$=5 | 2 → 5 |
|---|---|---|---|---|
| 3 | 69.152 | 63.746 | 44.260 | −36% |
| 45 | 1.235 | 1.003 | 0.482 | −61% |

So the direction is unambiguous: **the production coefficients are systematically too large.**

**But 6.39 eV lies outside the 2–5 eV range P12 scanned, so the magnitude cannot be stated.**
Extrapolating P12's "21% per eV" across 3.39 eV would be arithmetic, not measurement. Closing
that needs one 152-point bands run at $U \approx 6.4$ -- the same cost as any other single bands
run, about 40 minutes.

### Two reservations that must be quoted with the number

1. **The q-mesh is 2×2×1 (2 q-points).** `hp.x` results depend on it, and this is coarse.
2. **`conv_thr_chi` was loosened to 1e-4.** At 1e-6 the response would not converge -- the
   residue oscillated at 1e-4 with roughly 230 further iterations needed (recorded in P12's
   notes at the time).

So 6.39 eV is not a number to put in a paper. It is a strong signal that **the production $U$ is
badly low**, and since P12 established that $U$ dominates every other uncertainty this project
has measured, that makes it the largest open issue in the pipeline.

### Why it took four attempts, and what the failure was

Three of the four response matrices existed since 2026-09-22; the collection step died in 25
seconds and was not revisited. The error was `AB.chi.pert_1.dat does not exist`, and the cause
was a **missing stage**, not a broken one. QE's HP workflow is three steps:

```
dns.pert_i.q_j  --(sum_pertq, summing over q per atom)-->  chi.pert_i  --(compute_hp)-->  U
```

The runs had produced the `dns` files and then jumped straight to `compute_hp`, skipping
`sum_pertq` entirely. Adding that stage took 40 minutes of walltime, all of it queueing --
summing files already on disk is close to free. **Several days of DFPT was sitting one cheap
step short of an answer for two weeks.**

Fixture: `tests/data/hp_u/AB.Hubbard_parameters.dat`. Tests: `tests/test_p19_linear_response_u.py`.

---

## P19b — what the linear-response $U$ does to the numbers, measured

P19 left the magnitude open: 6.39 eV lies outside the 2–5 eV range P12 scanned, so the shift
could not be stated without extrapolating. Job 20796793 closes that — identical to P12's
`u2`/`u4`/`u5` in every respect (12×12×4, degauss 0.0150, `ecutwfc` 60, 152-point circle,
band 44), with $U = 6.3858$ (CrU) / $6.3880$ (CrD) from P19.

| $l$ | $U$=2 | $U$=3 | $U$=4 | $U$=5 | **$U$=6.39** | 3 → 6.39 |
|---|---|---|---|---|---|---|
| 3 | 69.152 | 63.746 | 54.660 | 44.260 | **31.403** | **−50.7%** |
| 9 | 45.094 | 42.921 | 36.080 | 28.892 | **20.671** | −51.8% |
| 15 | 0.000 | 3.533 | 5.358 | 5.940 | 5.526 | +56.4% |
| 21 | 9.550 | 8.644 | 6.889 | 5.141 | 3.132 | −63.8% |
| 27 | 1.098 | 1.837 | 1.984 | 1.762 | 1.203 | −34.5% |
| 33 | 2.962 | 2.614 | 2.069 | 1.457 | 0.697 | −73.3% |
| 39 | 1.244 | 1.168 | 0.889 | 0.602 | 0.258 | −77.9% |
| 45 | 1.235 | 1.003 | 0.761 | 0.482 | **0.156** | **−84.5%** |
| floor | 81.8 μeV | 24.3 | 6.5 | 1.3 | **1.30 μeV** | |

(meV unless marked. $l=15$ is the exception P12 already flagged: it sits near a zero crossing,
so its percentage is about that crossing rather than a sensitivity.)

**Median shift from the production $U$ to the linear-response one: −57.8%.** The leading
harmonic is halved; the weakest is down by a factor of six.

### The conclusion that matters: the ratios survive, the absolute values do not

| $U$ | next/leading ($l$=9 / $l$=3) | $l=3$ (meV) |
|---|---|---|
| 2 | 65.2% | 69.152 |
| 3 | 67.3% | 63.746 |
| 4 | 66.0% | 54.660 |
| 5 | 65.3% | 44.260 |
| **6.39** | **65.8%** | **31.403** |

**The next/leading ratio moves 1.5 percentage points across a range where the leading harmonic
changes by 55%.** So this project's headline result -- that the high-harmonic physics lives in
Cr2Se2As (67.6%) and not in the canonical altermagnets (MnTe 5.6%, CoF2 0.021%), P8 -- **is
untouched by the $U$ uncertainty.** All eight allowed harmonics also stay resolved over the
whole range.

What does *not* survive is any statement quoting an absolute amplitude in meV. Those depend on
$U$ at the factor-of-two level, and the right value of $U$ is not settled.

### Where this leaves the project

$U$ was already the largest measured uncertainty (P12: 21% per eV). P19b makes it larger than
that framing suggested, because the question is no longer sensitivity but which value is
correct. Two readings are consistent with everything measured:

1. **The linear-response $U$ is right**, and every absolute coefficient published from this
   project is roughly a factor of two too large.
2. **$U=3$ is defensible on other grounds** (agreement with experiment elsewhere, the coarse
   2×2×1 q-mesh and loosened `conv_thr_chi` behind the 6.39), and the discrepancy is a known
   limitation of DFPT $U$ for this system.

Nothing here decides between them -- that needs a comparison against something external, which
this project does not have. **What it does establish is that the choice does not threaten the
conclusions, only the numbers.**

### The decision: $U = 3.0$ eV (2026-10-07)

**Taken by the project owner, not derived here.** All reference values stay at $U = 3.0$ eV and
the published coefficients stand as they are.

This is a choice, and the record should not be read as saying 3.0 was validated. What is
established is:

- The linear-response value is **6.39 eV**, a factor of 2.13 away, with the reservations in P19
  (2×2×1 q-mesh, `conv_thr_chi` loosened to 1e-4).
- Had 6.39 been adopted instead, every absolute coefficient would move by a **median of −57.8%**
  (P19b) -- the leading harmonic halves, the weakest drops sixfold.
- **The conclusions are unaffected either way.** The next/leading ratio varies 1.5 percentage
  points over $U$ = 2 … 6.39, so P8's headline -- the high-harmonic physics is in Cr2Se2As, not
  in the canonical altermagnets -- holds at any $U$ in that range, as does the resolution of all
  eight allowed harmonics.

So the honest form of any statement built on these numbers is: *at $U = 3.0$ eV*, the
coefficients are X; the ratios are robust to $U$; a linear-response treatment would give
substantially smaller absolute values. P19 and P19b exist so that sentence can be written
without hand-waving.

Fixture: `tests/data/qe_cr2seas_u64/`. Tests: `tests/test_p19b_u_extended.py`.

---

## P20 — the 3D spherical route, on a real material for the first time

**Open since P2.** `SphericalBasis` and the spherical projector had unit tests and nothing else:
no DFT run had ever put a sphere of k-points through them. Job 20814237 does — 288 points
(12 × 24 Gauss-Legendre × uniform) on a sphere of radius 0.293879 tpiba, chosen to pass through
the existing 152-point circle so the two routes can be compared directly. Same `sm24` density.

### What it settles: the allowed set, and the $l=3$ / $l=4$ question

The symmetry projector's prediction, made without touching the data: **allowed $l = 4, 6, 8, 10$**.
Measured, band 44:

| $l$ | allowed? | largest $|$coefficient$|$ |
|---|---|---|
| 1, 3, 5, 7, 9 | no | **0.00000 meV** |
| 2 | no | 0.00090 meV |
| 4 | yes | 34.17 meV |
| 6 | yes | 28.53 meV |
| 8 | yes | 98.42 meV |
| 10 | yes | 23.62 meV |

Signal-to-forbidden ratio **110 000**. The first allowed $l$ is **4**.

**This resolves a real and recurring confusion.** "g-wave means $l=4$" and this project's own
"MnTe is g-wave with leading order $l=3$" are both correct, in different bases. Expanding the
g-wave form $\Delta E \propto k_z k_\perp^3 \cos 3\phi$ in real spherical harmonics gives
**exactly one term, $Y_4^3$** -- the other 48 coefficients below $10^{-6}$, relative residual
$6.2\times10^{-17}$. So:

- **spherical $l = 4$** -- which is the "g" of g-wave (s/p/d/f/g = 0/1/2/3/4), and what the
  nodal-surface count is really counting ($Y_4^3$ has three vertical nodes from $\cos3\phi$ plus
  one horizontal from $\cos\theta$);
- **circular $l = 3$** -- the azimuthal index left after $k_z$ is fixed.

Reading a circular result against "g-wave = $l$=4" judges a correct answer wrong. The 3D route
makes the correspondence explicit, which is something the 2D route cannot do. **That, rather
than being one more cross-check, is what this route is for.**

### What it does NOT settle: the amplitudes

The fit leaves a **10.2% relative residual** at $l_\mathrm{max}=10$, and $l=8$ comes out larger
than $l=4$, which is the signature of aliasing rather than of physics. Predicting the measured
circle from the spherical fit confirms it:

| $l_\mathrm{max}$ | spherical residual | predicted circle vs measured |
|---|---|---|
| 4 | 22.2% | 111.6% |
| 6 | 21.6% | 98.5% |
| 8 | 11.4% | 60.9% |
| 10 | 10.2% | **43.1%** |

Monotone but nowhere near converged. The cause is the grid, not the method: 12 × 24
Gauss-Legendre is exact only for $l < 12$, while Cr2Se2As carries a resolved ladder out to
$l = 45$ on the circle. Everything above $l=10$ folds into the coefficients below it -- P14's
aliasing, in three dimensions.

**So this run establishes the allowed set and the $l=3$ ↔ $l=4$ correspondence, and its
amplitudes must not be quoted.**

**What closing that would take, quantified.** Matching the circle's $l_\mathrm{max}=45$ needs
$n_\theta > 45$ and $n_\phi > 90$, i.e. about 46 × 92 = 4232 points. At this run's rate (288
points in 80 minutes) that is roughly **20 hours** -- inside a two-day walltime, but 14× the
cost, and worth doing only if the 3D amplitudes are wanted for their own sake rather than for
the symmetry statement, which is already settled.

Fixture: `tests/data/sphere/`. Tests: `tests/test_p20_spherical_route.py`.
