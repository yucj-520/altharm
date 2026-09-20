# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this
repository.

**This file is the rules. `PLAN.md` is the record.** How to work here — conventions,
refusals, what a finished phase looks like — is below. What was actually found — per-phase
validation numbers, the traps each phase turned up, what a noise floor was measured at — goes
in `PLAN.md`, one section per phase. Do not restate a phase's findings here; add them there.

## How to talk to the user

**The user is a physicist, so write to them as one, in plain English.** Say what a quantity
*is*, physically, before saying what the code does with it. Keep sentences short. When a
choice needs their input, put it in terms a physicist who writes some Python will recognise —
what a knob changes in the calculation, what it costs, what it gets wrong — not in
software-engineering vocabulary that assumes the question is about the code rather than the
physics.

This binds replies, questions and summaries in the session. It does **not** loosen `PLAN.md`,
docstrings or comments: the record stays exact and numeric, because the whole rule of this
project is stated below and is not optional.

## The rule everything else follows from

**A claim is a number, not an adjective.** "The noise floor is low" is not a finding; "the
forbidden-$l$ amplitude is 0.3 meV against a 4 meV allowed-$l$ coefficient" is. This applies
to code comments, docstrings, commit messages and conversation alike — if a sentence in this
project could be true or false depending on taste, it isn't finished yet.

Two consequences that are easy to skip past:

- **A test must check against a route independent of the implementation under test.** Calling
  the function and asserting it returns what it returns is not a test. Calling the function
  twice with the same input and asserting the two calls agree is not a test either — it
  catches nondeterminism and nothing else. What counts: a synthetic $\Delta E(\theta)$ built
  from $a_l, b_l$ chosen by hand and checked recovered by the fit (`PLAN.md` P2); a symmetry
  selection rule derived by hand from a character table and checked against the numerical
  projector's output (P1); two DFT codes, or DFT against Wannier, or either against defumat,
  agreeing on a coefficient neither one's implementation shares with the other (P6, P7); a
  known closed-form integral (an orthonormality relation, P0). If a proposed test shares code
  with the thing it's testing — the same projection matrix, the same symmetrizer, the same
  angle convention — it is not independent, however different the test file looks.
- **A missing term is refused by name, not silently approximated.** If `symmetry.spin_operations`
  is asked to close a set of operations with $\mathbf k \to -\mathbf k$ under an assumption
  that only holds without spin-orbit coupling, and the caller's calculation has spin-orbit
  coupling, that is a wrong answer produced with no diagnostic — it must raise, saying which
  assumption failed, rather than return a plausible-looking set of allowed $l$. Same for an
  unimplemented route: a `--method=wannier` flag before the Wannier reader exists raises
  `NotImplementedError("Wannier route not implemented; see PLAN.md P6")`, it doesn't fall back
  to whatever route happens to be written.

## What this project is

`altharm` fits the angular harmonics of an altermagnet's spin splitting $\Delta E(\mathbf k) =
\varepsilon_\uparrow(\mathbf k) - \varepsilon_\downarrow(\mathbf k)$ directly from plane-wave
DFT, on a circle (2D) or a sphere (3D) in k-space, and cross-checks the result against the
group's existing Wannier route and against defumat. `README.md` has the physics and the
organising idea — **fit every $l$ up to $l_\mathrm{max}$, not only the symmetry-allowed ones,
because the forbidden ones measure the calculation's own noise floor.** Read it before
touching the fitting or symmetry code; nothing here restates it.

**Status**, honestly: `altharm/units.py`, `basis.py`, `symmetry.py`, `kpath.py`,
`project.py`, `noise.py`, `bandtrack.py`, `planewaves.py`, `wannier.py` and `qe/` (writer +
reader) are written and tested (99 tests). The QE route has been run for real, twice: bcc iron
(interface validation) and the group's Cr2SeAs PBE+U structure (production). On Cr2SeAs, P3's
noise floor is measured at 0.106 μeV RMS over the forbidden harmonics — about five orders of
magnitude below the ~9 meV the group measured re-converging self-consistently on a k-mesh —
and P4 has a physical account of where the residual concentrates (the plane-wave count's own
angular spectrum peaks at the same $l$ the noise floor does). A second, independent
implementation exists too: `qe_native/harmwave.f90`, a new compiled QE executable (linked
against `PW/src`/`PP/src`, using QE's own `find_sym()` and LAPACK directly, no XML, no
spglib, no numpy) that agrees with the Python route to $10^{-9}$ eV on the same real data — see
`qe_native/README.md`. The Wannier90 route is written and cross-validated against the group's
real Cr2SeAs model. The defumat route is not started. `PLAN.md`'s phase list (P0-P9) is the
map of what's open and what each open thing is waiting on — P5 (ecutwfc/k-density convergence)
and P7 (defumat) are the largest gaps.

## Layout

- `altharm/units.py` — physical constants and conversions. **Single source of truth**: a
  conversion factor defined a second time somewhere else is a bug waiting to disagree with
  this one.
- `altharm/basis.py` — the angular bases (`CircularBasis` for a circle, `SphericalBasis` for
  real spherical harmonics on a sphere) that everything else is expressed in.
- `altharm/symmetry.py` — the numerical projector: which magnetic-space-group operations
  survive a given magnetic structure, tagged by whether they preserve or exchange the
  sublattices, and which angular harmonics that leaves allowed.
- `altharm/qe/` — the Quantum ESPRESSO interface: `writer.py` builds a `pw.x` input for a
  `bands` run on an explicit k-circle, `reader.py` parses the resulting `data-file-schema.xml`.
- `altharm/bandtrack.py` — following one band's identity around a k-path by continuity,
  for the "band index vs. band identity" convention below.
- `altharm/planewaves.py` — an independent plane-wave-count sphere counter, used to size the
  basis-set-discontinuity noise source (P4) and to cross-check `qe.reader` indirectly.
- `qe_native/` — `harmwave.f90`, a second, independent implementation of the harmonic fit and
  the magnetic-symmetry classification, compiled directly into Quantum ESPRESSO as a new
  `PP/src` executable rather than written in Python. See its own `README.md`.
- `examples/` — worked scripts against real materials (currently Cr2SeAs).
- `tests/data/` — committed real QE/Wannier output the tests check against; see the
  independence rule above for why these are real runs, not synthetic fixtures, wherever a
  synthetic one would not do.
- `PLAN.md` — the record: one numbered phase per question, what settles it, and the measured
  number once there is one.
- A git repository (`git log` is meaningful history, not an empty tree).

## Non-negotiable conventions

- **Directions are Cartesian, always.** `altharm.basis` only accepts Cartesian unit vectors,
  never fractional/crystal coordinates, because a fractional angle silently carries the
  lattice metric — on a hexagonal cell this turns a circle into an ellipse and mixes harmonics
  by $\sqrt3$. Anything that reads k-points out of a QE, Wannier90 or defumat output and hands
  them to a basis converts to Cartesian at that boundary, not inside the basis.
- **A band is identified by its character, not by its sorted position.** Eigenvalues come back
  sorted by energy; at a crossing the $n$-th band swaps identity and $\Delta E_n(\theta)$ jumps.
  Any code that follows one band around a circle or a sphere has to say how it tracks identity
  through a near-degeneracy, and that method is a first-class thing to test, not an assumed
  detail.
- **Self-consistency on the k-mesh the harmonics are fit on is a noise source, not a detail.**
  A frozen-density (non-self-consistent) run is the default assumption for the QE route; if a
  self-consistent alternative is ever used, say so next to the result, because the group has
  already measured a ≈9 meV noise floor from doing that on a related calculation, against a
  ≈1 meV splitting.
- **No conversion factor is a bare literal.** `27.211386245988` belongs in `units.py` once;
  everywhere else imports it by name.
- **Symmetry-forbidden coefficients are never dropped, only labelled.** The whole point of
  fitting all $l$ up to $l_\mathrm{max}$ is that the forbidden ones are the noise floor; code
  that only fits the allowed $l$ (because it's fewer numbers, or because the forbidden ones are
  "supposed to be zero") throws away the project's actual instrument.
- Pure NumPy/SciPy for now; nothing here depends on JAX or on defumat's code, only, eventually,
  on its *output*, once P7 is reached.

## What a finished phase in PLAN.md leaves behind

A phase is not done when the code runs. It needs, in this order:

1. **A number** — against a hand-derived case, a second code, or an identity that shares no
   machinery with the thing being checked. "It looks right" is not a phase result.
2. **The independent check named**, specifically enough that someone unfamiliar with the code
   could tell whether it was actually independent (see the rule above).
3. **The traps it ran into recorded**, even the ones that turned out to be nothing — a trap
   that was checked and ruled out is exactly as valuable a sentence as one that was found real,
   because the next person doesn't have to re-check it.

## Reference

The physics brief and the three routes (Quantum ESPRESSO non-self-consistent on an explicit
circle, Wannier90 from `hr.dat`, [defumat](https://github.com/joselado/defumat)'s differentiable
JAX plane-wave code) are in `README.md`. Numerical traps worth re-reading before writing
anything that touches real k-points: self-consistency noise, basis-set discontinuity from
`ecutwfc`'s discrete plane-wave count, the Cartesian-angle rule, and band index vs. band
identity — all four are in `README.md`'s "Known numerical traps" and are tracked with numbers
in `PLAN.md`.
