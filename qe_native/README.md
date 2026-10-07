# harmwave.x — a QE-native harmonic-fit tool

A new compiled Quantum ESPRESSO post-processing executable, built the way
`bands.x`/`dos.x` are: linked against QE's own `PP/src` and `PW/src`
libraries, reading a finished `calculation='bands'` run via QE's own
`read_file()`. It computes two things directly from QE's own in-memory state,
without ever touching an XML file by hand:

1. **The magnetic character of every symmetry operation.** The two (or more)
   magnetic sublattices are collapsed back to their true chemical species
   (read from the pseudopotential's own `psd` field, not from `ityp` — this
   group's convention gives each sublattice its own `ATOMIC_SPECIES` entry so
   LSDA can start them with opposite `starting_magnetization`, and QE's
   `find_sym()` would otherwise never let an operation map one onto the
   other). QE's own symmetry finder is called on that collapsed structure;
   each operation is then tagged `chi = +1` (preserves the sublattices) or
   `chi = -1` (exchanges them, using the `irt` atom-mapping array `find_sym`
   already computes), closed under time reversal, and filtered to the ones
   that leave the requested circle's plane fixed.

2. **A raw least-squares harmonic fit**, `dE(theta) = sum_l a_l cos(l theta)
   + b_l sin(l theta)`, computed by LAPACK's `DGELSS` directly on the
   Rydberg-unit eigenvalues QE holds in `wvfct`'s `et` array.

Both are written to a plain-text report. Which harmonics that operation list
*allows*, and the noise-floor separation between allowed and forbidden
amplitude, are deliberately left to the already-tested Python code in
`altharm.symmetry`/`altharm.noise` — see `tests/test_harmwave_cross_check.py`.

## Why this exists

The rest of `altharm` computes the same quantities by parsing QE's XML output
in Python and finding the crystal symmetry with spglib. `harmwave.x` is an
**independent second implementation**: different language, QE's own
symmetry engine instead of spglib, LAPACK instead of numpy, an in-memory
Rydberg-to-eV conversion instead of a Hartree-to-eV one parsed from XML.
Agreement between the two is a real check, not a tautology — and finding
where they *disagree* is exactly what a first version of this tool did.

## What building it caught

A first version compiled cleanly and ran to completion, but its
time-reversal closure only doubled operations that had **already** survived
the plane-preservation filter for the given circle. The operation that
actually makes Cr2Se2As an altermagnet — a horizontal mirror that exchanges
the two Cr sublattices — flips k_z on its own and so fails that filter
directly; its time-reversed partner (which restores k_z) needs the mirror
included in a closure taken over the *full* symmetry group, before
filtering, not after. Filtering first and closing second silently dropped
it, and the tool reported the allowed harmonics as {0, 3, 6, 9, 12} instead
of the correct {3, 9, 15, ...} — over-inclusive, not wrong in a way that
would have been obvious from the numbers alone, since l=0/6/12 still came
out numerically tiny (noise-floor-sized) in the actual fit. Only feeding the
tool's own operator dump into `altharm.symmetry.allowed_harmonics()` and
comparing the *predicted* selection rule against where the real signal
landed caught it. Fixed by closing under time reversal over the full
operation set first, filtering for the circle's plane second (see the
`keep_if_plane_preserving` internal subroutine and its docstring).

## Status

Validated once, on the group's real Cr2Se2As run (PBE+U, USPP, 56-point
circle at k_z = 0.25 b3, band 44 -- the one nearest E_F, at E-E_F = +0.04 eV):
eigenvalues agree with the Python XML route to 1e-9 eV, the fitted
coefficients agree to 1e-9 eV coefficient-by-coefficient, and the operation
list — fed into `altharm.symmetry` — reproduces the {3, 9} selection rule
exactly where the real signal is (l=3: 59.13 meV, l=9: 38.86 meV, agreeing
with the group's independent Wannier fit to under 1%). See
`tests/test_harmwave_cross_check.py` and `tests/data/harmwave_cr2seas/`.

A first run of this validation used band index 22 by mistake -- copied from
the Wannier route's 38-band numbering, where 22 *is* the band nearest E_F,
but this DFT run has 80 bands (22 more, matching the group's own "16
semicore + 6 Se/As-4s" comment in `nscf.in`), so band 22 here is a
dispersionless semicore state 4.4 eV below E_F. Caught by checking the
band's mean energy against E_F, which is worth doing on every new
`band_index` rather than assumed from a different model's numbering.

**Scope of this v1, honestly:** collinear LSDA only (no SOC, no
noncollinear); the circle must lie in the k_x k_y plane at a fixed k_z
(hardcoded, not general); the band index is given explicitly rather than
tracked across k-points (see `altharm.bandtrack` for why that matters and
is not yet ported here); `magmom` per species is given explicitly in the
tool's own namelist rather than read back from the run (a deliberate choice,
not a shortcut — see the comment in `harmwave.f90`). Symmetrization/noise
floor is not (yet) computed in Fortran; the Python route already does it
correctly and there is no need to duplicate that risk.

## Building

The Makefile change is committed as `harmwave-makefile.patch` rather than
described in prose, because "copy the bands.x target" is the kind of
instruction that is followed slightly wrong.  It adds five lines and applies
cleanly to **QE 7.2** (`qe-7.2`, commit `b4c5e8d`); on another release, apply it
by hand -- it is a verbatim copy of the `bands.x` rule with the name changed.

Linking needs FFTW on the path: without `module load fftw/3.3.10` (or the
equivalent) the compile succeeds and the link fails with
`cannot find -lfftw3`.

Requires a full QE 7.2 source tree (not just the compiled binaries a module
system provides) and the same toolchain used elsewhere in this project:

```bash
module load triton/2024.1-gcc openmpi/4.1.6 quantum-espresso/7.2
git clone --branch qe-7.2 --depth 1 https://gitlab.com/QEF/q-e.git
cd q-e
GF=$(mpif90 -show | grep -oE '^\S+gfortran\S*' || true)   # or locate it directly
./configure MPIF90=mpif90 CC=mpicc FC="$GF" F90="$GF" \
    BLAS_LIBS="-lopenblas" LAPACK_LIBS="-lopenblas" SCALAPACK_LIBS="-lscalapack"
make -j8 pw pp
cp /path/to/altharm/qe_native/harmwave.f90 PP/src/
git apply /path/to/altharm/qe_native/harmwave-makefile.patch
make harmwave.x
```

`quantum_espresso/` (the vendored QE source tree) is gitignored, like
`defumat`'s own vendored QE copy — it does not belong in this repository's
history. This file is the one QE-native source file that *is* committed,
since it is the actual deliverable; the surrounding QE tree is reference
material to build against.

## Running

```
&inputpp
  prefix      = 'AB'          ! same prefix as the bands.in that produced it
  outdir      = '../tmp'      ! points at prefix.save
  lmax        = 12
  band_index  = 45             ! 1-based; Python's 0-based index + 1.
                                ! Check the band's mean energy against E_F --
                                ! do not carry an index over from a different
                                ! model's band count (see Status above).
  magmom(1)   = 1.0            ! one entry per ATOMIC_SPECIES, same order/signs
  magmom(2)   = -1.0           ! as the bands.in that produced the run
  magmom(3)   = 0.0
  magmom(4)   = 0.0
  moment_tol  = 1.0d-3
  report_file = 'harmwave_report.txt'
/
```

```bash
harmwave.x -in harm.in
```
