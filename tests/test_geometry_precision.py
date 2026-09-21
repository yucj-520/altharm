"""How accurately do the atoms have to sit before the harmonics stop meaning anything?

The structure came out of a VASP relaxation.  A relaxation converges forces, not
coordinates, and the usual criteria leave the atoms good to something like
1e-3 A.  But the operation that makes Cr2SeAs an altermagnet -- the horizontal
mirror M_z, which maps each Cr onto its opposite-spin partner -- only holds to
1.6e-5 A in the structure as delivered.  There are two orders of magnitude
between "as accurate as the relaxation guarantees" and "as accurate as the
symmetry needs", and nothing in the pipeline had measured what lives in that
gap.  This file measures it, on both sides:

**Side one, free: what the symmetry classification does.**  Displace every atom
by a Cartesian random vector of fixed seed and per-component RMS `amp`, and ask
spglib (through `altharm.symmetry`) what group it sees.

    amp (A)     group at symprec=1e-4   allowed l
    0           P-6m2, 6 exchanging     {3,9,15,21,27,33,39,45}   (8 of them)
    1e-5        P-6m2, 6 exchanging     the same 8
    2e-5        Cm,    0 exchanging     all 49 up to lmax=48
    >= 5e-5     P1,    0 exchanging     all 49

The exact breakdown is at **amp = 1.99e-5 A**, and it is a cliff, not a slope:
the allowed set never loses two or three members, it goes from 8 straight to
everything.  `spin_operations` refuses rather than returning the 49, which is
the one thing that stops this being silent.  The tolerance that survives an
amplitude `amp` is symprec ~ 1.4 x (worst atom-mapping residual) ~ 5.4 x amp, so
a structure good to 1e-3 A needs symprec >= 5.4e-3 A -- 500x spglib's default.

**Side two, DFT: what the coefficients do.**  Three isotropic amplitudes plus
two anisotropic controls, each a full SCF + 152-point bands run whose input
differs from the sm12 baseline in the ATOMIC_POSITIONS block and nowhere else.
See `test_the_perturbed_inputs_differ_only_in_the_positions`.

**The DFT half, on the same five structures** (12x12x4, degauss 0.0150, band 44,
152-point circle; every input differs from the sm12 baseline in ATOMIC_POSITIONS
and nowhere else):

    amp (A)   worst allowed-l shift   noise floor   l = 0
    0         --                        24.3 ueV     10.7 ueV
    1e-4      0.36% (l=15)              29.4 ueV     76.3 ueV
    1e-3      3.70% (l=15)             166.4 ueV    676.0 ueV
    1e-2     34.8%  (l=15)            1695.3 ueV   6708.4 ueV

**The answer to the question this file asks:** relaxation precision (~1e-3 A) is
good enough for the coefficients -- the strong harmonics move under 0.3% and
even the weakest moves under 4% -- and **not good enough for the error bar**,
which degrades 7x.  The floor is what tells you how much to trust a coefficient,
so a structure at 1e-3 A gives you numbers you can mostly believe and no honest
way to say so.  Getting the floor back to its intrinsic ~24 ueV needs ~1e-4 A,
which a relaxation does not deliver: the structure has to be **symmetrised**
before the harmonics are fitted.

And l = 0 is the sharpest single number here.  It is the angular average of the
splitting, which a compensated altermagnet must have at exactly zero.  At 1e-3 A
it reads 676 ueV.  Nothing in the fit objects; only comparing it against zero does.

**The forbidden spectrum says WHICH symmetry broke**, which is more than an
error bar.  The two anisotropic controls separate the two generators cleanly:

    displacement     worst forbidden   even-multiples-of-3   non-multiples-of-3   l = 0
    z only,   1e-3   l=12, 340 ueV          167.8 ueV              29.9 ueV      536 ueV
    in-plane, 1e-3   l=2,  783 ueV            2.2 ueV             180.3 ueV       32 ueV

Moving atoms along z breaks M_z -- the mirror that exchanges the magnetic
sublattices -- so the 'l is odd' half of the rule fails and even multiples of 3
leak in, l = 0 among them.  Moving them in-plane breaks C3 instead, so the
'multiple of 3' half fails and non-multiples leak in, while the sublattice
compensation survives (l = 0 stays 17x smaller).  This mirrors the spglib result
exactly: z-noise gives P3m1 with six preserving and no exchanging operation
(17 allowed l), in-plane noise gives P1 (49).  The two halves of the selection
rule are broken by orthogonal displacements, and the fit can tell them apart.

A corollary worth stating: in-plane noise leaves every allowed coefficient
within 0.18% of baseline while multiplying the floor by 6.7.  Coefficient
stability is not evidence that a structure is symmetric enough.

**Independent check.**  The symmetry side does not take spglib's word for where
the cliff is: the test computes, by hand and with no spglib involved, the worst
Cartesian distance between g(atom) and its intended partner for each of the
twelve operations, and confirms that spglib's bisected symprec threshold equals
that residual (ratio 1.00 for the ideal structure, 1.37 once perturbed).  The
DFT side is independent by construction -- pw.x's eigenvalues know nothing about
`altharm.symmetry`, and the allowed set they are scored against is P1's
hand-derived {3,9,...,45}, hard-coded here rather than recomputed.
"""

from pathlib import Path

import numpy as np
import pytest

from altharm.basis import CircularBasis
from altharm.noise import noise_report
from altharm.project import project
from altharm.qe.reader import read_bands
from altharm.symmetry import (
    SymmetryToleranceError,
    allowed_l,
    spin_operations,
    symprec_scan,
)

HERE = Path(__file__).parent

# ------------------------------------------------------------------ the structure
# The group's relaxed Cr2SeAs cell, exactly as it goes into scf.in.
CELL = np.array([
    [3.7258396285, 0.0000002055, 0.0000000080],
    [-1.8629196361, 3.2266718606, -0.0000000204],
    [0.0000000476, -0.0000000907, 16.9179784926],
])
POSITIONS = np.array([
    [0.3333333331, 0.6666666714, 0.5896431867],
    [0.3333333340, 0.6666666619, 0.4103569463],
    [0.6666666651, 0.3333333405, 0.0896259946],
    [0.6666666703, 0.3333333398, 0.9103736386],
    [0.0, 0.0, 0.3332384847],
    [0.0, 0.0, 0.6667617671],
    [0.3333333245, 0.6666666877, 0.8334101415],
    [0.3333333154, 0.6666666914, 0.1665899648],
    [0.6666666980, 0.3333332982, 0.5000002984],
    [0.0, 0.0, 0.0],
])
NUMBERS = np.array([24, 24, 24, 24, 34, 34, 34, 34, 33, 33])
MAGMOMS = np.array([+3.0, -3.0, -3.0, +3.0, 0, 0, 0, 0, 0, 0])

SEED = 20250921          # the displacements below are reproducible from this alone
WORKING_SYMPREC = 1e-4   # what tests/test_cr2seas_regression.py uses

BAND, LMAX = 44, 45
# P1's hand-derived selection rule, NOT recomputed here: C_3z preserves the
# sublattices so l must be a multiple of 3; M_z x time reversal gives C_2z with
# character -1, so l must be odd.
ALLOWED = (3, 9, 15, 21, 27, 33, 39, 45)


def displacements(amp, mode="iso", seed=SEED):
    """Cartesian displacements in Angstrom, per-component RMS exactly ``amp``.

    The full 10x3 normal sample is always drawn, so ``mode`` selects a subspace
    of the *same* random numbers rather than a different random sequence.  The
    mean is removed because a rigid translation breaks no symmetry -- spglib
    chooses its own origin -- and would otherwise waste amplitude on a null
    direction.
    """
    rng = np.random.default_rng(seed)
    u = rng.normal(size=(len(POSITIONS), 3))
    if mode == "z":
        u[:, :2] = 0.0
    elif mode == "xy":
        u[:, 2] = 0.0
    elif mode != "iso":
        raise ValueError(mode)
    u -= u.mean(axis=0)
    live = u[:, {"iso": slice(0, 3), "z": slice(2, 3), "xy": slice(0, 2)}[mode]]
    return u * (amp / np.sqrt((live ** 2).mean()))


def perturbed(amp, mode="iso", seed=SEED):
    """Fractional positions with that displacement applied."""
    if amp == 0.0:
        return POSITIONS.copy()
    return (POSITIONS @ CELL + displacements(amp, mode, seed)) @ np.linalg.inv(CELL)


# ------------------------------------------------- an independent geometric ruler

def _parent_operations():
    """(R, t) of the twelve operations of the ideal structure's magnetic parent."""
    return [(op.rot_frac, op.trans) for op in spin_operations(
        CELL, POSITIONS, NUMBERS, MAGMOMS,
        symprec=WORKING_SYMPREC, include_time_reversal=False)]


def mapping_residual(positions, ops):
    """Worst Cartesian distance, in Angstrom, from g(atom) to its nearest partner.

    Plain geometry: apply the operation, wrap into the cell, take the distance to
    the closest atom of the same species, keep the worst over atoms and
    operations.  spglib is not consulted; this is the quantity a tolerance has to
    be larger than for the operation to be accepted at all.
    """
    worst = 0.0
    for rot, trans in ops:
        moved = positions @ np.asarray(rot, dtype=float).T + np.asarray(trans)
        for i in range(len(positions)):
            d = moved[i] - positions
            d -= np.rint(d)
            dist = np.linalg.norm(d @ CELL, axis=1)
            worst = max(worst, float(np.where(NUMBERS == NUMBERS[i], dist, np.inf).min()))
    return worst


def _true_group(positions, symprec):
    row = symprec_scan(CELL, positions, NUMBERS, MAGMOMS, symprecs=(symprec,))[0]
    return row["spacegroup"] == "P-6m2" and row["n_exchanging"] == 6


def _threshold(predicate, lo=1e-8, hi=1.0, rising=False):
    """Bisect in log space for where ``predicate`` flips."""
    for _ in range(50):
        mid = np.sqrt(lo * hi)
        if predicate(mid) is rising:
            hi = mid
        else:
            lo = mid
    return hi if rising else lo


# ----------------------------------------------------------- the symmetry cliff

def test_the_structure_own_mirror_residual_is_1_6e_5_angstrom():
    """What the geometry is actually good to, measured without spglib.

    PLAN.md quotes 1.0e-5 A for this; that is the As atom at z = 0.5000002984.
    The worst pair is the Cr at z = 0.0896259946 / 0.9103736386, at 1.63e-5 A.
    Either way it sits on top of spglib's 1e-5 default, which is the trap.
    """
    r = mapping_residual(POSITIONS, _parent_operations())
    assert r == pytest.approx(1.630e-5, rel=0.02)


def test_spglib_threshold_is_exactly_that_residual():
    """The independent check: a hand-computed distance predicts spglib's cliff.

    For the ideal structure the bisected symprec at which P-6m2 appears equals
    the hand-measured residual to within 0.1%.  Once the structure is perturbed
    spglib needs about 1.4x the residual, because it refines the lattice and the
    origin as well as the atom mapping.
    """
    ops = _parent_operations()
    ideal = _threshold(lambda sp: _true_group(POSITIONS, sp), rising=True)
    assert ideal / mapping_residual(POSITIONS, ops) == pytest.approx(1.00, abs=0.05)

    for amp in (1e-4, 1e-3, 1e-2):
        pos = perturbed(amp)
        got = _threshold(lambda sp: _true_group(pos, sp), rising=True)
        assert got / mapping_residual(pos, ops) == pytest.approx(1.37, abs=0.06), amp


def test_the_allowed_set_collapses_at_2e_5_angstrom():
    """At the repo's working symprec = 1e-4, the largest displacement that keeps
    the true group is 1.99e-5 A -- fifty times finer than a relaxation delivers."""
    amp = _threshold(lambda a: _true_group(perturbed(a), WORKING_SYMPREC))
    assert amp == pytest.approx(1.99e-5, rel=0.05)

    assert _true_group(perturbed(1e-5), WORKING_SYMPREC)
    assert not _true_group(perturbed(2e-5), WORKING_SYMPREC)


def test_loosening_symprec_buys_back_exactly_what_the_ruler_predicts():
    """symprec = 1e-3 tolerates 1.85e-4 A; the default 1e-5 tolerates nothing.

    The default fails on the *unperturbed* structure -- its own 1.63e-5 A
    residual already exceeds it -- which is why the bisection bottoms out.
    """
    tolerated = {sp: _threshold(lambda a: _true_group(perturbed(a), sp))
                 for sp in (1e-5, 5e-5, 1e-4, 1e-3)}
    assert tolerated[1e-5] < 1e-7
    assert tolerated[5e-5] == pytest.approx(1.06e-5, rel=0.05)
    assert tolerated[1e-4] == pytest.approx(1.99e-5, rel=0.05)
    assert tolerated[1e-3] == pytest.approx(1.85e-4, rel=0.05)


def _n_allowed(positions, symprec=WORKING_SYMPREC, lmax=48):
    basis = CircularBasis(lmax)
    b = 2.0 * np.pi * np.linalg.inv(CELL).T
    from altharm.kpath import circle
    pts = circle(0.25 * np.linalg.norm(b[0]), 4 * lmax + 8,
                 center=(0.0, 0.0, 0.25 * np.linalg.norm(b[2])))
    ops = spin_operations(CELL, positions, NUMBERS, MAGMOMS,
                          symprec=symprec, require_exchange=False)
    return len(allowed_l(basis, ops, pts))


def test_isotropic_noise_collapses_the_set_in_one_step():
    """8 allowed l or 49, never anything between, when the noise is isotropic.

    There is no regime where a slightly-wrong geometry quietly admits two extra
    harmonics.  The six sublattice-exchanging operations go together, the three
    rotations go with them, and `allowed_l` returns the whole of l = 0..lmax.
    """
    counts = [_n_allowed(perturbed(a))
              for a in (0.0, 1e-5, 1.8e-5, 2e-5, 1e-4, 1e-3, 1e-2)]
    assert counts == [8, 8, 8, 49, 49, 49, 49]


def test_z_only_noise_fails_differently_and_worse():
    """Displacing only z keeps C_3z and kills only M_z: 8 allowed l become 17.

    M_z is the operation that exchanges the magnetic sublattices, and it acts on
    z, so a purely out-of-plane error is the one thing that destroys it while
    leaving the three-fold rotation untouched.  What survives is 'l is a multiple
    of 3' without 'l is odd', i.e. every multiple of 3 from 0 to lmax.

    **l = 0 becoming allowed is the whole danger**: l = 0 is the angular average
    of the splitting, so a compensated altermagnet is now formally permitted to
    carry a net moment's worth of splitting, and nothing in the fit flags it.
    This is the same wrong answer the symprec cliff produces on the unperturbed
    structure -- P3m1, six preserving operations, no exchanging one -- reached
    here through the geometry instead of through the tolerance.
    """
    assert _n_allowed(perturbed(1e-3, mode="z")) == 17
    assert _n_allowed(perturbed(1e-3, mode="xy")) == 49

    row = symprec_scan(CELL, perturbed(1e-3, mode="z"), NUMBERS, MAGMOMS,
                       symprecs=(WORKING_SYMPREC,))[0]
    assert (row["spacegroup"], row["n_preserving"], row["n_exchanging"]) == ("P3m1", 6, 0)


def test_the_collapse_is_refused_by_name_rather_than_returned():
    """The guard that keeps this from being silent."""
    assert sorted(allowed_l(
        CircularBasis(LMAX),
        spin_operations(CELL, perturbed(1e-5), NUMBERS, MAGMOMS,
                        symprec=WORKING_SYMPREC),
        _circle_dirs())) == list(ALLOWED)

    with pytest.raises(SymmetryToleranceError, match="symprec_scan"):
        spin_operations(CELL, perturbed(1e-3), NUMBERS, MAGMOMS,
                        symprec=WORKING_SYMPREC)


def _circle_dirs():
    b = 2.0 * np.pi * np.linalg.inv(CELL).T
    from altharm.kpath import circle
    return circle(0.25 * np.linalg.norm(b[0]), 4 * LMAX + 8,
                  center=(0.0, 0.0, 0.25 * np.linalg.norm(b[2])))


# ------------------------------------------------------------- the DFT half

_GEOM = {
    "baseline": "qe_cr2seas_sm12",
    "1e-4": "qe_cr2seas_geom1e4",
    "1e-3": "qe_cr2seas_geom1e3",
    "1e-2": "qe_cr2seas_geom1e2",
    "z1e-3": "qe_cr2seas_geomz1e3",
    "xy1e-3": "qe_cr2seas_geomxy1e3",
}
_DATA = Path(__file__).parent / "data"
BAND = 44

geomtests = pytest.mark.skipif(
    not all((_DATA / d / "data-file-schema.xml").exists() for d in _GEOM.values()),
    reason="the perturbed-geometry fixtures are not present",
)


def _harm(tag):
    q = read_bands(_DATA / _GEOM[tag] / "data-file-schema.xml")
    assert len(q.kpoints_tpiba) == 152
    th = np.arctan2(q.kpoints_tpiba[:, 1], q.kpoints_tpiba[:, 0])
    dirs = np.stack([np.cos(th), np.sin(th), 0.093 * np.ones_like(th)], axis=1)
    f = project(CircularBasis(LMAX), dirs,
                q.energies_ev[0, :, BAND] - q.energies_ev[1, :, BAND])
    return f.amplitude_by_l(), noise_report(f, set(ALLOWED))


@geomtests
def test_the_perturbed_inputs_differ_only_in_the_positions():
    """Guard the premise: one knob, not several."""
    base = (_DATA / "qe_cr2seas_sm12" / "scf.in").read_text().splitlines()
    for tag, folder in _GEOM.items():
        if tag == "baseline":
            continue
        other = (_DATA / folder / "scf.in").read_text().splitlines()
        differing = [(a, b) for a, b in zip(base, other) if a != b]
        assert len(base) == len(other)
        # every differing line is a coordinate line: species + 3 floats
        for a, b in differing:
            assert a.split()[0] in ("CrU", "CrD", "Se", "As"), f"{folder}: {a!r}"
            assert len(a.split()) == 4
        assert "degauss     = 0.0150" in "\n".join(other)
        assert "12 12 4" in "\n".join(other)


@geomtests
def test_one_angstrom_in_ten_thousand_is_harmless_to_the_coefficients():
    base, _ = _harm("baseline")
    a, _ = _harm("1e-4")
    worst = max(abs(a[l] - base[l]) / base[l] for l in ALLOWED)
    assert worst < 0.01, f"1e-4 A moved a coefficient by {worst:.3f}"


@geomtests
def test_relaxation_precision_is_fine_for_the_coefficients():
    """1e-3 A: strong harmonics under 0.3%, the weakest under 4%."""
    base, _ = _harm("baseline")
    a, _ = _harm("1e-3")
    assert abs(a[3] - base[3]) / base[3] < 0.005
    assert abs(a[9] - base[9]) / base[9] < 0.005
    worst_l = max(ALLOWED, key=lambda l: abs(a[l] - base[l]) / base[l])
    assert worst_l == 15
    assert abs(a[15] - base[15]) / base[15] == pytest.approx(0.037, abs=0.01)


@geomtests
def test_relaxation_precision_destroys_the_error_bar():
    """The finding: the floor degrades 7x while the coefficients barely move.

    This is the practical conclusion of the whole file.  A structure straight
    out of a relaxation gives coefficients you can mostly believe and no honest
    way to say how much.
    """
    _, base = _harm("baseline")
    _, nr = _harm("1e-3")
    assert base.floor == pytest.approx(24.3e-6, rel=0.05)
    assert nr.floor == pytest.approx(166e-6, rel=0.05)
    assert nr.floor / base.floor > 5


@geomtests
def test_l0_is_the_sharpest_single_indicator():
    """l = 0 must be exactly zero for a compensated altermagnet.  It is not."""
    base, _ = _harm("baseline")
    assert base[0] < 20e-6
    for tag, expect in (("1e-4", 76e-6), ("1e-3", 676e-6), ("1e-2", 6708e-6)):
        a, _ = _harm(tag)
        assert a[0] == pytest.approx(expect, rel=0.08), tag


@geomtests
def test_ten_times_worse_geometry_breaks_everything():
    base, _ = _harm("baseline")
    a, nr = _harm("1e-2")
    assert abs(a[15] - base[15]) / base[15] > 0.30
    assert nr.floor > 1e-3


@geomtests
def test_the_forbidden_spectrum_identifies_which_symmetry_broke():
    """The result that is more than an error bar.

    z-displacements break M_z (sublattice exchange) -> the 'l odd' half of the
    rule fails -> even multiples of 3 leak in, l = 0 among them.
    In-plane displacements break C3 -> the 'multiple of 3' half fails ->
    non-multiples leak in, and the compensation survives.
    """
    def split(tag):
        _, nr = _harm(tag)
        even3 = [l for l in nr.forbidden if l % 3 == 0]
        rest = [l for l in nr.forbidden if l % 3 != 0]
        rms = lambda ls: float(np.sqrt(np.mean([nr.forbidden[l] ** 2 for l in ls])))
        return rms(even3), rms(rest), max(nr.forbidden, key=nr.forbidden.get)

    z3, zr, z_worst = split("z1e-3")
    x3, xr, x_worst = split("xy1e-3")

    # z-noise: the multiples-of-3 channel dominates, and the worst is even
    assert z3 / zr > 4
    assert z_worst == 12 and z_worst % 2 == 0
    # in-plane noise: the other channel dominates, by two orders of magnitude
    assert xr / x3 > 50
    assert x_worst == 2
    # and only z-noise lifts l = 0 appreciably: the compensation is M_z's doing
    az, _ = _harm("z1e-3")
    ax, _ = _harm("xy1e-3")
    assert az[0] / ax[0] > 10


@geomtests
def test_stable_coefficients_do_not_prove_a_symmetric_structure():
    """In-plane noise: every coefficient within 0.18%, floor up 6.7x."""
    base, nbase = _harm("baseline")
    a, nr = _harm("xy1e-3")
    assert max(abs(a[l] - base[l]) / base[l] for l in ALLOWED) < 0.003
    assert nr.floor / nbase.floor > 5
