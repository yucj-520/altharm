"""Regression on the group's real Cr2SeAs Wannier model.

Independence rule, in its strongest form here: the symmetry module predicts
which harmonics may be non-zero from the structure alone, and the Wannier
Hamiltonian -- fitted to PBE+U bands, knowing nothing about that prediction --
supplies the numbers.  Neither was derived from the other.

These numbers are the record.  If one moves, either the physics changed or
something broke, and either way it must be explained rather than re-baselined.
"""

from pathlib import Path

import numpy as np
import pytest

from altharm.basis import CircularBasis
from altharm.kpath import circle
from altharm.project import project
from altharm.noise import noise_report
from altharm.symmetry import spin_operations, allowed_l, symprec_scan
from altharm import wannier

ROOT = "/scratch/work/yuc11/HIGH_THROUGHT_WANNIERcr2X3/1.cr2se2as/Cr2SeAs_wannier_v03"
HR_UP = f"{ROOT}/hamiltonian/c.up_hr.dat"
HR_DN = f"{ROOT}/hamiltonian/c.dn_hr.dat"
E_FERMI = 8.1169

pytestmark = pytest.mark.skipif(
    not __import__("os").path.exists(HR_UP),
    reason="the Cr2SeAs Wannier model is not on this machine",
)

CELL = np.array([
    [3.7258396285, 0.0000002055, 0.0000000080],
    [-1.8629196361, 3.2266718606, -0.0000000204],
    [0.0000000476, -0.0000000907, 16.9179784926],
])
POSITIONS = np.array([
    [1 / 3, 2 / 3, 0.5896431867], [1 / 3, 2 / 3, 0.4103569463],
    [2 / 3, 1 / 3, 0.0896259946], [2 / 3, 1 / 3, 0.9103736386],
    [0.0, 0.0, 0.3332384847], [0.0, 0.0, 0.6667617671],
    [1 / 3, 2 / 3, 0.8334101415], [1 / 3, 2 / 3, 0.1665899648],
    [2 / 3, 1 / 3, 0.5000002984], [0.0, 0.0, 0.0],
])
NUMBERS = np.array([24, 24, 24, 24, 34, 34, 34, 34, 33, 33])
MAGMOMS = np.array([+3.0, -3.0, -3.0, +3.0, 0, 0, 0, 0, 0, 0])
SYMPREC = 1e-4


def reciprocal():
    return 2.0 * np.pi * np.linalg.inv(CELL).T


@pytest.fixture(scope="module")
def model():
    return (wannier.read_hr(HR_UP, cell=CELL), wannier.read_hr(HR_DN, cell=CELL))


@pytest.fixture(scope="module")
def ops():
    return spin_operations(CELL, POSITIONS, NUMBERS, MAGMOMS, symprec=SYMPREC)


def geometry(lmax, kz_frac):
    B = reciprocal()
    pts = circle(0.25 * np.linalg.norm(B[0]), 4 * lmax + 8,
                 center=(0.0, 0.0, kz_frac * np.linalg.norm(B[2])))
    return CircularBasis(lmax), pts


def band_at_fermi(up, frac):
    return int(np.argmin(np.abs(wannier.eigenvalues(up, frac).mean(axis=0) - E_FERMI)))


# ------------------------------------------------------------ the symprec cliff

def test_the_parent_group_changes_between_1e_5_and_5e_5():
    """Measured: P3m1 with no exchanging operation below, P-6m2 with six above.

    The horizontal mirror maps each atom onto its partner to within 1.0e-5 A,
    which is spglib's default tolerance, so the default silently returns a
    structure that cannot be an altermagnet.
    """
    rows = {r["symprec"]: r for r in symprec_scan(CELL, POSITIONS, NUMBERS, MAGMOMS)}
    assert rows[1e-5]["spacegroup"] == "P3m1"
    assert rows[1e-5]["n_exchanging"] == 0
    assert rows[5e-5]["spacegroup"] == "P-6m2"
    assert rows[5e-5]["n_exchanging"] == 6


def test_operation_count_at_the_working_tolerance(ops):
    assert len(ops) == 24
    assert sum(o.chi == +1 for o in ops) == 12
    assert sum(o.chi == -1 for o in ops) == 12


# --------------------------------------------------------------- the selection

def test_nothing_is_allowed_in_the_k_z_zero_plane(ops):
    """M_z exchanges the sublattices and acts as k_z -> -k_z, so at k_z = 0 it
    forces dE = -dE.  The splitting vanishes identically in that plane."""
    basis, pts = geometry(12, 0.0)
    assert allowed_l(basis, ops, pts) == set()


def test_allowed_harmonics_are_the_odd_multiples_of_three(ops):
    """C_3z preserves the sublattices, so l must be a multiple of 3.  M_z
    combined with time reversal gives C_2z with character -1, so l must be odd.
    What survives is 3, 9, 15, ..."""
    basis, pts = geometry(48, 0.25)
    assert sorted(allowed_l(basis, ops, pts)) == [3, 9, 15, 21, 27, 33, 39, 45]


# ------------------------------------------------------------------ the numbers

def test_the_k_z_zero_plane_measures_the_wannier_model_own_error(model, ops):
    """Symmetry says zero; what the model returns instead is its angular error.

    Measured: every harmonic below 0.07 meV, consistent with the 0.06 meV RMS
    the model reports against DFT over its main window.  This is the floor the
    DFT route will have to beat to say anything new.
    """
    up, dn = model
    basis, pts = geometry(12, 0.0)
    frac = wannier.cart_to_frac_k(pts, CELL)
    band = band_at_fermi(up, frac)
    dE = wannier.spin_splitting(up, dn, frac, band)

    assert np.abs(dE).max() < 0.25e-3, f"{np.abs(dE).max() / 1e-3:.3f} meV at k_z = 0"
    rep = noise_report(project(basis, pts, dE), allowed_l(basis, ops, pts))
    assert rep.floor_max < 0.07e-3, f"worst harmonic {rep.floor_max / 1e-3:.4f} meV"


def test_allowed_harmonics_dominate_the_forbidden_ones(model, ops):
    """The cross-check: the two largest harmonics must be the two symmetry allows.

    Measured at k_z = 0.25 b3 on the band nearest E_F: l = 3 at 59.5 meV and
    l = 9 at 39.1 meV, with every forbidden harmonic below 1.2 meV.
    """
    up, dn = model
    basis, pts = geometry(12, 0.25)
    frac = wannier.cart_to_frac_k(pts, CELL)
    band = band_at_fermi(up, frac)
    dE = wannier.spin_splitting(up, dn, frac, band)

    rep = noise_report(project(basis, pts, dE), allowed_l(basis, ops, pts))
    assert set(rep.allowed) == {3, 9}
    assert np.isclose(rep.allowed[3], 59.5e-3, rtol=0.02)
    assert np.isclose(rep.allowed[9], 39.1e-3, rtol=0.02)
    assert max(rep.forbidden.values()) < 1.2e-3
    assert min(rep.snr().values()) > 50.0


def test_high_harmonics_need_long_in_plane_hoppings(model):
    """What a short-ranged Wannier model would lose, measured.

    Measured, in-plane cut with the c axis left intact, on the l = 21 harmonic:
    1.33 meV at 4 A against 8.44 meV in the full model, a factor of six.  The
    low harmonic is wrong too, by about a third.  This is the quantitative form
    of "a Wannier model does not capture the high harmonics".
    """
    up, dn = model
    lmax = 48
    basis, pts = geometry(lmax, 0.25)
    frac = wannier.cart_to_frac_k(pts, CELL)
    band = band_at_fermi(up, frac)

    def amps(rin=None):
        u, d = (up, dn) if rin is None else (
            wannier.truncate(up, rin, plane_normal=[0, 0, 1]),
            wannier.truncate(dn, rin, plane_normal=[0, 0, 1]),
        )
        dE = wannier.spin_splitting(u, d, frac, band)
        return project(basis, pts, dE).amplitude_by_l()

    full, short, long_ = amps(), amps(4.0), amps(15.0)

    assert np.isclose(full[3], 59.5e-3, rtol=0.02)
    assert np.isclose(full[9], 39.1e-3, rtol=0.02)
    assert np.isclose(full[21], 8.44e-3, rtol=0.03)

    assert short[21] < full[21] / 5.0, "a 4 A model should lose most of l = 21"
    assert abs(short[3] - full[3]) > 0.25 * full[3], "and get l = 3 badly wrong"

    # by 15 A everything is converged to better than 1 %
    for l in (3, 9, 21):
        assert abs(long_[l] - full[l]) < 0.01 * full[l], (
            f"l={l} not converged at 15 A: {long_[l] / 1e-3:.3f} vs {full[l] / 1e-3:.3f} meV"
        )


# ------------------------------------------------------ P3: the noise floor

def test_p3_noise_floor_on_the_real_frozen_density_run():
    """The headline number of the whole project -- corrected once already.

    **This test replaces an earlier version that measured band index 22.**
    That was wrong: it is the band nearest E_F in the group's Wannier model,
    which has 38 bands, but this DFT run has 80 (16 semicore + 6 Se/As-4s +
    38 p-d + 20 buffer, per the group's own nscf.in comment), so band 22
    there is a deep, nearly flat semicore state 4.4 eV below E_F -- not the
    band the project cares about, and its splitting is trivially near zero,
    which is why the noise floor measured on it looked spectacular (0.106
    ueV). Caught by checking the band's mean energy against E_F, which the
    first version never did.

    Measured correctly, on band 44 (E - E_F = +0.04 eV, matching the
    Wannier route's own reported band to the same energy): the allowed
    harmonics agree with the independent Wannier-route fit to under 1%
    (l=3: 59.13 vs 59.48 meV; l=9: 38.86 vs 39.06 meV) -- a genuine,
    parameter-free cross-check between DFT and Wannier this project did not
    have before. But the "noise floor" picture is now much less clean: with
    only 56 k-points, the forbidden harmonics are NOT uniformly tiny. Every
    forbidden **even** l is at the ~0.1-1 ueV level (as clean as the
    wrong-band measurement suggested), but every forbidden **odd,
    non-multiple-of-three** l (1, 5, 7, 11) sits at 0.4-1.2 meV -- three
    orders of magnitude worse, and not negligible next to the smaller real
    signals this band carries at higher l (l=15 is 3.2 meV, l=21 is 8.5 meV,
    both confirmed against the same Wannier route). Whether that asymmetry
    is genuine DFT noise or an artefact of only having 56 points to resolve
    up to l=12 is P3's next open question, not yet answered.
    """
    from altharm.qe.reader import read_bands
    from altharm.symmetry import SpinOp, allowed_l

    report = Path(__file__).parent / "data" / "harmwave_cr2seas" / "harmwave_report.txt"
    xml = Path(__file__).parent / "data" / "qe_cr2seas_p3" / "data-file-schema.xml"
    if not report.exists() or not xml.exists():
        pytest.skip("real Cr2SeAs fixtures not present")

    ops = []
    for line in report.read_text().splitlines():
        if line.startswith("OP"):
            p = line.split()
            ops.append(SpinOp(rot_frac=None, trans=None, chi=int(p[1]),
                               rot_cart=np.array([float(x) for x in p[2:11]]).reshape(3, 3)))

    qb = read_bands(xml)
    band0 = 44
    assert abs((qb.energies_ev[0, :, band0].mean() - qb.fermi_ev) - 0.0415) < 0.001, (
        "band 44 should be the one at E - E_F = +0.04 eV; check nothing upstream moved"
    )
    dE = qb.energies_ev[0, :, band0] - qb.energies_ev[1, :, band0]
    theta = np.arctan2(qb.kpoints_tpiba[:, 1], qb.kpoints_tpiba[:, 0])
    dirs = np.stack([np.cos(theta), np.sin(theta), 0.093 * np.ones_like(theta)], axis=1)

    basis = CircularBasis(12)
    allowed = allowed_l(basis, ops, dirs)
    rep = noise_report(project(basis, dirs, dE), allowed)

    assert set(rep.allowed) == {3, 9}
    assert np.isclose(rep.allowed[3], 59.13e-3, rtol=0.01)
    assert np.isclose(rep.allowed[9], 38.86e-3, rtol=0.01)
    # the even forbidden harmonics: clean, ~ueV level
    even_forbidden = [rep.forbidden[l] for l in (2, 4, 6, 8, 10, 12)]
    assert max(even_forbidden) < 5e-6, f"even forbidden harmonics drifted: {even_forbidden}"
    # the odd, non-multiple-of-3 forbidden harmonics: not clean, ~meV level
    odd_forbidden = [rep.forbidden[l] for l in (1, 5, 7, 11)]
    assert 0.3e-3 < max(odd_forbidden) < 2e-3, (
        f"odd-non-3 forbidden harmonics moved outside their measured range: {odd_forbidden}"
    )


# ---------------------------------------------- P4: basis-set discontinuity

def test_p4_planewave_discontinuity_is_even_in_l_and_peaks_at_l6():
    """What the basis-set discontinuity can, and provably cannot, explain.

    **This test previously claimed the l=6 peak here explained P3's noise
    peak.** That claim was built on P3's since-corrected band (see
    test_p3_noise_floor_on_the_real_frozen_density_run's docstring): the
    wrong band's one substantial forbidden harmonic happened to sit at l=6,
    matching this spectrum's own l=6 peak by coincidence, not causation.

    On the corrected band, the real problem is the opposite parity: every
    forbidden EVEN l is clean (~ueV) and every forbidden ODD,
    non-multiple-of-three l (1, 5, 7, 11) sits at 0.4-1.2 meV. The
    measurement below -- the plane-wave count's own harmonic content is
    EXACTLY zero at every odd l -- therefore *rules out* the basis-set
    discontinuity as the explanation for that odd-l noise, rather than
    supplying one. That is still worth recording: a candidate mechanism
    checked and excluded is as much a result as one confirmed. What
    actually causes the odd-l noise is P3's open question.
    """
    from altharm.qe.reader import read_bands
    from altharm.planewaves import count_planewaves, discontinuity_spectrum
    from altharm.units import BOHR_ANG

    xml = Path(__file__).parent / "data" / "qe_cr2seas_p3" / "data-file-schema.xml"
    if not xml.exists():
        pytest.skip("real Cr2SeAs fixtures not present")

    qb = read_bands(xml)
    k_invang = qb.kpoints_tpiba * (2 * np.pi / (qb.alat_bohr * BOHR_ANG))
    npw = count_planewaves(k_invang, qb.cell_ang, ecutwfc_ry=60.0)

    assert npw.min() == 10770 and npw.max() == 10801

    theta = np.arctan2(k_invang[:, 1], k_invang[:, 0])
    basis = CircularBasis(12)
    dirs = np.stack([np.cos(theta), np.sin(theta), np.zeros_like(theta)], axis=1)
    amp = discontinuity_spectrum(basis, dirs, npw).amplitude_by_l()

    for l in (1, 3, 5, 7, 9, 11):
        assert amp[l] < 1e-8, f"odd l={l} should vanish exactly, got {amp[l]:.2e}"
    even_others = [amp[l] for l in (2, 4, 8, 10, 12)]
    assert amp[6] > 10 * max(even_others), (
        f"l=6 ({amp[6]:.2f}) should dominate the other even harmonics {even_others}"
    )
