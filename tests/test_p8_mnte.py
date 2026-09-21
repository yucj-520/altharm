"""P8, first real system beyond Cr2SeAs: MnTe, and the harmonic ladder to l = 45.

MnTe (NiAs-type, P6_3/mmc, Mn moments along c on alternating layers) is the
canonical bulk altermagnet, so unlike Cr2SeAs it has a literature classification
to be checked against rather than only an internal one.

**The symmetry side, before any DFT.** No symprec cliff -- P6_3/mmc with 12
sublattice-preserving and 12 exchanging operations at every tolerance from 1e-6
to 1e-2, where Cr2SeAs flips at 1e-5. Allowed set {3, 9, 15, 21, ...}, the odd
multiples of 3, and **nothing at all at k_z = 0**. That reproduces the textbook
g-wave form for MnTe, dE ∝ k_z k_x(k_x^2 - 3k_y^2) = k_z k_perp^3 cos(3θ):
l = 3 in plane, odd in k_z, and four nodal surfaces (three vertical planes plus
k_z = 0) which is what "g-wave" counts. Note that g-wave here means **l = 3**,
not l = 4 -- the labels count nodal surfaces, and the k_z = 0 plane is one.

**The DFT side** (12x12x8, degauss 0.0150, 152-point circle at k_z = 0.25 b3,
U = 4.0 eV on Mn-3d), band 20 = valence-band maximum:

    l        3       9      15      21      27      33      39      45
    meV  949.38   53.45   12.18    3.96    1.466   0.582   0.243   0.105
    /floor 4e7    2e6     5e5     2e5     6e4     2e4     1e4     4493

**All eight allowed harmonics are resolved**, the weakest at 4500x the noise
floor. Three things make that possible and each is worth separating:

1. MnTe is a **semiconductor**, so the Fermi-surface sampling problem that
   dominated Cr2SeAs (P3's smearing result) does not arise. Floor 0.02-0.03
   μeV against Cr2SeAs's 24 μeV -- three orders of magnitude.
2. Its atoms sit at **exact Wyckoff positions**, so there is no relaxation noise
   breaking the symmetry. The cos/sin ratio at l = 3 is **1e-13 to 4e-12**
   across both U values and both bands, against **8.6e-10** for relaxed Cr2SeAs
   at the same l -- 200x to 7000x cleaner. This is the structure-precision trap
   of `test_geometry_precision.py` seen from the good side.

   (An earlier draft of this docstring quoted 5e-4 for Cr2SeAs, which is its
   *worst* l at U = 2, not its l = 3. Comparing MnTe's l = 3 against another
   system's worst l overstated the gap by six orders of magnitude.)
3. Its splitting is large: 949 meV at l = 3, against Cr2SeAs's 64 meV.

**Band 20's ladder is exponential in l**: A ∝ exp(-l/6.35), R^2 = 0.997 over
four decades, fitting better than a power law (R^2 = 0.989). And it passes a
self-consistency check that does not involve the fit: the residual left after
l_max = 45 is 0.0374 meV, while summing the same exponential over l = 51..75
predicts 0.0425 meV -- 12% apart.

**This is a per-band property, not a law.** Band 21 (conduction-band minimum) is
not a clean exponential: its step ratios run 5.4, 10.1, 3.1, 1.25, and a power
law fits it marginally better. Whatever sets the decay constant is a feature of
the individual band, and one band is not enough to say what.
"""

from pathlib import Path

import numpy as np
import pytest

from altharm.basis import CircularBasis
from altharm.kpath import circle
from altharm.noise import noise_report
from altharm.project import project
from altharm.qe.reader import read_bands
from altharm.symmetry import allowed_l, spin_operations, symprec_scan

HERE = Path(__file__).parent
# One directory per circle, because they are separate bands runs over one density.
# The first attempt kept only `tmp/MT.save/`, so the second run silently
# overwrote the first's eigenvalues -- see `test_the_two_circles_are_separate_files`.
KZ0 = HERE / "data" / "mnte_kz0"                      # the nodal plane, U = 0
KZ0_U4 = HERE / "data" / "mnte_kz0_u4"                # the same plane at U = 4
RUNS = {0.0: HERE / "data" / "mnte_kz025_u0", 4.0: HERE / "data" / "mnte_kz025_u4"}

kz025 = pytest.mark.skipif(
    not all((d / "data-file-schema.xml").exists() for d in RUNS.values()),
    reason="the k_z != 0 MnTe fixtures are being recomputed",
)
kz0 = pytest.mark.skipif(not (KZ0 / "data-file-schema.xml").exists(),
                         reason="the k_z = 0 MnTe fixture is not present")

CELL = np.array([[4.2200465202, 0.0, 0.0],
                 [-2.1100232601, 3.6546674917, 0.0],
                 [0.0, 0.0, 5.9460678101]])
POS = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.5],
                [0.666666723, 0.333333334, 0.75], [0.333333334, 0.666666668, 0.25]])
NUM = np.array([25, 25, 52, 52])
MAG = np.array([+1.0, -1.0, 0.0, 0.0])
ALLOWED = (3, 9, 15, 21, 27, 33, 39, 45)
VBM, CBM, LMAX = 20, 21, 45


def _ops():
    return spin_operations(CELL, POS, NUM, MAG, symprec=1e-4)


def _fit(U, band, lmax=LMAX):
    q = read_bands(RUNS[U] / "data-file-schema.xml")
    th = np.arctan2(q.kpoints_tpiba[:, 1], q.kpoints_tpiba[:, 0])
    rad = np.linalg.norm(q.kpoints_tpiba[0][:2])
    kz = q.kpoints_tpiba[0, 2]
    dirs = np.stack([np.cos(th), np.sin(th), (kz / rad) * np.ones_like(th)], axis=1)
    basis = CircularBasis(lmax)
    f = project(basis, dirs, q.energies_ev[0, :, band] - q.energies_ev[1, :, band])
    floor = noise_report(f, set(allowed_l(basis, _ops(), dirs))).floor
    return f, floor, q, dirs


# ------------------------------------------------ the symmetry side, free

def test_mnte_has_no_symprec_cliff():
    """Where Cr2SeAs flips at 1e-5, MnTe is stable over seven decades."""
    rows = symprec_scan(CELL, POS, NUM, MAG,
                        symprecs=(1e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 1e-2))
    for r in rows:
        assert r["spacegroup"] == "P6_3/mmc"
        assert r["n_preserving"] == 12 and r["n_exchanging"] == 12


def test_the_allowed_set_is_the_textbook_g_wave_form():
    """l = 3 in plane, and identically zero on the k_z = 0 plane."""
    ops = _ops()
    b = 2 * np.pi * np.linalg.inv(CELL).T
    b1, b3 = np.linalg.norm(b[0]), np.linalg.norm(b[2])
    basis = CircularBasis(21)
    for frac in (0.1, 0.25, 0.4, 0.5):
        pts = circle(0.25 * b1, 200, center=(0, 0, frac * b3))
        assert sorted(allowed_l(basis, ops, pts)) == [3, 9, 15, 21], f"k_z = {frac} b3"
    pts = circle(0.25 * b1, 200, center=(0, 0, 0.0))
    assert allowed_l(basis, ops, pts) == set(), "the k_z = 0 plane must be nodal"


# ------------------------------------------------------------ the DFT side

@kz025
def test_the_splitting_is_essentially_a_pure_sine_of_three_theta():
    """rms should equal A_3 / sqrt(2) if nothing else contributes."""
    for U in (0.0, 4.0):
        for band in (VBM, CBM):
            f, _, q, _ = _fit(U, band)
            s = q.energies_ev[0, :, band] - q.energies_ev[1, :, band]
            a3 = f.amplitude_by_l()[3]
            assert np.sqrt((s ** 2).mean()) / (a3 / np.sqrt(2)) == pytest.approx(1.0, abs=0.01)


@kz025
def test_the_cosine_channel_is_far_cleaner_than_a_relaxed_structure():
    """Exact Wyckoff positions: no relaxation noise to break the symmetry.

    Compared at the SAME l, which a first version of this test did not do: the
    observed MnTe range is 1.2e-13 to 4.0e-12 over both U and both bands, and
    Cr2SeAs's l = 3 on its reference run is 8.6e-10.  The threshold below is set
    from that measured range rather than from one sample -- the first attempt
    used 1e-12, which the k_z != 0 re-run then exceeded at 1.06e-12.

    The ratio is independent of l_max (21 and 45 agree to every printed digit),
    so this is a property of the eigenvalues and not of the fit.
    """
    worst = 0.0
    for U in (0.0, 4.0):
        for band in (VBM, CBM):
            f, _, _, _ = _fit(U, band)
            c = dict(zip(f.labels, f.coeffs))
            worst = max(worst, abs(c[(3, "c")]) / abs(c[(3, "s")]))
    assert worst < 1e-11, f"worst cos/sin {worst:.2e}"
    assert worst > 1e-14, "if this ever gets below 1e-14 the comparison below is stale"

    cr = HERE / "data" / "qe_cr2seas_sm24" / "data-file-schema.xml"
    if cr.exists():
        q = read_bands(cr)
        th = np.arctan2(q.kpoints_tpiba[:, 1], q.kpoints_tpiba[:, 0])
        dirs = np.stack([np.cos(th), np.sin(th), 0.093 * np.ones_like(th)], axis=1)
        f = project(CircularBasis(45), dirs,
                    q.energies_ev[0, :, 44] - q.energies_ev[1, :, 44])
        c = dict(zip(f.labels, f.coeffs))
        relaxed = abs(c[(3, "c")]) / abs(c[(3, "s")])
        assert relaxed == pytest.approx(8.6e-10, rel=0.3)
        assert relaxed / worst > 100


@kz025
def test_a_semiconductor_has_no_fermi_surface_sampling_problem():
    """Floor 0.02-0.03 ueV, three orders of magnitude below Cr2SeAs's 24 ueV.

    Independent support for P3's diagnosis: that trap needed two bands within
    ~130 meV of E_F, and MnTe is gapped.
    """
    for band in (VBM, CBM):
        _, floor, _, _ = _fit(4.0, band)
        assert floor < 0.1e-6, f"band {band}: {floor*1e6:.3f} ueV"


@kz025
def test_all_eight_allowed_harmonics_are_resolved():
    """The result P8 exists for: the full ladder to l = 45, weakest at 4500x."""
    f, floor, _, _ = _fit(4.0, VBM)
    a = f.amplitude_by_l()
    margins = {l: a[l] / floor for l in ALLOWED}
    assert min(margins.values()) > 1000, f"weakest margin {min(margins.values()):.0f}x"
    assert margins[45] == pytest.approx(4493, rel=0.15)
    assert a[3] * 1e3 == pytest.approx(949.4, rel=0.01)
    assert a[45] * 1e3 == pytest.approx(0.105, rel=0.05)


@kz025
def test_the_valence_band_ladder_is_exponential_in_l():
    """A ∝ exp(-l/6.35), R^2 = 0.997, beating a power law."""
    f, _, _, _ = _fit(4.0, VBM)
    a = f.amplitude_by_l()
    L = np.array(ALLOWED, dtype=float)
    A = np.array([a[l] for l in ALLOWED])
    m = L >= 15                                   # the leading two are not asymptotic
    p = np.polyfit(L[m], np.log(A[m]), 1)
    r2 = lambda pred: 1 - np.sum((np.log(A[m]) - np.log(pred)) ** 2) / \
                          np.sum((np.log(A[m]) - np.log(A[m]).mean()) ** 2)
    assert -1 / p[0] == pytest.approx(6.35, rel=0.05)
    assert r2(np.exp(np.polyval(p, L[m]))) > 0.99
    q = np.polyfit(np.log(L[m]), np.log(A[m]), 1)
    assert r2(np.exp(np.polyval(q, np.log(L[m])))) < r2(np.exp(np.polyval(p, L[m])))


@kz025
def test_the_residual_matches_the_extrapolated_tail():
    """A check that does not go through the fit: what is left over above l = 45
    should be what the same exponential predicts for l = 51..75."""
    f, _, _, _ = _fit(4.0, VBM)
    a = f.amplitude_by_l()
    L = np.array(ALLOWED, dtype=float)
    A = np.array([a[l] for l in ALLOWED])
    m = L >= 15
    p = np.polyfit(L[m], np.log(A[m]), 1)
    tail = sum(np.exp(np.polyval(p, l)) for l in (51, 57, 63, 69, 75)) / np.sqrt(2)
    assert f.residual_rms == pytest.approx(tail, rel=0.25)


@kz025
def test_the_exponential_is_a_band_property_not_a_law():
    """Band 21 is not a clean exponential: one band does not make a law."""
    f, _, _, _ = _fit(4.0, CBM)
    a = f.amplitude_by_l()
    ratios = [a[ALLOWED[i]] / a[ALLOWED[i + 1]] for i in range(len(ALLOWED) - 1)]
    assert not all(ratios[i] >= ratios[i + 1] for i in range(len(ratios) - 1)), (
        "band 21's step ratios are non-monotonic, unlike band 20's"
    )
    L = np.array(ALLOWED, dtype=float)
    A = np.array([a[l] for l in ALLOWED])
    m = L >= 15
    p = np.polyfit(L[m], np.log(A[m]), 1)
    r2 = 1 - np.sum((np.log(A[m]) - np.polyval(p, L[m])) ** 2) / \
             np.sum((np.log(A[m]) - np.log(A[m]).mean()) ** 2)
    assert r2 < 0.99, "if this ever exceeds band 20's, revisit the claim"


@kz025
def test_u_matters_here_too_and_differently_per_band():
    """VBM falls 19% from U=0 to 4, CBM rises 146%."""
    a0v = _fit(0.0, VBM)[0].amplitude_by_l()[3]
    a4v = _fit(4.0, VBM)[0].amplitude_by_l()[3]
    a0c = _fit(0.0, CBM)[0].amplitude_by_l()[3]
    a4c = _fit(4.0, CBM)[0].amplitude_by_l()[3]
    assert (a4v - a0v) / a0v == pytest.approx(-0.19, abs=0.03)
    assert (a4c - a0c) / a0c == pytest.approx(+1.46, abs=0.10)


# --------------------------------------------- the nodal plane, and a mistake


@kz0
def test_the_kz0_plane_is_nodal_to_four_decimal_orders():
    """The other half of the prediction, and the one that is hard to fake.

    Symmetry forbids EVERY harmonic at k_z = 0, so the splitting there must
    vanish identically.  Measured on band 20: 19.25 ueV, against 833 meV on the
    same band at k_z = 0.25 b3 -- a ratio of 4.3e4, i.e. nodal to 0.0023%.

    CoF2 must come out the opposite way on the same plane (its d-wave form
    k_x k_y carries no k_z), which is why the two systems are run together:
    a projector with a systematic error would struggle to satisfy both.
    """
    q = read_bands(KZ0 / "data-file-schema.xml")
    assert abs(q.kpoints_tpiba[0, 2]) < 1e-12, "this fixture must be the k_z = 0 circle"
    s = q.energies_ev[0, :, VBM] - q.energies_ev[1, :, VBM]
    rms = float(np.sqrt((s ** 2).mean()))
    assert rms == pytest.approx(19.25e-6, rel=0.10)
    assert rms < 1e-4, "the k_z = 0 plane must be nodal"
    # against the same band away from the plane
    assert 833e-3 / rms > 2e4


@kz0
def test_the_nodal_plane_residual_is_not_the_forbidden_harmonic_floor():
    """An open question this run raised, recorded rather than smoothed over.

    At k_z != 0 the forbidden-harmonic floor on this band is 0.02-0.03 ueV.  At
    k_z = 0, where symmetry forbids everything, the total splitting is 19 ueV --
    a thousand times larger.  So the two are not the same quantity: the k_z != 0
    floor measures aliasing of the signal's own high-l tail (P3/P12), while the
    k_z = 0 residual is whatever survives when there is no signal to alias.
    **It is not a fixed noise level.** The same measurement at U = 4 gives
    0.20 ueV, a hundred times smaller, so whatever produces it is U-dependent --
    which points at near-degeneracy in the U = 0 band structure rather than at a
    solver noise floor.  At U = 4 the plane is nodal to 0.20 ueV against 950 meV
    off-plane, a ratio of 5e6.
    """
    q = read_bands(KZ0 / "data-file-schema.xml")
    rms = {b: float(np.sqrt(((q.energies_ev[0, :, b] - q.energies_ev[1, :, b]) ** 2).mean()))
           for b in (19, 20, 21, 22)}
    assert all(1e-5 < v < 1e-4 for v in rms.values()), rms
    # a thousand times the k_z != 0 floor, so it is a different effect
    assert min(rms.values()) > 100 * 0.03e-6


def test_the_two_circles_are_separate_files():
    """Guards the mistake that produced this split.

    Both circles are bands runs over one converged density, and the first
    version of the MnTe job wrote both into `tmp/MT.save/` and copied the
    fixture afterwards -- so the k_z = 0 run silently overwrote the k_z != 0
    eigenvalues, and the committed fixture was the wrong circle under the right
    name.  Three tests here failed on it, which is the only reason it was
    caught.  The CoF2 job had it right: copy out between runs.
    """
    assert KZ0.name.endswith("kz0")
    for d in RUNS.values():
        assert "kz025" in d.name
        if (d / "data-file-schema.xml").exists():
            q = read_bands(d / "data-file-schema.xml")
            assert abs(q.kpoints_tpiba[0, 2] - 0.177430) < 1e-4, f"{d.name} is not the k_z != 0 circle"


@pytest.mark.skipif(not (KZ0_U4 / "data-file-schema.xml").exists(),
                    reason="the U=4 nodal-plane fixture is not present")
def test_the_nodal_plane_is_cleaner_at_u4():
    """0.20 ueV against 19.25 at U = 0: the residual is U-dependent, not a floor."""
    q = read_bands(KZ0_U4 / "data-file-schema.xml")
    assert abs(q.kpoints_tpiba[0, 2]) < 1e-12
    s = q.energies_ev[0, :, VBM] - q.energies_ev[1, :, VBM]
    rms = float(np.sqrt((s ** 2).mean()))
    assert rms == pytest.approx(0.20e-6, rel=0.15)
    q0 = read_bands(KZ0 / "data-file-schema.xml")
    s0 = q0.energies_ev[0, :, VBM] - q0.energies_ev[1, :, VBM]
    assert float(np.sqrt((s0 ** 2).mean())) / rms > 50
    # nodal to 5e6 against the off-plane splitting
    assert 949e-3 / rms > 1e6
