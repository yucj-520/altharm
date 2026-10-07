"""P8: CoF2, the rutile d-wave, and the k_z = 0 test that MnTe must fail.

CoF2 is the third system and the one that makes the symmetry projector hard to
fool. Rutile's d-wave form is dE ∝ k_x k_y, which carries **no** k_z, so the
k_z = 0 plane is NOT nodal -- the opposite of MnTe and Cr2Se2As, where the
sublattice-exchanging mirror forces dE = -dE there. One projector, one analysis
pipeline, two systems, two opposite required answers.

**Measured** (U = 0, 10x10x14, degauss 0.0150, 152-point circle of radius
0.25|b1|, band 30 = valence-band maximum):

| | rms splitting on band 30 |
|---|---|
| CoF2 at k_z = 0 | **86.6 meV** |
| MnTe at k_z = 0 | **0.019 meV** |

A ratio of 4500, from symmetry alone with no DFT input on either side.

**And l = 2 is the right one of the two.** l = 2 has two functions, k_x k_y and
k_x^2 - k_y^2, and they are physically different. Rutile must give the first:
the x = y diagonal plane is a mirror through the Co-F bonds (F sits at (u,u,0))
that maps each Co to itself, so chi = +1, and under it k_x k_y is even while
k_x^2 - k_y^2 is odd. Measured cos/sin ratio at l = 2: **9.5e-15**, machine
precision. The symmetrizer predicted this before the DFT ran -- it keeps 100% of
sin2θ and 0.00% of cos2θ.

**CoF2 is nearly a single harmonic**, which turned out to be the interesting
part. Its ladder on band 30 at k_z = 0 is 122.5 / 0.021 / 0.0014 / 0.00004 meV
for l = 2, 6, 10, 14: the second allowed harmonic is **0.021%** of the first.
Against MnTe's 5.6% and Cr2Se2As's 67.6%. So of the three, the two canonical
altermagnets are essentially single-wave and **the group's own Cr2Se2As is the
one with genuinely multi-harmonic structure** -- l = 9 at two thirds of l = 3.
That is where the "high harmonic part would be genuinely new" lives.

The comparison survives a caveat that runs the other way: CoF2's allowed set
steps by 4 (2, 6, 10, 14) while the other two step by 6, and a shorter step
should decay *less*, so 0.021% over Δl = 4 is more extreme than it looks.
"""

from pathlib import Path

import numpy as np
import pytest

from altharm.basis import CircularBasis
from altharm.kpath import circle
from altharm.noise import noise_report
from altharm.project import project
from altharm.qe.reader import read_bands
from altharm.symmetry import allowed_l, spin_operations, symmetrizer, symprec_scan

HERE = Path(__file__).parent
KZ0 = HERE / "data" / "cof2_u0_kz0"
KZ025 = HERE / "data" / "cof2_u0_kz025"
MNTE_KZ0 = HERE / "data" / "mnte_kz0"

pytestmark = pytest.mark.skipif(
    not all((d / "data-file-schema.xml").exists() for d in (KZ0, KZ025)),
    reason="the CoF2 fixtures are not present",
)

CELL = np.array([[4.7467970848, 0.0, 0.0],
                 [0.0, 4.7467970848, 0.0],
                 [0.0, 0.0, 3.2185559273]])
POS = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5],
                [0.695713003, 0.695713003, 0.0], [0.195713003, 0.804286997, 0.5],
                [0.304286997, 0.304286997, 0.0], [0.804286997, 0.195713003, 0.5]])
NUM = np.array([27, 27, 9, 9, 9, 9])
MAG = np.array([+1.0, -1.0, 0.0, 0.0, 0.0, 0.0])
ALLOWED = (2, 6, 10, 14)
VBM, LMAX = 30, 14


def _ops():
    return spin_operations(CELL, POS, NUM, MAG, symprec=1e-4)


def _fit(path, band=VBM):
    q = read_bands(path / "data-file-schema.xml")
    th = np.arctan2(q.kpoints_tpiba[:, 1], q.kpoints_tpiba[:, 0])
    rad = np.linalg.norm(q.kpoints_tpiba[0][:2])
    kz = q.kpoints_tpiba[0, 2]
    dirs = np.stack([np.cos(th), np.sin(th), (kz / rad) * np.ones_like(th)], axis=1)
    basis = CircularBasis(LMAX)
    f = project(basis, dirs, q.energies_ev[0, :, band] - q.energies_ev[1, :, band])
    floor = noise_report(f, set(allowed_l(basis, _ops(), dirs))).floor
    return f, floor, q


def test_cof2_has_no_symprec_cliff():
    rows = symprec_scan(CELL, POS, NUM, MAG, symprecs=(1e-6, 1e-5, 1e-4, 1e-3, 1e-2))
    for r in rows:
        assert r["spacegroup"] == "P4_2/mnm"
        assert r["n_preserving"] == 8 and r["n_exchanging"] == 8


def test_the_projector_picks_kxky_over_kx2_minus_ky2_before_any_dft():
    """The discriminating prediction, made from symmetry alone."""
    b = 2 * np.pi * np.linalg.inv(CELL).T
    pts = circle(0.25 * np.linalg.norm(b[0]), 400, center=(0, 0, 0.0))
    basis = CircularBasis(6)
    S, kept = symmetrizer(basis, _ops(), pts)
    th = np.arctan2(pts[:, 1], pts[:, 0])
    keep = lambda v: float(np.linalg.norm(S @ np.asarray(project(basis, pts, v).coeffs))
                           / np.linalg.norm(np.asarray(project(basis, pts, v).coeffs)))
    assert keep(np.sin(2 * th)) > 0.999, "k_x k_y must survive"
    assert keep(np.cos(2 * th)) < 0.001, "k_x^2 - k_y^2 must be projected out"


def test_the_kz0_plane_is_not_nodal_unlike_mnte():
    """The test a systematically wrong projector would fail on one system or the other."""
    f, _, q = _fit(KZ0)
    assert abs(q.kpoints_tpiba[0, 2]) < 1e-12
    s = q.energies_ev[0, :, VBM] - q.energies_ev[1, :, VBM]
    rms_cof2 = float(np.sqrt((s ** 2).mean()))
    assert rms_cof2 == pytest.approx(86.6e-3, rel=0.05)
    assert allowed_l(CircularBasis(LMAX), _ops(),
                     circle(0.25 * np.linalg.norm(2 * np.pi * np.linalg.inv(CELL).T[0]),
                            400, center=(0, 0, 0.0))) == {2, 6, 10, 14}

    if (MNTE_KZ0 / "data-file-schema.xml").exists():
        qm = read_bands(MNTE_KZ0 / "data-file-schema.xml")
        sm = qm.energies_ev[0, :, 20] - qm.energies_ev[1, :, 20]
        rms_mnte = float(np.sqrt((sm ** 2).mean()))
        assert rms_cof2 / rms_mnte > 2000, (
            f"CoF2 {rms_cof2*1e3:.1f} meV vs MnTe {rms_mnte*1e6:.1f} ueV on the same plane"
        )


def test_l2_is_pure_kxky_to_machine_precision():
    f, _, _ = _fit(KZ0)
    c = dict(zip(f.labels, f.coeffs))
    assert abs(c[(2, "c")]) / abs(c[(2, "s")]) < 1e-13


def test_the_splitting_is_essentially_a_pure_sine_of_two_theta():
    for path in (KZ0, KZ025):
        f, _, q = _fit(path)
        s = q.energies_ev[0, :, VBM] - q.energies_ev[1, :, VBM]
        a2 = f.amplitude_by_l()[2]
        assert np.sqrt((s ** 2).mean()) / (a2 / np.sqrt(2)) == pytest.approx(1.0, abs=0.01)


def test_cof2_is_nearly_a_single_harmonic():
    """0.021% -- the contrast that makes Cr2Se2As the interesting one."""
    f, floor, _ = _fit(KZ0)
    a = f.amplitude_by_l()
    assert a[2] * 1e3 == pytest.approx(122.5, rel=0.02)
    assert a[6] / a[2] == pytest.approx(1.7e-4, rel=0.3)
    assert a[6] / a[2] < 1e-3
    # and it is still above the floor, so this is a measurement not a limit
    assert a[6] > 100 * floor


def test_high_harmonic_richness_spans_three_orders_across_the_three_systems():
    """The headline of P8, and it did not go the way one would guess.

    next-allowed over leading, on the band nearest E_F of each:
      CoF2    2.1e-4   (l=6 / l=2)
      MnTe    5.6e-2   (l=9 / l=3)
      Cr2Se2As 6.8e-1   (l=9 / l=3)

    The two canonical altermagnets are essentially single-wave; the group's own
    material is the multi-harmonic one.  Note CoF2's set steps by 4 rather than
    6, and a shorter step should decay LESS, so its ratio is if anything
    understated relative to the others.
    """
    f, _, _ = _fit(KZ025)
    a = f.amplitude_by_l()
    cof2 = a[6] / a[2]

    cr = HERE / "data" / "qe_cr2seas_sm24" / "data-file-schema.xml"
    if not cr.exists():
        pytest.skip("the Cr2Se2As reference fixture is not present")
    q = read_bands(cr)
    th = np.arctan2(q.kpoints_tpiba[:, 1], q.kpoints_tpiba[:, 0])
    dirs = np.stack([np.cos(th), np.sin(th), 0.093 * np.ones_like(th)], axis=1)
    ac = project(CircularBasis(45), dirs,
                 q.energies_ev[0, :, 44] - q.energies_ev[1, :, 44]).amplitude_by_l()
    cr2seas = ac[9] / ac[3]

    assert cof2 < 1e-3
    assert cr2seas == pytest.approx(0.676, rel=0.05)
    assert cr2seas / cof2 > 1000, "Cr2Se2As should be far richer in high harmonics"


# ------------------------- U = 4 reorders the bands, so the index stops meaning
#                           the same thing -- a failure mode, recorded as one

U4_KZ0 = HERE / "data" / "cof2_u4_kz0"
U4_KZ025 = HERE / "data" / "cof2_u4_kz025"

u4 = pytest.mark.skipif(
    not all((d / "data-file-schema.xml").exists() for d in (U4_KZ0, U4_KZ025)),
    reason="the CoF2 U=4 fixtures are not present",
)


@u4
def test_u4_reorders_the_bands_so_band_30_is_not_the_same_state():
    """Why this file quotes CoF2 at U = 0 and not as a U-robustness pair.

    At U = 0, band 30 is the valence-band maximum at -0.256 eV and band 31 the
    conduction minimum at +0.283 eV: a 540 meV gap.  At U = 4 both sit within
    21 meV of E_F and are **4.07 meV apart**.  Comparing coefficient-by-index
    across that is meaningless, and it shows up as an apparent collapse of
    l = 2 from 46.7 to 6.9 meV with the noise floor jumping 150x.

    It is not a failed calculation: both runs converged, the compensation holds
    (total magnetization 0), the moments agree (5.65 vs 5.72 uB) and the Hubbard
    occupations are symmetric between the two Co (4.987/2.232 and its mirror).
    U really does push these two bands together here.

    So this is the band-identity trap of `CLAUDE.md` in a form the project had
    not guarded: not across two models with different band counts, but **across
    two values of U in the same code**.  A U scan cannot be read off a fixed
    band index without checking that the index still labels the same state.
    """
    e, gap = {}, {}
    for tag, d in (("u0", KZ025), ("u4", U4_KZ025)):
        q = read_bands(d / "data-file-schema.xml")
        e[tag] = (q.energies_ev[0, :, 30].mean() - q.fermi_ev,
                  q.energies_ev[0, :, 31].mean() - q.fermi_ev)
        gap[tag] = float((q.energies_ev[0, :, 31] - q.energies_ev[0, :, 30]).min())

    assert e["u0"][0] == pytest.approx(-0.256, abs=0.01)
    assert e["u0"][1] == pytest.approx(+0.283, abs=0.01)
    assert gap["u0"] == pytest.approx(0.517, abs=0.03)

    assert abs(e["u4"][0]) < 0.01 and abs(e["u4"][1]) < 0.03
    assert gap["u4"] == pytest.approx(0.00407, abs=0.001)
    assert gap["u0"] / gap["u4"] > 50, "the gap collapses by two orders of magnitude"


@u4
def test_the_symmetry_survives_the_reordering_even_though_the_index_does_not():
    """What is still true at U = 4: l = 2 is pure k_x k_y on whatever band it is.

    The selection rule is a statement about the crystal, not about a band, so it
    has to hold regardless of which state the index happens to name.
    """
    for d in (U4_KZ0, U4_KZ025):
        for band in (30, 31):
            q = read_bands(d / "data-file-schema.xml")
            th = np.arctan2(q.kpoints_tpiba[:, 1], q.kpoints_tpiba[:, 0])
            rad = np.linalg.norm(q.kpoints_tpiba[0][:2])
            kz = q.kpoints_tpiba[0, 2]
            dirs = np.stack([np.cos(th), np.sin(th), (kz / rad) * np.ones_like(th)], axis=1)
            f = project(CircularBasis(LMAX), dirs,
                        q.energies_ev[0, :, band] - q.energies_ev[1, :, band])
            c = dict(zip(f.labels, f.coeffs))
            assert abs(c[(2, "c")]) / abs(c[(2, "s")]) < 1e-11, f"{d.name}, band {band}"


@u4
def test_the_kz0_plane_stays_non_nodal_at_u4():
    """The discriminating prediction is U-independent, as a symmetry claim must be."""
    q = read_bands(U4_KZ0 / "data-file-schema.xml")
    assert abs(q.kpoints_tpiba[0, 2]) < 1e-12
    rms = [float(np.sqrt(((q.energies_ev[0, :, b] - q.energies_ev[1, :, b]) ** 2).mean()))
           for b in (29, 30, 31, 32)]
    assert max(rms) > 1e-3, "CoF2's k_z = 0 plane must carry splitting at any U"
