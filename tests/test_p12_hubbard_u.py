"""The Hubbard U sensitivity: the largest uncertainty in the whole pipeline.

U = 3.0 eV is an inherited choice.  Nothing in this project had tested what it
costs, and the forbidden-harmonic floor cannot see it -- U is an input to the
Hamiltonian, and the floor measures the fit.  This is the third and worst member
of the family P3 (smearing) and the geometry traps opened.

Four SCF + 152-point bands runs at U = 2, 3, 4, 5 eV, all at the corrected
settings (12x12x4, degauss 0.0150), band 44:

    l        U=2      U=3      U=4      U=5      2 -> 5
    3     69.152   63.746   54.660   44.260     -36%
    9     45.094   42.921   36.080   28.892     -36%
    15     0.000    3.533    5.358    5.940    switched on
    21     9.550    8.644    6.889    5.141     -46%
    27     1.098    1.837    1.984    1.762     +61%
    33     2.962    2.614    2.069    1.457     -51%
    39     1.244    1.168    0.889    0.602     -52%
    45     1.235    1.003    0.761    0.482     -61%

**About 21% per eV at the median and 24% at worst** (excluding l=15, which is near
zero at U=2 and on its own moves 52% from U=3 to U=4) -- against 12.9% for the
smearing error P3 found, under 4% for a 1e-3 A geometry perturbation, and 2.4%
of l=45 for the noise floor itself.  U dominates every other uncertainty here.

**And it changes the shape, not the scale.**  Six harmonics shrink with U while
l=27 rises and falls and l=15 grows from nothing, so no rescaling maps one U's
spectrum onto another's.  l=15 is the sharpest case: at U=2 it is 0.3 ueV, which
is 250x *below* that run's own noise floor -- the harmonic is absent, not small.
Its coefficient is positive at U=3,4,5 and approaching zero from above at U=2,
so it crosses zero near U ~ 2.  A selection rule says which harmonics *may* be
non-zero; which ones *are* is a question about the interaction strength.

**The floor tracks U too**, 81.8 / 24.3 / 6.5 / 1.3 ueV, a factor of 65 across
the range, while band 44's distance from E_F grows 0.017 -> 0.066 eV.  That is
P3's Fermi-surface-sampling mechanism seen from the other side: the closer the
band sits to E_F, the worse the sampling, and U is what moves it.  **The quoted
noise floor is therefore a property of (settings, U), not of the material.**

**Consequence.**  Every coefficient this project reports is a coefficient *at
U = 3.0 eV*, and that has to be said.  The principled fix is not a scan but a
first-principles U: QE's `hp.x` computes it by linear response, which is the
obvious next step on this route.
"""

from pathlib import Path

import numpy as np
import pytest

from altharm.basis import CircularBasis
from altharm.noise import noise_report
from altharm.project import project
from altharm.qe.reader import read_bands
from altharm.symmetry import SpinOp, allowed_l

HERE = Path(__file__).parent
RUNS = {2.0: "qe_cr2seas_u2", 3.0: "qe_cr2seas_sm12",
        4.0: "qe_cr2seas_u4", 5.0: "qe_cr2seas_u5"}
REPORT = HERE / "data" / "harmwave_cr2seas" / "harmwave_report.txt"

pytestmark = pytest.mark.skipif(
    not (all((HERE / "data" / d / "data-file-schema.xml").exists() for d in RUNS.values())
         and REPORT.exists()),
    reason="the Hubbard-U scan fixtures are not present",
)

BAND, LMAX = 44, 45
ALLOWED = (3, 9, 15, 21, 27, 33, 39, 45)


def _ops():
    ops = []
    for line in REPORT.read_text().splitlines():
        if line.startswith("OP"):
            p = line.split()
            ops.append(SpinOp(rot_frac=None, trans=None, chi=int(p[1]),
                              rot_cart=np.array([float(x) for x in p[2:11]]).reshape(3, 3)))
    return ops


@pytest.fixture(scope="module")
def scan():
    out = {}
    for U, folder in RUNS.items():
        q = read_bands(HERE / "data" / folder / "data-file-schema.xml")
        th = np.arctan2(q.kpoints_tpiba[:, 1], q.kpoints_tpiba[:, 0])
        dirs = np.stack([np.cos(th), np.sin(th), 0.093 * np.ones_like(th)], axis=1)
        basis = CircularBasis(LMAX)
        f = project(basis, dirs, q.energies_ev[0, :, BAND] - q.energies_ev[1, :, BAND])
        floor = noise_report(f, allowed_l(basis, _ops(), dirs)).floor
        out[U] = (f.amplitude_by_l(), floor, q, dict(zip(f.labels, f.coeffs)))
    return out


def test_only_u_differs_between_the_runs():
    """Guard the premise."""
    for U, folder in RUNS.items():
        text = (HERE / "data" / folder / "scf.in").read_text()
        assert "degauss     = 0.0150" in text and "12 12 4" in text
        assert "ecutwfc     = 60" in text and "vdw_corr" not in text
        assert f"U CrU-3d {U:g}" in text or f"U CrU-3d {U:.1f}" in text


def test_u_is_the_largest_uncertainty_measured_anywhere_here(scan):
    """21% per eV at the median, 24% worst -- above smearing, geometry, floor.

    l = 15 is excluded from the spread: it is near zero at U = 2, so a relative
    change there is meaningless.  Counted in, it alone moves 52% from U = 3 to
    U = 4, which is the largest single number in the scan but says more about
    its zero crossing than about the typical sensitivity.
    """
    a3, a4 = scan[3.0][0], scan[4.0][0]
    per_ev = {l: abs(a4[l] - a3[l]) / a3[l] for l in ALLOWED if l != 15}
    assert np.median(list(per_ev.values())) == pytest.approx(0.21, abs=0.05)
    assert max(per_ev.values()) == pytest.approx(0.24, abs=0.03)
    # l=15 on its own, the outlier the exclusion is about
    assert abs(a4[15] - a3[15]) / a3[15] == pytest.approx(0.52, abs=0.05)
    # bigger than the smearing error P3 found (12.9% median)
    assert np.median(list(per_ev.values())) > 0.129


def test_u_changes_the_shape_not_the_scale(scan):
    """No rescaling maps one U's spectrum onto another's."""
    a2, a5 = scan[2.0][0], scan[5.0][0]
    ratios = np.array([a5[l] / a2[l] for l in ALLOWED if l != 15])
    assert ratios.max() / ratios.min() > 1.5, "a uniform rescale would make these equal"
    # six shrink, l=27 does not
    assert a5[3] < a2[3] and a5[9] < a2[9] and a5[45] < a2[45]
    assert a5[27] > a2[27]


def test_l15_is_switched_on_by_u(scan):
    """At U=2 it is 250x BELOW that run's own floor: absent, not small."""
    a2, floor2, _, c2 = scan[2.0]
    assert a2[15] == pytest.approx(0.3e-6, abs=0.3e-6)
    assert a2[15] < floor2 / 100
    for U in (3.0, 4.0, 5.0):
        a, floor, _, _ = scan[U]
        assert a[15] > 20 * floor
    # positive at 3,4,5 and approaching zero from above at 2: it crosses near U~2
    for U in (2.0, 3.0, 4.0, 5.0):
        assert scan[U][3][(15, "s")] > 0
    assert scan[2.0][0][15] < scan[3.0][0][15] < scan[4.0][0][15] < scan[5.0][0][15]


def test_every_allowed_harmonic_is_pure_sine_at_every_u(scan):
    """The symmetry holds regardless of U -- a control on the whole scan."""
    for U in RUNS:
        a, floor, _, c = scan[U]
        for l in ALLOWED:
            assert abs(c[(l, "c")]) < max(3 * floor, 0.02 * abs(c[(l, "s")])), f"U={U}, l={l}"


def test_the_noise_floor_is_itself_a_function_of_u(scan):
    """81.8 -> 1.3 ueV across U = 2..5."""
    floors = {U: scan[U][1] for U in RUNS}
    assert floors[2.0] == pytest.approx(81.8e-6, rel=0.05)
    assert floors[5.0] == pytest.approx(1.3e-6, rel=0.15)
    assert floors[2.0] / floors[5.0] > 40
    assert floors[2.0] > floors[3.0] > floors[4.0] > floors[5.0]


def test_the_floor_tracks_the_unresolved_tail_not_the_distance_from_ef():
    """The correction U = 0 forced.

    On U = 2..5 alone the floor falls monotonically as band 44 moves away from
    E_F, and this file first read that as the Fermi-surface-sampling mechanism
    of P3 seen from the other side.  Plain PBE breaks it: band 44 sits at
    +0.0332 eV, between U=2's +0.0166 and U=3's +0.0357, yet the floor is
    0.92 ueV -- the smallest of all five, 90x below U=2's.

    What it does track is the fit residual, i.e. the part of the signal living
    above l = 45.  That is what P3 identified as the thing that aliases in, so
    the corrected statement is stronger than the original: a five-point
    correlation of +0.91 against the residual, versus -0.77 against the gap.
    """
    import numpy as np
    from altharm.basis import CircularBasis
    from altharm.noise import noise_report
    from altharm.project import project
    from altharm.qe.reader import read_bands

    folders = {0.0: "qe_cr2seas_u0", 2.0: "qe_cr2seas_u2", 3.0: "qe_cr2seas_sm12",
               4.0: "qe_cr2seas_u4", 5.0: "qe_cr2seas_u5"}
    if not (HERE / "data" / folders[0.0] / "data-file-schema.xml").exists():
        pytest.skip("the plain-PBE fixture is not present")

    gaps, floors, resid = [], [], []
    for U, folder in folders.items():
        q = read_bands(HERE / "data" / folder / "data-file-schema.xml")
        th = np.arctan2(q.kpoints_tpiba[:, 1], q.kpoints_tpiba[:, 0])
        dirs = np.stack([np.cos(th), np.sin(th), 0.093 * np.ones_like(th)], axis=1)
        basis = CircularBasis(LMAX)
        f = project(basis, dirs, q.energies_ev[0, :, BAND] - q.energies_ev[1, :, BAND])
        gaps.append(q.energies_ev[0, :, BAND].mean() - q.fermi_ev)
        floors.append(noise_report(f, allowed_l(basis, _ops(), dirs)).floor)
        resid.append(f.residual_rms)

    # plain PBE: middling gap, smallest floor, smallest tail
    assert gaps[0] == pytest.approx(0.0332, abs=0.002)
    assert gaps[1] < gaps[0] < gaps[2], "the gap is NOT monotonic in U once U=0 is in"
    assert floors[0] == min(floors) and resid[0] == min(resid)
    assert floors[1] / floors[0] > 50

    c = lambda x, y: float(np.corrcoef(x, y)[0, 1])
    assert c(resid, floors) > 0.85, "the floor should track the unresolved tail"
    assert c(gaps, floors) > -0.85, "and not the distance from E_F"


def test_plain_pbe_is_a_different_physics_not_a_neutral_baseline():
    """Removing U is another arbitrary choice, and a qualitatively different one.

    l=9 drops by 2.8x, l=21 by 4.5x, and l=15 comes out with the OPPOSITE sign
    (-5.540 meV against +3.533 at U=3), confirming the zero crossing near U ~ 2
    from the other side.  l=3 is not the U=2 trend continued either: 64.2 meV at
    U=0 against 69.2 at U=2, so the behaviour between 0 and 2 is not monotonic.
    """
    import numpy as np
    from altharm.basis import CircularBasis
    from altharm.project import project
    from altharm.qe.reader import read_bands

    if not (HERE / "data" / "qe_cr2seas_u0" / "data-file-schema.xml").exists():
        pytest.skip("the plain-PBE fixture is not present")
    out = {}
    for tag, folder in (("u0", "qe_cr2seas_u0"), ("u2", "qe_cr2seas_u2"),
                        ("u3", "qe_cr2seas_sm12")):
        q = read_bands(HERE / "data" / folder / "data-file-schema.xml")
        th = np.arctan2(q.kpoints_tpiba[:, 1], q.kpoints_tpiba[:, 0])
        dirs = np.stack([np.cos(th), np.sin(th), 0.093 * np.ones_like(th)], axis=1)
        f = project(CircularBasis(LMAX), dirs,
                    q.energies_ev[0, :, BAND] - q.energies_ev[1, :, BAND])
        out[tag] = (f.amplitude_by_l(), dict(zip(f.labels, f.coeffs)))

    assert out["u3"][0][9] / out["u0"][0][9] == pytest.approx(2.8, abs=0.3)
    assert out["u3"][0][21] / out["u0"][0][21] == pytest.approx(4.5, abs=0.5)
    # the sign flip
    assert out["u0"][1][(15, "s")] < 0 < out["u3"][1][(15, "s")]
    assert out["u0"][1][(15, "s")] * 1e3 == pytest.approx(-5.54, abs=0.2)
    # and l=3 is not monotonic from 0 to 2
    assert out["u0"][0][3] < out["u2"][0][3]
