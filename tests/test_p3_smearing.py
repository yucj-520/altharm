"""P3, resolved: the non-convergence was Fermi-surface sampling, and the fix works.

Bands 43 and 44 straddle E_F (-0.13 and +0.04 eV) while the production Gaussian
smearing is only 0.0037 Ry = 50 meV.  Their occupations therefore track how
finely the Fermi surface is sampled, and the density -- converged by moment,
Hubbard occupation and total energy -- keeps moving the valence eigenvalues.

Widening the smearing to 0.015 Ry makes the k-convergence fast.  Two knobs, four
runs, and the point is that three of the four land on the same answer:

    degauss / mesh      l=3 (meV)
    0.0037 / 12x12x4      59.17     <- the production setting
    0.0150 / 12x12x4      63.75
    0.0037 / 24x24x8      64.36
    0.0150 / 24x24x8      64.29

The last three agree to ~1%; the production one is 8% below them.  Since the two
24x24x8 runs differ by only 0.1% although their smearing differs 4x, the
smearing's own bias is negligible there and 64.3 meV is the converged value.

**Consequence for the group's setup:** 12x12x4 at degauss 0.0037 underestimates
l=3 by 8% and the weak harmonics by up to 27%.  The cheap fix is degauss = 0.015
at the same mesh, which lands within 1-2% of converged.
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
RUNS = {
    ("small", 12): HERE / "data" / "qe_cr2seas_dense" / "data-file-schema.xml",
    ("small", 24): HERE / "data" / "qe_cr2seas_scf2488" / "data-file-schema.xml",
    ("wide", 12): HERE / "data" / "qe_cr2seas_sm12" / "data-file-schema.xml",
    ("wide", 24): HERE / "data" / "qe_cr2seas_sm24" / "data-file-schema.xml",
}
REPORT = HERE / "data" / "harmwave_cr2seas" / "harmwave_report.txt"

pytestmark = pytest.mark.skipif(
    not (all(p.exists() for p in RUNS.values()) and REPORT.exists()),
    reason="the smearing-scan fixtures are not present",
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
def runs():
    out = {}
    for key, xml in RUNS.items():
        q = read_bands(xml)
        th = np.arctan2(q.kpoints_tpiba[:, 1], q.kpoints_tpiba[:, 0])
        dirs = np.stack([np.cos(th), np.sin(th), 0.093 * np.ones_like(th)], axis=1)
        s = q.energies_ev[0, :, BAND] - q.energies_ev[1, :, BAND]
        fit = project(CircularBasis(LMAX), dirs, s)
        floor = noise_report(fit, allowed_l(CircularBasis(LMAX), _ops(), dirs)).floor
        out[key] = (fit.amplitude_by_l(), s, floor, q)
    return out


def _rms(x):
    return float(np.sqrt((np.asarray(x) ** 2).mean()))


def test_the_inputs_differ_only_in_smearing_and_mesh(runs):
    """Guard the premise."""
    for tag, folder in (("wide", "qe_cr2seas_sm12"), ("wide", "qe_cr2seas_sm24")):
        text = (HERE / "data" / folder / "scf.in").read_text()
        assert "degauss     = 0.0150" in text
        assert "ecutwfc     = 60" in text and "U CrU-3d 3.0" in text
        assert "vdw_corr" not in text


def test_wider_smearing_makes_the_mesh_converge(runs):
    """The headline: a 9x smaller change per mesh doubling."""
    small = _rms(runs[("small", 24)][1] - runs[("small", 12)][1])
    wide = _rms(runs[("wide", 24)][1] - runs[("wide", 12)][1])
    assert small == pytest.approx(4.90e-3, rel=0.05)
    assert wide == pytest.approx(0.541e-3, rel=0.08)
    assert small / wide > 6, "widening the smearing should fix the k-convergence"


def test_the_worst_harmonic_stops_moving(runs):
    """26.5% -> 2.2% on the most fragile coefficient."""
    def worst(tag):
        a12, a24 = runs[(tag, 12)][0], runs[(tag, 24)][0]
        return max(abs(a24[l] - a12[l]) / a12[l] for l in ALLOWED)
    assert worst("small") == pytest.approx(0.265, rel=0.10)
    assert worst("wide") < 0.035
    assert worst("small") / worst("wide") > 7


def test_three_settings_agree_so_the_converged_value_is_known(runs):
    """Two independent knobs landing together is what makes 64.3 meV the answer."""
    a = {k: runs[k][0] for k in RUNS}
    # the two 24x24x8 runs differ 4x in smearing yet agree: smearing bias is small there
    for l in (3, 9, 21, 27):
        rel = abs(a[("wide", 24)][l] - a[("small", 24)][l]) / a[("small", 24)][l]
        assert rel < 0.015, f"l={l} smearing bias {rel:.3f} at the fine mesh"
    # and the wide/12 run is already close to both
    assert a[("wide", 12)][3] == pytest.approx(a[("small", 24)][3], rel=0.02)
    assert a[("small", 24)][3] == pytest.approx(64.3e-3, rel=0.01)


def test_the_production_setting_is_the_outlier(runs):
    """What this costs the group: 8% on l=3, more on the weak harmonics."""
    prod = runs[("small", 12)][0]
    conv = runs[("wide", 24)][0]
    assert abs(conv[3] - prod[3]) / prod[3] == pytest.approx(0.087, abs=0.02)
    worst = max(abs(conv[l] - prod[l]) / prod[l] for l in ALLOWED)
    assert worst > 0.20, "the weak harmonics are off by more than 20%"


def test_the_noise_floor_barely_moves(runs):
    """A control: none of this is the floor changing under us."""
    floors = [runs[k][2] for k in RUNS]
    assert all(20e-6 < f < 30e-6 for f in floors)
    assert max(floors) / min(floors) < 1.25
