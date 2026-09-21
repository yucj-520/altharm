"""The corrected reference values, and what the old settings cost.

P3's smearing result invalidated the absolute coefficients everywhere else in
this project: every number in P3, P6 and P7 was computed at degauss = 0.0037,
which under-samples the Fermi surface (bands 43 and 44 straddle E_F) and
under-converges the harmonics.  This file records the corrected set and pins the
size of the correction, so no later phase quietly picks up an old number.

**The reference is 24x24x8 at degauss = 0.0150.**  The production-cost stand-in
(12x12x4 at the same smearing) is within 1.2% median / 2.1% worst of it, which
is what makes the correction affordable.

The third column is the part that matters for the group's published work: the
Wannier model was fitted to DFT at the OLD settings, so it inherits the same
error.  P6's "<1% agreement at all eight harmonics" was therefore agreement on a
*shared* error, not evidence of absolute accuracy -- a distinction that only
became visible once the DFT side was corrected independently.
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
OLD = HERE / "data" / "qe_cr2seas_dense" / "data-file-schema.xml"      # 12x12x4, 0.0037
NEW12 = HERE / "data" / "qe_cr2seas_sm12" / "data-file-schema.xml"     # 12x12x4, 0.0150
NEW24 = HERE / "data" / "qe_cr2seas_sm24" / "data-file-schema.xml"     # 24x24x8, 0.0150
REPORT = HERE / "data" / "harmwave_cr2seas" / "harmwave_report.txt"

pytestmark = pytest.mark.skipif(
    not (OLD.exists() and NEW12.exists() and NEW24.exists() and REPORT.exists()),
    reason="the corrected-settings fixtures are not present",
)

BAND, LMAX = 44, 45
ALLOWED = (3, 9, 15, 21, 27, 33, 39, 45)

# The group's Wannier model, fitted to DFT at the OLD settings.
WANNIER_MEV = {3: 59.4812, 9: 39.0640, 15: 3.1750, 21: 8.4434,
               27: 1.4946, 33: 2.8382, 39: 0.9964, 45: 1.1423}

# THE REFERENCE VALUES: 24x24x8, degauss 0.0150, band 44, 152-point circle.
REFERENCE_MEV = {3: 64.2928, 9: 43.4467, 15: 3.6033, 21: 8.6617,
                 27: 1.8763, 33: 2.5832, 39: 1.1819, 45: 0.9853}


def _ops():
    ops = []
    for line in REPORT.read_text().splitlines():
        if line.startswith("OP"):
            p = line.split()
            ops.append(SpinOp(rot_frac=None, trans=None, chi=int(p[1]),
                              rot_cart=np.array([float(x) for x in p[2:11]]).reshape(3, 3)))
    return ops


def _fit(xml):
    q = read_bands(xml)
    th = np.arctan2(q.kpoints_tpiba[:, 1], q.kpoints_tpiba[:, 0])
    dirs = np.stack([np.cos(th), np.sin(th), 0.093 * np.ones_like(th)], axis=1)
    basis = CircularBasis(LMAX)
    fit = project(basis, dirs, q.energies_ev[0, :, BAND] - q.energies_ev[1, :, BAND])
    floor = noise_report(fit, allowed_l(basis, _ops(), dirs)).floor
    return fit.amplitude_by_l(), floor, q


def test_band_44_is_still_the_right_band_under_the_new_smearing():
    """Guard first: a different smearing moves E_F, so re-check the band."""
    _, _, q = _fit(NEW24)
    assert q.fermi_ev == pytest.approx(8.1186, abs=1e-3)
    assert q.energies_ev[0, :, BAND].mean() - q.fermi_ev == pytest.approx(0.0373, abs=0.002)


def test_the_reference_values():
    """The numbers later work should quote."""
    a, _, _ = _fit(NEW24)
    for l in ALLOWED:
        assert a[l] * 1e3 == pytest.approx(REFERENCE_MEV[l], rel=0.002), f"l={l}"


def test_the_cheap_setting_is_close_enough_to_be_usable():
    """12x12x4 at degauss 0.015 costs a normal run and lands within ~2%."""
    a12, _, _ = _fit(NEW12)
    dev = {l: abs(a12[l] * 1e3 - REFERENCE_MEV[l]) / REFERENCE_MEV[l] for l in ALLOWED}
    assert max(dev.values()) < 0.03, f"worst l={max(dev, key=dev.get)}: {max(dev.values()):.3f}"
    assert np.median(list(dev.values())) < 0.02


def test_the_size_of_the_correction():
    """What the old settings cost: median 13%, worst 26%."""
    old, _, _ = _fit(OLD)
    shift = {l: abs(REFERENCE_MEV[l] - old[l] * 1e3) / (old[l] * 1e3) for l in ALLOWED}
    assert np.median(list(shift.values())) == pytest.approx(0.129, abs=0.02)
    assert max(shift.values()) == pytest.approx(0.263, abs=0.02)
    assert max(shift, key=shift.get) == 27
    # and it is not a uniform rescaling -- two harmonics move the other way
    signed = [REFERENCE_MEV[l] - old[l] * 1e3 for l in ALLOWED]
    assert any(s > 0 for s in signed) and any(s < 0 for s in signed)


def test_the_wannier_model_inherits_the_same_error():
    """P6's <1% agreement was agreement on a shared error.

    The Wannier model was fitted to DFT at the old settings, so it sits with the
    old DFT numbers, not with the corrected ones.  This is the single most
    important consequence of P3's smearing result for the group's published work.
    """
    old, _, _ = _fit(OLD)
    # Wannier still agrees with the OLD DFT to ~1%, as P6 reported
    for l in ALLOWED:
        assert WANNIER_MEV[l] == pytest.approx(old[l] * 1e3, rel=0.011)
    # but it is 12-26% away from the corrected reference
    dev = {l: abs(REFERENCE_MEV[l] - WANNIER_MEV[l]) / WANNIER_MEV[l] for l in ALLOWED}
    assert np.median(list(dev.values())) == pytest.approx(0.124, abs=0.02)
    assert max(dev.values()) > 0.20


def test_the_noise_floor_did_not_pay_for_any_of_this():
    """The floor is unchanged, which is exactly why it never warned us."""
    _, f_old, _ = _fit(OLD)
    _, f_12, _ = _fit(NEW12)
    _, f_24, _ = _fit(NEW24)
    assert f_old == pytest.approx(27.49e-6, rel=0.02)
    assert f_12 == pytest.approx(24.30e-6, rel=0.02)
    assert f_24 == pytest.approx(23.61e-6, rel=0.02)
    # a 26% swing in the signal moved the floor by less than 20%
    assert max(f_old, f_12, f_24) / min(f_old, f_12, f_24) < 1.2
    # and the weakest allowed harmonic still sits far above it
    a24, _, _ = _fit(NEW24)
    assert a24[45] / f_24 > 35
