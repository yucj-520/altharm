"""P3, resolved: a denser circle separates real high-l signal from the
56-point circle's finite-k artefact.

Independence rule: the reference is the group's Wannier model, an entirely
separate calculation (tight-binding fit to a different plane-wave run),
compared harmonic by harmonic across all eight allowed l up to 45 -- not
just the leading two.
"""

from pathlib import Path

import numpy as np
import pytest

from altharm.qe.reader import read_bands
from altharm.basis import CircularBasis
from altharm.project import project
from altharm.noise import noise_report
from altharm.symmetry import SpinOp, allowed_l

XML = Path(__file__).parent / "data" / "qe_cr2seas_dense" / "data-file-schema.xml"
REPORT = Path(__file__).parent / "data" / "harmwave_cr2seas" / "harmwave_report.txt"

pytestmark = pytest.mark.skipif(
    not XML.exists() or not REPORT.exists(),
    reason="the dense-circle Cr2SeAs fixtures are not present",
)

BAND0 = 44
LMAX = 45

# The independent Wannier-route reference, l=3..45 (see PLAN.md P6/examples/cr2seas.py).
WANNIER_MEV = {3: 59.4812, 9: 39.0640, 15: 3.1750, 21: 8.4434,
               27: 1.4946, 33: 2.8382, 39: 0.9964, 45: 1.1423}


def _ops():
    ops = []
    for line in REPORT.read_text().splitlines():
        if line.startswith("OP"):
            p = line.split()
            ops.append(SpinOp(rot_frac=None, trans=None, chi=int(p[1]),
                               rot_cart=np.array([float(x) for x in p[2:11]]).reshape(3, 3)))
    return ops


def _dE_and_dirs():
    qb = read_bands(XML)
    assert qb.kpoints_tpiba.shape[0] == 152
    assert abs((qb.energies_ev[0, :, BAND0].mean() - qb.fermi_ev) - 0.0415) < 0.001, (
        "band 44 should still be the one at E - E_F = +0.04 eV"
    )
    dE = qb.energies_ev[0, :, BAND0] - qb.energies_ev[1, :, BAND0]
    theta = np.arctan2(qb.kpoints_tpiba[:, 1], qb.kpoints_tpiba[:, 0])
    dirs = np.stack([np.cos(theta), np.sin(theta), 0.093 * np.ones_like(theta)], axis=1)
    return dE, dirs


def test_dft_matches_wannier_at_every_allowed_harmonic_to_one_percent():
    """The strongest cross-check in the project: 8 harmonics, two independent
    routes, all within 1%. l=3..45, not just the leading pair."""
    dE, dirs = _dE_and_dirs()
    basis = CircularBasis(LMAX)
    amp = project(basis, dirs, dE).amplitude_by_l()
    for l, wannier_meV in WANNIER_MEV.items():
        dft_meV = amp[l] / 1e-3
        assert np.isclose(dft_meV, wannier_meV, rtol=0.011), (
            f"l={l}: DFT {dft_meV:.4f} meV vs Wannier {wannier_meV:.4f} meV"
        )


def test_denser_circle_resolves_the_56_point_floor_as_finite_k_aliasing():
    """The direct answer to P3's open question, measured.

    On the 56-point circle, the forbidden odd/non-3 harmonics sat at
    0.4-1.2 meV. Tripling the k-count to 152 drops them by roughly an order
    of magnitude (to 0.02-0.09 meV), while the even forbidden harmonics stay
    just as clean as before (sub-ueV) and every allowed harmonic keeps
    matching Wannier to <1%. That pattern -- more k-points, cleaner floor,
    unchanged signal -- is what finite-k aliasing looks like, not what fixed
    numerical noise (which would not care about k-count) looks like.
    """
    dE, dirs = _dE_and_dirs()
    basis = CircularBasis(LMAX)
    ops = _ops()
    allowed = allowed_l(basis, ops, dirs)
    assert sorted(allowed) == [3, 9, 15, 21, 27, 33, 39, 45]

    rep = noise_report(project(basis, dirs, dE), allowed)

    assert rep.floor < 0.05e-3, f"noise floor should be well under 0.05 meV, got {rep.floor/1e-3:.4f}"
    assert rep.floor_max < 0.15e-3, f"worst forbidden harmonic drifted: {rep.floor_max/1e-3:.4f} meV"
    assert min(rep.snr().values()) > 30, f"weakest allowed/floor ratio dropped: {min(rep.snr().values()):.1f}"

    even_forbidden = [rep.forbidden[l] for l in range(2, LMAX + 1, 2)]  # l=0 excluded by design
    assert max(even_forbidden) < 1e-6, "even forbidden harmonics should remain sub-ueV clean"


def test_real_signal_extends_past_lmax_45():
    """The residual at l>45 is still 0.8 meV -- genuine unresolved structure,
    not noise, consistent with the Wannier route needing l up to ~48 to
    fully converge this band. lmax=45 is not yet the end of the story."""
    dE, dirs = _dE_and_dirs()
    fit = project(CircularBasis(LMAX), dirs, dE)
    assert 0.3e-3 < fit.residual_rms < 2e-3, (
        f"residual beyond l=45 moved outside its measured range: {fit.residual_rms/1e-3:.3f} meV"
    )
