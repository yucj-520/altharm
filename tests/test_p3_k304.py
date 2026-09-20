"""The aliasing prediction, confirmed: doubling the k-count cuts the floor 27x.

P3 narrowed the noise floor's origin to finite-k aliasing by excluding
everything else (basis size, cutoff, diagonalisation threshold, density
convergence) and by showing the floor rises when the circle is subsampled.
This is the prediction that follows, tested directly: if aliasing is the
mechanism, a denser circle must lower the floor while leaving the physics
alone.

Independence rule: the prediction was made and recorded before the run
existed (PLAN.md P3, commit 42c8817), and the allowed harmonics act as the
control -- a mechanism that moved them too would not be aliasing.
"""

from pathlib import Path

import numpy as np
import pytest

from altharm.qe.reader import read_bands
from altharm.basis import CircularBasis
from altharm.project import project
from altharm.noise import noise_report
from altharm.symmetry import SpinOp, allowed_l

DATA = Path(__file__).parent / "data"
XML_152 = DATA / "qe_cr2seas_dense" / "data-file-schema.xml"
XML_304 = DATA / "qe_cr2seas_k304" / "data-file-schema.xml"
XML_NOD3 = DATA / "qe_cr2seas_nod3" / "data-file-schema.xml"
REPORT = DATA / "harmwave_cr2seas" / "harmwave_report.txt"

pytestmark = pytest.mark.skipif(
    not all(p.exists() for p in (XML_152, XML_304, XML_NOD3, REPORT)),
    reason="the k-density fixtures are not present",
)
LMAX = 45


def _ops():
    ops = []
    for line in REPORT.read_text().splitlines():
        if line.startswith("OP"):
            p = line.split()
            ops.append(SpinOp(rot_frac=None, trans=None, chi=int(p[1]),
                               rot_cart=np.array([float(x) for x in p[2:11]]).reshape(3, 3)))
    return ops


def _report(xml):
    qb = read_bands(xml)
    band = int(np.argmin(np.abs(qb.energies_ev[0].mean(axis=0) - qb.fermi_ev)))
    dE = qb.energies_ev[0, :, band] - qb.energies_ev[1, :, band]
    th = np.arctan2(qb.kpoints_tpiba[:, 1], qb.kpoints_tpiba[:, 0])
    dirs = np.stack([np.cos(th), np.sin(th), 0.093 * np.ones_like(th)], axis=1)
    basis = CircularBasis(LMAX)
    return noise_report(project(basis, dirs, dE), allowed_l(basis, _ops(), dirs)), len(th)


def test_doubling_the_k_count_cuts_the_floor_by_more_than_an_order_of_magnitude():
    """Measured: 27.5 -> 1.0 ueV, a factor of 27, for twice the k-points.

    Far steeper than the naive expectation. What aliases into the fitted
    range comes from the true signal's own high-l tail, so pushing the
    Nyquist limit from 76 to 152 replaces content at l~30-76 with content
    from l>152, where the spectrum is very much smaller.
    """
    r152, n152 = _report(XML_152)
    r304, n304 = _report(XML_304)
    assert (n152, n304) == (152, 304)

    assert r152.floor == pytest.approx(27.5e-6, rel=0.05)
    assert r304.floor == pytest.approx(1.01e-6, rel=0.10)
    assert r152.floor / r304.floor > 15, (
        f"only {r152.floor/r304.floor:.1f}x improvement, expected ~27x"
    )
    assert r304.floor_max == pytest.approx(2.81e-6, rel=0.10)


def test_the_physics_does_not_move_when_the_floor_does():
    """The control. A mechanism that changed the allowed harmonics too would
    not be aliasing -- it would mean the signal itself was undersampled."""
    r152, _ = _report(XML_152)
    r304, _ = _report(XML_304)
    for l in (3, 9, 15, 21, 27, 33, 39, 45):
        shift = abs(r152.allowed[l] - r304.allowed[l])
        assert shift < 2e-6, f"l={l} moved {shift/1e-6:.2f} ueV when the sampling doubled"


def test_the_weakest_harmonic_gains_a_factor_of_thirty_in_signal_to_noise():
    """What the floor is for. At 152 points the weakest allowed harmonic
    (l=39, 0.99 meV) sits ~36x above the floor; at 304 points, ~1000x. That
    is the headroom needed to go after harmonics beyond l=45."""
    r152, _ = _report(XML_152)
    r304, _ = _report(XML_304)
    assert r152.snr()[39] > 30
    assert r304.snr()[39] > 500


def test_the_dispersion_correction_does_not_touch_the_eigenvalues():
    """Grimme's correction depends only on nuclear positions, so at fixed
    geometry it cannot enter the Kohn-Sham Hamiltonian. Measured: removing
    `vdw_corr = 'grimme-d3'` leaves every eigenvalue bit-identical.

    This is what allows the defumat cross-check (P7) to run without D3 --
    which defumat does not implement -- instead of rebuilding every result
    under D2 to match it.
    """
    a = read_bands(XML_152)
    b = read_bands(XML_NOD3)
    assert np.array_equal(a.energies_ev, b.energies_ev), (
        f"max difference {np.abs(a.energies_ev - b.energies_ev).max():.3e} eV"
    )
    assert a.fermi_ev == b.fermi_ev
