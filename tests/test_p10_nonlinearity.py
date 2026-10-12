"""P10: the high harmonics are made by band repulsion, not by long-range hopping.

This addresses the question that motivated the DFT route: whether a Wannier
model omits terms that are present in the underlying DFT calculation.  The
concern is well founded *for a linear response*: expand a plane wave on a circle
of radius k0,

    e^{i k.R} = sum_l i^l J_l(k0 |R_parallel|) e^{i l (theta - phi_R)},

so any quantity LINEAR in H(R) has its l-th harmonic weighted by J_l(k0 |R|),
which is exponentially small until k0|R| ~ l.  For Cr2Se2As, k0 = 0.4868 1/Ang,
so l = 45 would need hoppings out to ~92 Ang -- 25 unit cells.  No Wannier model
reaches that; the Cr2Se2As model used here reaches 25.8 Ang.

But dE(k) is an EIGENVALUE, not a matrix element, and eigenvalues are not linear
in H.  Level repulsion at avoided crossings is a sharp feature in theta, and a
sharp feature has high harmonics no matter how short-ranged H(R) is.

The test below is the control that separates the two channels, on one and the
same truncated Hamiltonian: two manifestly linear quantities (a Wannier-gauge
diagonal element, and the gauge-invariant trace) against the eigenvalue
splitting.  If the linear ones did NOT collapse, the Bessel argument would be
wrong; if the eigenvalue one DID collapse, the Wannier route would indeed be
blind to high harmonics and only DFT could see them.  Neither happens.
"""

import os
from pathlib import Path

import numpy as np
import pytest

from altharm.basis import CircularBasis
from altharm.kpath import circle
from altharm.project import project
from altharm import wannier

ROOT = "/scratch/work/yuc11/HIGH_THROUGHT_WANNIERcr2X3/1.cr2se2as/Cr2SeAs_wannier_v03"
HR_UP = f"{ROOT}/hamiltonian/c.up_hr.dat"
HR_DN = f"{ROOT}/hamiltonian/c.dn_hr.dat"

pytestmark = pytest.mark.skipif(
    not os.path.exists(HR_UP),
    reason="the Cr2Se2As Wannier model is not on this machine",
)

CELL = np.array([
    [3.7258396285, 0.0000002055, 0.0000000080],
    [-1.8629196361, 3.2266718606, -0.0000000204],
    [0.0000000476, -0.0000000907, 16.9179784926],
])
BAND = 22          # the Wannier-model band matching DFT band 44
LMAX = 45
ALLOWED = (3, 9, 15, 21, 27, 33, 39, 45)


@pytest.fixture(scope="module")
def setup():
    up = wannier.read_hr(HR_UP, cell=CELL)
    dn = wannier.read_hr(HR_DN, cell=CELL)
    B = 2 * np.pi * np.linalg.inv(CELL).T
    k0 = 0.25 * np.linalg.norm(B[0])
    pts = circle(k0, 152, center=(0, 0, 0.25 * np.linalg.norm(B[2])))
    frac = wannier.cart_to_frac_k(pts, CELL)
    th = np.arctan2(pts[:, 1], pts[:, 0])
    dirs = np.stack([np.cos(th), np.sin(th), 0.093 * np.ones_like(th)], axis=1)
    return up, dn, frac, dirs, CircularBasis(LMAX), k0


def _amps(basis, dirs, values):
    return project(basis, dirs, values).amplitude_by_l()


def test_circle_radius_puts_high_l_far_outside_the_bessel_reach(setup):
    """The premise: the model is far too short-ranged for a linear l=45."""
    up, _, _, _, _, k0 = setup
    rmax_inplane = np.linalg.norm(up.rcart()[:, :2], axis=1).max()
    assert 25.0 < rmax_inplane < 26.5, "the model's in-plane reach moved"
    # Bessel turns on around k0|R| ~ l; l = 45 would need ~3.5x this range.
    assert k0 * rmax_inplane < 13.0
    assert 45 / k0 > 3.0 * rmax_inplane


def test_linear_quantities_obey_the_bessel_cutoff(setup):
    """Truncated to 8 Ang, anything linear in H(R) is dead above l ~ 15."""
    up, dn, frac, dirs, basis, _ = setup
    u = wannier.truncate(up, 8.0, plane_normal=(0, 0, 1))
    d = wannier.truncate(dn, 8.0, plane_normal=(0, 0, 1))
    Hu, Hd = wannier.hk(u, frac), wannier.hk(d, frac)

    elem = _amps(basis, dirs, np.real(Hu[:, BAND, BAND] - Hd[:, BAND, BAND]))
    trace = _amps(basis, dirs, np.real(np.trace(Hu, axis1=1, axis2=2)
                                       - np.trace(Hd, axis1=1, axis2=2)))
    # l = 3 survives (k0 * 8 Ang ~ 3.9, comfortably above l = 3) ...
    assert elem[3] > 1e-4
    assert trace[3] > 1e-3
    # ... and by l = 21 both have fallen to the numerical floor.
    for l in (21, 27, 33, 39, 45):
        assert elem[l] < 1e-11, f"linear element at l={l} is not Bessel-dead"
        assert trace[l] < 1e-10, f"linear trace at l={l} is not Bessel-dead"


def test_eigenvalue_splitting_ignores_the_bessel_cutoff(setup):
    """The same 8 Ang model still carries l = 45 at ~0.9 meV: 1e10 x the linear channel."""
    up, dn, frac, dirs, basis, _ = setup
    u = wannier.truncate(up, 8.0, plane_normal=(0, 0, 1))
    d = wannier.truncate(dn, 8.0, plane_normal=(0, 0, 1))
    dE = np.linalg.eigvalsh(wannier.hk(u, frac))[:, BAND] \
        - np.linalg.eigvalsh(wannier.hk(d, frac))[:, BAND]
    amp = _amps(basis, dirs, dE)
    assert amp[45] == pytest.approx(0.8690e-3, rel=0.02)
    assert amp[39] == pytest.approx(0.7035e-3, rel=0.02)
    # every allowed harmonic is >= 8 orders of magnitude above the linear channel
    for l in (21, 27, 33, 39, 45):
        assert amp[l] > 1e-4


def test_truncation_converges_smoothly_not_exponentially(setup):
    """l = 45 grows 0.15 -> 1.14 meV over the range: a slow approach, not a switch.

    This is the quantitative form of the answer.  Truncation DOES matter -- about
    15 Ang is needed for the weakest harmonics to be within a few percent -- but
    it costs a factor of ~7, not the 40 orders of magnitude the linear picture
    predicts.
    """
    up, dn, frac, dirs, basis, _ = setup
    seen = []
    for rmax in (6.0, 8.0, 11.5, 15.0):
        u = wannier.truncate(up, rmax, plane_normal=(0, 0, 1))
        d = wannier.truncate(dn, rmax, plane_normal=(0, 0, 1))
        dE = np.linalg.eigvalsh(wannier.hk(u, frac))[:, BAND] \
            - np.linalg.eigvalsh(wannier.hk(d, frac))[:, BAND]
        seen.append(_amps(basis, dirs, dE)[45])
    full = _amps(basis, dirs, wannier.spin_splitting(up, dn, frac, BAND))[45]

    assert seen == sorted(seen), "l=45 should grow monotonically with range"
    assert seen[0] == pytest.approx(0.1546e-3, rel=0.02)   # R < 6 Ang, 29 hoppings
    assert full == pytest.approx(1.1426e-3, rel=0.02)
    assert seen[0] / full > 0.1, "even 29 hoppings give >10% of the converged l=45"
    assert seen[-1] / full > 0.95, "15 Ang should be within 5%"
