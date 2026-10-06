"""P18: why the Taylor-vs-Fourier cross-check cannot be run where P7 wanted it.

P7 proposed checking the fit by comparing a Taylor expansion in k against the
Fourier coefficients from a sampled circle.  The content of that check needs no
defumat: if dE(k) were analytic at the circle's centre, the lowest-degree
polynomial producing cos(l theta) would be degree l, so a_l ~ k0**l as k0 -> 0.

**Measured** on the v03 Wannier model, radii 0.003-0.034 1/A at k_z = 0.25 b3:

    band   l   fitted exponent   Taylor says   R^2
      20   3         1.42             3       0.995
      21   3         1.33             3       0.996
      22   3         1.11             3       0.993
      21   9        -0.22             9       0.27

l = 3 scales roughly linearly, not cubically, and l = 9 is not a power law.

The reason is symmetry, not noise: k_x = k_y = 0 lies on the Gamma-A line,
whose little group has two-dimensional irreps, so bands are degenerate there --
measured gaps of 0.005, 0.009 and 0.010 meV for bands 19/20, 21/22, 24/25.
dE has a square-root singularity at the centre, and an exponent near 1 is that
singularity's fingerprint.

Which corroborates P10 from the other side: P10 found the high harmonics come
from the non-linearity of the eigenvalue problem rather than hopping range, and
that same non-analyticity is what defeats the Taylor expansion here.
"""

from pathlib import Path

import numpy as np
import pytest

from altharm.wannier import cart_to_frac_k, eigenvalues, read_hr, spin_splitting

WANN = Path("/scratch/work/yuc11/HIGH_THROUGHT_WANNIERcr2X3/1.cr2se2as/"
            "Cr2SeAs_wannier_v03/hamiltonian")
CELL = np.array([[3.7258396285, 2.055e-7, 8.0e-9],
                 [-1.8629196361, 3.2266718606, -2.04e-8],
                 [4.76e-8, -9.07e-8, 16.9179784926]])

pytestmark = pytest.mark.skipif(
    not (WANN / "c.up_hr.dat").exists(),
    reason="the v03 Wannier hr.dat lives outside the repo and is not present")


def _models():
    return (read_hr(WANN / "c.up_hr.dat", cell=CELL),
            read_hr(WANN / "c.dn_hr.dat", cell=CELL))


def _kz():
    return 0.25 * (2 * np.pi * np.linalg.inv(CELL).T)[2]


def _a_l(up, dn, k0, band, l, n=256):
    th = 2 * np.pi * np.arange(n) / n
    kc = np.stack([k0 * np.cos(th), k0 * np.sin(th), np.zeros(n)], axis=1) + _kz()
    dE = spin_splitting(up, dn, cart_to_frac_k(kc, CELL), band)
    return 2 * float(np.mean(dE * np.sin(l * th)))


def test_the_circle_centre_sits_on_a_degenerate_line():
    """The symmetry fact that makes the Taylor route inapplicable."""
    up, _ = _models()
    e = eigenvalues(up, cart_to_frac_k(_kz()[None, :], CELL))[0]
    pairs = [(19, 20), (21, 22), (24, 25)]
    for lo, hi in pairs:
        assert abs(e[hi] - e[lo]) < 2e-5, (lo, hi, e[hi] - e[lo])


def test_a_l_does_not_scale_as_k0_to_the_l():
    """If it did, the fit could be cross-checked against a Taylor expansion."""
    up, dn = _models()
    ks = np.array([0.003, 0.005, 0.008, 0.013, 0.021, 0.034])
    for band in (20, 21, 22):
        v = np.array([abs(_a_l(up, dn, k, band, 3)) for k in ks])
        slope = np.polyfit(np.log(ks), np.log(v), 1)[0]
        assert slope < 2.0, (band, slope)        # nowhere near 3


def test_the_leading_harmonic_scales_about_linearly():
    """The fingerprint of a square-root singularity at the centre."""
    up, dn = _models()
    ks = np.array([0.003, 0.005, 0.008, 0.013, 0.021, 0.034])
    for band in (20, 21, 22):
        v = np.array([abs(_a_l(up, dn, k, band, 3)) for k in ks])
        p, _ = np.polyfit(np.log(ks), np.log(v), 1)
        assert 0.8 < p < 1.8, (band, p)


def test_the_ninth_harmonic_is_not_a_power_law_at_all():
    up, dn = _models()
    ks = np.array([0.003, 0.005, 0.008, 0.013, 0.021, 0.034])
    v = np.array([abs(_a_l(up, dn, k, 21, 9)) for k in ks])
    p = np.polyfit(np.log(ks), np.log(v), 1)
    pred = np.polyval(p, np.log(ks))
    r2 = 1 - np.sum((np.log(v) - pred) ** 2) / np.sum(
        (np.log(v) - np.log(v).mean()) ** 2)
    assert r2 < 0.6, r2
