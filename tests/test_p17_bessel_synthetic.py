"""P17: the Bessel suppression, on a synthetic hr.dat where the answer is known.

P10's truncation study asked for this check in so many words: "build a
*synthetic* hr.dat containing a known harmonic at a controlled radius, truncate
it at a series of R_max values, and confirm the measured damping matches the
J_l(k0 R_max) prediction rather than some other falloff" (PLAN.md P10).

The argument is Jacobi-Anger.  A hopping at in-plane displacement R, angle phi,
contributes exp(i k0 R cos(theta - phi)) to H(k) on a circle of radius k0, whose
l-th Fourier coefficient is i^l J_l(k0 R) exp(-i l phi).  So for a LINEAR
quantity -- a matrix element, a trace -- the l-th harmonic from a shell at
radius R is exactly proportional to J_l(k0 R).

**Measured, single shell** (square lattice a = 3 A, one orbital, hoppings on one
shell with angular pattern cos(l phi), k0 = 0.9 1/A).  The ratio of the measured
harmonic to J_l(k0 R) is constant across ten shells to **4e-15**, for each of
l = 2, 4, 6, 8.  That is not a falloff that merely resembles a Bessel function;
it is the Bessel function, zeros and sign changes included.

One lattice artifact had to be excluded first.  A shell whose radius is a
Pythagorean hypotenuse is not a single angular orbit: |R| = 5 holds 12 sites in
3 distinct orbits ((5,0), (3,4), (4,3)), so cos(l phi) on it does not excite a
single l.  Those shells deviate; every other shell is exact.  Measuring without
noticing that would have read as a 5% failure of the Bessel law.

**Measured, truncation.**  For a model with many shells the harmonic is a sum
over everything retained, so it does NOT track J_l(k0 R_max) as a smooth
envelope -- the ratio wanders over a factor of 70.  What survives is the
threshold: J_l(x) is negligible for x << l, so a model truncated below
R ~ l/k0 cannot represent harmonic l at all.  At half that radius the harmonic
is suppressed by 10^3 to 10^5 relative to the converged value, and the
suppression is of the ORDER of J_l(k0 R_max), within about a decade.

This is consistent with P10's note that the falloff "is not monotonic -- it
overshoots before settling", and it is why a 4-7 A model gets l = 21 wrong by a
factor of six while a 15 A one does not.
"""

import numpy as np
import pytest

from altharm.wannier import WannierHam, cart_to_frac_k, hk, truncate

jv = pytest.importorskip("scipy.special").jv


def _cell(a, c=10.0):
    return np.array([[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, c]])


def _harmonic(wh, cell, k0, l, n=1024):
    """l-th cosine harmonic of a linear quantity -- here the single element."""
    th = 2 * np.pi * np.arange(n) / n
    kc = np.stack([k0 * np.cos(th), k0 * np.sin(th), np.zeros(n)], axis=1)
    h = hk(wh, cart_to_frac_k(kc, cell))[:, 0, 0].real
    return 2 * float(np.mean(h * np.cos(l * th)))


def _shell_sites(ns):
    return [(i, j) for i in range(-ns, ns + 1) for j in range(-ns, ns + 1)
            if abs(np.hypot(i, j) - ns) < 1e-9]


def _single_orbit(ns):
    """True when |R| = ns is one angular orbit, i.e. ns is not a hypotenuse."""
    pts = _shell_sites(ns)
    return len({round(np.arctan2(j, i) % (np.pi / 2), 6) for i, j in pts}) == 1


def _shell_model(ns, l, a):
    rv, val = [], []
    for i, j in _shell_sites(ns):
        rv.append([i, j, 0])
        val.append(np.cos(l * np.arctan2(j, i)))
    rv = np.array(rv, dtype=int)
    return WannierHam(rvecs=rv, deg=np.ones(len(rv), dtype=int),
                      ham=np.array(val, float).reshape(-1, 1, 1).astype(complex),
                      cell=_cell(a))


def _decaying_model(rmax_cells, a, scale, seed=0):
    rng = np.random.default_rng(seed)
    rv, val = [], []
    for i in range(-rmax_cells, rmax_cells + 1):
        for j in range(-rmax_cells, rmax_cells + 1):
            r = np.hypot(i, j)
            if r == 0 or r > rmax_cells:
                continue
            rv.append([i, j, 0])
            val.append(np.exp(-r / scale) * rng.normal())
    rv = np.array(rv, dtype=int)
    return WannierHam(rvecs=rv, deg=np.ones(len(rv), dtype=int),
                      ham=np.array(val, float).reshape(-1, 1, 1).astype(complex),
                      cell=_cell(a))


# --------------------------------------------------------------------------

def test_a_pythagorean_shell_is_not_one_angular_orbit():
    """The artifact that must be excluded before the law can be seen."""
    assert not _single_orbit(5)          # (5,0), (3,4), (4,3)
    assert not _single_orbit(10)
    assert len(_shell_sites(5)) == 12
    for ns in (1, 2, 3, 4, 6, 7, 8, 9):
        assert _single_orbit(ns), ns


@pytest.mark.parametrize("l", [2, 4, 6, 8])
def test_a_single_shell_gives_exactly_the_bessel_function(l):
    a, k0 = 3.0, 0.9
    ratios = []
    for ns in range(1, 13):
        if not _single_orbit(ns):
            continue
        b = jv(l, k0 * ns * a)
        if abs(b) < 1e-6:                # skip near a zero of J_l
            continue
        ratios.append(_harmonic(_shell_model(ns, l, a), _cell(a), k0, l) / b)
    ratios = np.array(ratios)
    assert len(ratios) >= 8, len(ratios)
    assert np.std(ratios) / abs(np.mean(ratios)) < 1e-9, ratios


def test_truncation_below_the_bessel_threshold_kills_the_harmonic():
    """R ~ l/k0 is a threshold, not a smooth envelope."""
    a, k0 = 1.0, 1.0
    full = _decaying_model(40, a, scale=25.0)
    for l in (8, 12, 16, 20):
        thr = l / k0
        below = _harmonic(truncate(full, thr / 2, plane_normal=(0, 0, 1)),
                          _cell(a), k0, l)
        above = _harmonic(truncate(full, 2 * thr, plane_normal=(0, 0, 1)),
                          _cell(a), k0, l)
        assert abs(above) / abs(below) > 100.0, (l, below, above)


def test_the_truncated_model_does_not_track_jl_as_an_envelope():
    """Because it sums over every retained shell -- so do not claim it does."""
    a, k0 = 3.0, 0.9
    full = _decaying_model(12, a, scale=8.0)
    ratios = []
    for rm in (6, 9, 12, 15, 18, 21, 24):
        h = abs(_harmonic(truncate(full, rm, plane_normal=(0, 0, 1)),
                          _cell(a), k0, 6))
        b = abs(jv(6, k0 * rm))
        if b > 1e-9:
            ratios.append(h / b)
    assert max(ratios) / min(ratios) > 10.0, ratios
