"""P20: the 3D spherical route, on a real material for the first time.

SphericalBasis and the spherical projector had unit tests and nothing else --
no DFT run had ever put a sphere of k-points through them.  Job 20814237 does:
288 points (12 x 24 Gauss-Legendre x uniform) on a sphere of radius 0.293879
tpiba, chosen to pass through the existing 152-point circle so the two routes
can be compared.  Same sm24 density, band 44.

**What it settles.**  The symmetry projector predicts allowed l = 4, 6, 8, 10
without touching the data, and the fit agrees: every forbidden l sits at
0.00000 meV except l = 2 at 0.00090, against allowed coefficients of tens of
meV -- a signal-to-forbidden ratio of 110 000.  The first allowed l is 4.

That resolves a recurring confusion.  "g-wave means l=4" and this project's own
"MnTe is g-wave with leading order l=3" are both right, in different bases:
expanding dE ~ k_z k_perp^3 cos(3 phi) in real spherical harmonics gives
EXACTLY ONE term, Y_4^3, with the other 48 coefficients below 1e-6 and a
relative residual of 6.2e-17.  Spherical l=4 is the "g"; circular l=3 is the
azimuthal index left after k_z is fixed.  The 3D route makes that explicit,
which the 2D route cannot -- and that, rather than being one more cross-check,
is what the route is for.

**What it does not settle: the amplitudes.**  The fit leaves 10.2% residual at
lmax=10 and puts l=8 above l=4, which is aliasing, not physics.  Predicting the
measured circle from the spherical fit confirms it:

    lmax   spherical residual   predicted circle vs measured
      4          22.2%                   111.6%
      6          21.6%                    98.5%
      8          11.4%                    60.9%
     10          10.2%                    43.1%

Monotone, nowhere near converged.  The grid is the cause: 12 x 24 is exact only
for l < 12, while Cr2Se2As carries a ladder out to l = 45 on the circle, and
everything above l=10 folds down into what is below.  Closing it needs about
46 x 92 = 4232 points, roughly 20 hours -- 14x this run.
"""

from pathlib import Path

import numpy as np
import pytest

from altharm import SphericalBasis, evaluate, project, sphere_grid
from altharm.qe.reader import read_bands
from altharm.symmetry import allowed_harmonics, spin_operations

DATA = Path(__file__).parent / "data"
SPHERE = DATA / "sphere" / "data-file-schema.xml"
WEIGHTS = DATA / "sphere" / "sphere_w.npy"
CIRCLE = DATA / "qe_cr2seas_sm24" / "data-file-schema.xml"
BAND = 44

pytestmark = pytest.mark.skipif(
    not (SPHERE.exists() and WEIGHTS.exists()),
    reason="the spherical fixture is not present")


def _sphere():
    q = read_bands(SPHERE)
    k = q.kpoints_tpiba
    d = k / np.linalg.norm(k, axis=1, keepdims=True)
    dE = q.energies_ev[0, :, BAND] - q.energies_ev[1, :, BAND]
    return d, dE, np.load(WEIGHTS)


def _ops():
    txt = (Path("/scratch/work/yuc11/altharm_runs/hw_yz/bands.in")).read_text().splitlines()
    i = next(n for n, l in enumerate(txt) if "CELL_PARAMETERS" in l)
    cell = np.array([[float(x) for x in txt[i + 1 + j].split()] for j in range(3)]) * 7.040816482142
    j = next(n for n, l in enumerate(txt) if "ATOMIC_POSITIONS" in l)
    sym = [l.split()[0] for l in txt[j + 1:j + 11]]
    pos = np.array([[float(x) for x in l.split()[1:4]] for l in txt[j + 1:j + 11]])
    z = {"CrU": 24, "CrD": 24, "Se": 34, "As": 33}
    mag = [1.0 if s == "CrU" else (-1.0 if s == "CrD" else 0.0) for s in sym]
    return spin_operations(cell, pos, [z[s] for s in sym], mag, symprec=1e-3)


def test_the_g_wave_form_is_exactly_one_spherical_harmonic():
    """Analytic, no DFT: dE ~ k_z k_perp^3 cos(3 phi) is Y_4^3 and nothing else.

    This is what reconciles "g-wave means l=4" with the circular l=3.
    """
    b = SphericalBasis(6)
    pts, w = sphere_grid(60, 120, radius=1.0)
    d = pts / np.linalg.norm(pts, axis=1, keepdims=True)
    x, y, z = d[:, 0], d[:, 1], d[:, 2]
    f = z * np.hypot(x, y) ** 3 * np.cos(3 * np.arctan2(y, x))
    fit = project(b, d, f, weights=w)
    big = [(l, m, c) for (l, m), c in zip(fit.labels, fit.coeffs) if abs(c) > 1e-6]
    assert len(big) == 1, big
    assert (big[0][0], big[0][1]) == (4, 3)
    assert fit.residual_rms / np.sqrt((f ** 2).mean()) < 1e-14


def test_the_sphere_passes_through_the_circle():
    """Otherwise the two routes are not comparable."""
    qs = read_bands(SPHERE)
    qc = read_bands(CIRCLE)
    rs = np.linalg.norm(qs.kpoints_tpiba[0])
    rc = np.linalg.norm(qc.kpoints_tpiba[0])
    assert abs(rs - rc) < 1e-6, (rs, rc)
    assert len(qs.kpoints_tpiba) == 288


def test_the_allowed_set_is_four_six_eight_ten():
    d, _, _ = _sphere()
    assert sorted(allowed_harmonics(SphericalBasis(10), _ops(), d)) == [4, 6, 8, 10]


def test_forbidden_channels_are_five_orders_below_the_signal():
    d, dE, w = _sphere()
    b = SphericalBasis(10)
    allowed = set(allowed_harmonics(b, _ops(), d))
    f = project(b, d, dE, weights=w)
    by_l = {}
    for (l, m), c in zip(f.labels, f.coeffs):
        by_l.setdefault(l, []).append(abs(c))
    sig = max(max(by_l[l]) for l in allowed)
    forb = max(max(by_l[l]) for l in range(1, 11) if l not in allowed)
    assert sig / forb > 1e4, (sig, forb)


def test_the_amplitudes_are_not_converged_and_must_not_be_quoted():
    """10% residual at lmax=10, and the circle prediction is 43% off."""
    d, dE, w = _sphere()
    qc = read_bands(CIRCLE)
    kc = qc.kpoints_tpiba
    dc = kc / np.linalg.norm(kc, axis=1, keepdims=True)
    dEc = qc.energies_ev[0, :, BAND] - qc.energies_ev[1, :, BAND]

    rel = {}
    for lmax in (4, 6, 8, 10):
        b = SphericalBasis(lmax)
        f = project(b, d, dE, weights=w)
        pred = evaluate(b, f, dc)
        rel[lmax] = (f.residual_rms / np.sqrt((dE ** 2).mean()),
                     np.sqrt(((pred - dEc) ** 2).mean()) / np.sqrt((dEc ** 2).mean()))

    assert rel[10][0] > 0.05, rel          # still 10% unexplained
    assert rel[10][1] > 0.3, rel           # circle prediction still 43% off
    # both improve monotonically with lmax -- it is resolution, not a bug
    assert [rel[l][1] for l in (4, 6, 8, 10)] == sorted(
        [rel[l][1] for l in (4, 6, 8, 10)], reverse=True)


def test_l8_exceeding_l4_is_the_aliasing_signature():
    """Not physics: the grid cannot reach the l > 10 content that folds down."""
    d, dE, w = _sphere()
    b = SphericalBasis(10)
    f = project(b, d, dE, weights=w)
    by_l = {}
    for (l, m), c in zip(f.labels, f.coeffs):
        by_l.setdefault(l, []).append(abs(c))
    assert max(by_l[8]) > max(by_l[4])
