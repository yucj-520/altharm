"""P16: the floor follows the tail's decay, and plateaus at the chi=-1 error.

P3 left open whether the forbidden-harmonic floor falls as some power of the
sample count or saturates at a genuine noise level.  P14 makes the question
answerable without a new run: the floor IS the unresolved tail folded down, so
an independent measurement of the tail predicts it.

The 56-, 152- and 304-point Cr2SeAs circles share one density (60 Ry, degauss
0.0037), and the 304-point fit resolves l to 151 -- which is exactly the tail
the shorter circles fold.  Reading the floor only from the channels aliasing
feeds (l' forbidden, N-l' allowed), as P14.3 requires:

    N    channels   measured floor   predicted from 304 tail   ratio
    56       3       1050.4 ueV            986.2 ueV           1.065
   152       7         54.67 ueV            53.84 ueV          1.015

So the floor tracks the tail's own decay -- a property of the material's
harmonic ladder -- not a power of N.

It stops falling where the chi = -1 constraints stop being exact, since those
are what make a forbidden harmonic forbidden.  On this density chi = -1 is
0.766 ueV and chi = +1 is 0.072 ueV; the 304-point floor is 1.865 ueV, 2.4x
above.  Cr2SeAs at 304 points is approaching that crossover, not yet at it.
"""

from pathlib import Path

import numpy as np
import pytest

from altharm.qe.reader import read_bands, read_structure
from altharm.symmetry import spin_operations

DATA = Path(__file__).parent / "data"
FIX = {56: DATA / "qe_cr2seas_p3", 152: DATA / "qe_cr2seas_dense",
       304: DATA / "qe_cr2seas_k304"}
BAND = 44
ALLOWED = lambda l: l % 6 == 3          # noqa: E731

pytestmark = pytest.mark.skipif(
    not all((d / "data-file-schema.xml").exists() for d in FIX.values()),
    reason="the 56/152/304-point Cr2SeAs fixtures are not present")


def _fit(n, lmax):
    q = read_bands(FIX[n] / "data-file-schema.xml")
    k = q.kpoints_tpiba
    th = np.arctan2(k[:, 1], k[:, 0]) % (2 * np.pi)
    order = np.argsort(th)
    th = th[order]
    dE = (q.energies_ev[0, :, BAND] - q.energies_ev[1, :, BAND])[order]
    d = np.column_stack([np.ones(len(th))] +
                        [g(l * th) for l in range(1, lmax + 1)
                         for g in (np.cos, np.sin)])
    c, *_ = np.linalg.lstsq(d, dE, rcond=None)
    return {(l, kk): c[1 + 2 * (l - 1) + j]
            for l in range(1, lmax + 1) for j, kk in enumerate("cs")}


def _fed_channels(n, lmax):
    """Forbidden l' that aliasing actually feeds, paired with their source.

    No cap on the source: for N = 304 the sources are 259..303, which are what
    that circle's own floor is made of.  The 151 cap belongs only to the
    prediction test, where the source has to be something the 304-point fit
    can actually resolve.
    """
    return [(lp, n - lp) for lp in range(1, lmax + 1)
            if not ALLOWED(lp) and n - lp >= 1 and ALLOWED(n - lp)]


def _resolvable(pairs):
    """Those pairs whose source the 304-point fit reaches (l <= 151)."""
    return [(lp, src) for lp, src in pairs if src <= 151]


def test_the_three_circles_share_one_density():
    for n in (56, 304):
        txt = (FIX[n] / "bands.in").read_text()
        assert "degauss = 0.00370000" in txt, n
        assert "ecutwfc = 60.000000" in txt, n


def test_the_floor_is_predicted_by_the_independently_measured_tail():
    hi = _fit(304, 151)
    for n, lmax, tol in ((56, 20, 0.10), (152, 45, 0.05)):
        co = _fit(n, lmax)
        pairs = _resolvable(_fed_channels(n, lmax))
        assert len(pairs) >= 3, (n, len(pairs))
        got = np.median([abs(co[(lp, "s")]) for lp, _ in pairs])
        pred = np.median([abs(hi[(src, "s")]) for _, src in pairs])
        assert abs(got / pred - 1.0) < tol, (n, got, pred)


def test_the_floor_falls_by_far_more_than_any_power_of_n_would_give():
    """56 -> 152 is 2.7x in N but 19x in the floor: not a 1/N law."""
    f = {}
    for n, lmax in ((56, 20), (152, 45)):
        co = _fit(n, lmax)
        f[n] = float(np.median([abs(co[(lp, "s")])
                                for lp, _ in _fed_channels(n, lmax)]))
    assert f[56] / f[152] > 15.0
    assert (152 / 56) ** 2 < f[56] / f[152]       # steeper than 1/N^2 too


def test_the_plateau_is_the_chi_minus_one_self_consistency_error():
    q = read_bands(FIX[304] / "data-file-schema.xml")
    st = read_structure(FIX[304] / "data-file-schema.xml")
    k = q.kpoints_tpiba
    dE = q.energies_ev[0, :, BAND] - q.energies_ev[1, :, BAND]
    d = k / np.linalg.norm(k, axis=1, keepdims=True)
    ops = spin_operations(st.cell_ang, st.scaled_positions, st.numbers,
                          st.magmoms, symprec=1e-3)
    by_chi = {}
    for o in ops:
        rot = d @ np.linalg.inv(o.rot_cart).T
        dist = np.linalg.norm(rot[:, None, :] - d[None, :, :], axis=2)
        j = dist.argmin(axis=1)
        if dist[np.arange(len(d)), j].max() > 1e-6:
            continue
        r = dE[j] - o.chi * dE
        by_chi.setdefault(o.chi, []).append(
            float(np.sqrt(((r - r.mean()) ** 2).mean())))

    chi_minus = max(by_chi[-1])
    chi_plus = max(by_chi[+1])
    assert chi_minus > 10 * chi_plus                   # the two channels differ

    co = _fit(304, 45)
    floor = float(np.median([abs(co[(lp, "s")])
                             for lp, _ in _fed_channels(304, 45)]))
    # approaching the crossover, not yet at it: a factor of a few, not orders
    assert 1.0 < floor / chi_minus < 10.0, (floor, chi_minus)
