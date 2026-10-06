"""P15: the k_z = 0 residual is a rigid offset, not an angular error.

MnTe's k_z = 0 plane is nodal -- symmetry forbids every l -- so whatever
splitting is measured there is pure residual.  P8 recorded it as 19 ueV at
U = 0 falling to 0.20 ueV at U = 4, and noted only that it is "not a fixed
floor".  Decomposing it says what it is.

**Measured** (band 20, 152-point circle, allowed set empty at both U):

    U    total rms     l=0 rigid     angular      gap (20->21)
    0    19.252 ueV    19.245 ueV    0.500 ueV      759 meV
    4     0.197 ueV     0.197 ueV    0.009 ueV     1813 meV

99.96% of it is a rigid, angle-independent shift between the spin channels.
A constant in theta is l = 0, and l = 0 is not among the coefficients anyone
reports -- so the alarming-looking number was never a harmonic error.  The part
that could corrupt a harmonic is 0.5 ueV at U = 0 and 9 neV at U = 4.

By P14's channel split it sits overwhelmingly in chi = -1:

    U    chi = +1        chi = -1        ratio
    0    1.5e-11 eV      1.0e-06 eV      65000
    4    5.8e-12 eV      1.7e-08 eV       3000

which is why it moves with U rather than sitting at a fixed floor: it
measures how well the self-consistent density respects sublattice exchange, and
opening the gap 2.4x improves that 98x.

The practical consequence: a nodal plane is a good SCF diagnostic, but its
residual has to be quoted with l = 0 removed, or it overstates the angular
error by 40x (U=0) to 22x (U=4).
"""

from pathlib import Path

import numpy as np
import pytest

from altharm import CircularBasis
from altharm.qe.reader import read_bands
from altharm.symmetry import allowed_l, spin_operations

HERE = Path(__file__).parent
RUNS = {0.0: HERE / "data" / "mnte_kz0", 4.0: HERE / "data" / "mnte_kz0_u4"}
VBM = 20

CELL = np.array([[4.2200465202, 0.0, 0.0],
                 [-2.1100232601, 3.6546674915, 0.0],
                 [0.0, 0.0, 5.9460678141]])
POS = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.5],
                [0.666666723, 0.333333334, 0.75], [0.333333334, 0.666666668, 0.25]])
NUM = np.array([25, 25, 52, 52])
MAG = np.array([+1.0, -1.0, 0.0, 0.0])

pytestmark = pytest.mark.skipif(
    not all((d / "data-file-schema.xml").exists() for d in RUNS.values()),
    reason="the MnTe k_z = 0 fixtures are not present")


def _load(u):
    q = read_bands(RUNS[u] / "data-file-schema.xml")
    k = q.kpoints_tpiba
    th = np.arctan2(k[:, 1], k[:, 0])
    dirs = np.stack([np.cos(th), np.sin(th), np.zeros_like(th)], axis=1)
    dE = q.energies_ev[0, :, VBM] - q.energies_ev[1, :, VBM]
    return q, dirs, dE


def _ops():
    return spin_operations(CELL, POS, NUM, MAG, symprec=1e-4)


def test_the_plane_is_nodal_at_both_u():
    """Every l forbidden -- so the whole signal there is residual."""
    for u in RUNS:
        _, dirs, _ = _load(u)
        assert allowed_l(CircularBasis(45), _ops(), dirs) == set()


def test_the_residual_is_almost_entirely_a_rigid_offset():
    for u, floor in ((0.0, 0.99), (4.0, 0.99)):
        _, _, dE = _load(u)
        total = float(np.sqrt((dE ** 2).mean()))
        rigid = abs(float(dE.mean()))
        assert rigid / total > floor, (u, rigid / total)


def test_the_angular_part_is_tiny_and_shrinks_with_u():
    ang = {}
    for u in RUNS:
        _, _, dE = _load(u)
        ang[u] = float(np.sqrt(((dE - dE.mean()) ** 2).mean()))
    assert ang[0.0] < 1e-6          # 0.5 ueV
    assert ang[4.0] < 2e-8          # 9 neV
    assert ang[0.0] / ang[4.0] > 10


def test_the_residual_tracks_the_gap_not_a_fixed_floor():
    out = {}
    for u in RUNS:
        q, _, dE = _load(u)
        e = q.energies_ev[0]
        out[u] = (float(np.sqrt((dE ** 2).mean())),
                  float(np.min(e[:, VBM + 1] - e[:, VBM])))
    (r0, g0), (r4, g4) = out[0.0], out[4.0]
    assert g4 > g0                   # U opens the gap
    assert r4 < r0                   # and the residual falls with it
    assert r0 / r4 > 50              # by far more than the gap ratio


def test_the_residual_lives_entirely_in_the_chi_minus_one_channel():
    """Same split as P14.2 on Cr2SeAs, and the reason it is U-dependent."""
    ops = _ops()
    for u in RUNS:
        _, dirs, dE = _load(u)
        by_chi = {}
        for o in ops:
            rot = dirs @ np.linalg.inv(o.rot_cart).T
            dist = np.linalg.norm(rot[:, None, :] - dirs[None, :, :], axis=2)
            j = dist.argmin(axis=1)
            if dist[np.arange(len(dirs)), j].max() > 1e-6:
                continue
            r = dE[j] - o.chi * dE
            by_chi.setdefault(o.chi, []).append(
                float(np.sqrt(((r - r.mean()) ** 2).mean())))
        # not machine zero -- 1.5e-11 eV at U=0, 5.8e-12 at U=4 -- but three
        # to five orders below the chi = -1 channel, which is the point
        assert max(by_chi[+1]) < 1e-10, (u, max(by_chi[+1]))
        assert max(by_chi[-1]) / max(by_chi[+1]) > 1000, (u, by_chi)
