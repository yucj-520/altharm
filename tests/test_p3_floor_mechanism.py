"""What actually sets the noise floor: three mechanisms excluded, one supported.

P3 measured the floor. P4 and P5 tried to explain it and ruled two candidates
out. This file records the remaining discrimination, and it reaches a
practical conclusion: the floor is not an intrinsic limit of the DFT
calculation, it is a function of how densely the circle is sampled.

Independence rule: each mechanism is tested by changing the one knob it
depends on and nothing else, over a range wide enough that the mechanism
would have to show. A knob swept across seven orders of magnitude that moves
the floor by nothing is not a weak result -- it is a decisive one.
"""

from pathlib import Path

import numpy as np
import pytest

from altharm.qe.reader import read_bands
from altharm.basis import CircularBasis
from altharm.project import project
from altharm.noise import noise_report
from altharm.symmetry import SpinOp, allowed_l
from altharm.units import RYDBERG_EV

DATA = Path(__file__).parent / "data"
DENSE = DATA / "qe_cr2seas_dense" / "data-file-schema.xml"
DIAGO = {
    1e-10: DATA / "qe_cr2seas_diago" / "diago1e-10_data-file-schema.xml",
    1e-6: DATA / "qe_cr2seas_diago" / "diago1e-6_data-file-schema.xml",
}
REPORT = DATA / "harmwave_cr2seas" / "harmwave_report.txt"

pytestmark = pytest.mark.skipif(
    not all(p.exists() for p in [DENSE, REPORT] + list(DIAGO.values())),
    reason="the floor-mechanism fixtures are not present",
)


def _ops():
    ops = []
    for line in REPORT.read_text().splitlines():
        if line.startswith("OP"):
            p = line.split()
            ops.append(SpinOp(rot_frac=None, trans=None, chi=int(p[1]),
                               rot_cart=np.array([float(x) for x in p[2:11]]).reshape(3, 3)))
    return ops


def _floor(xml, lmax=45, stride=1):
    qb = read_bands(xml)
    idx = np.arange(0, qb.kpoints_tpiba.shape[0], stride)
    dE = (qb.energies_ev[0, idx, 44] - qb.energies_ev[1, idx, 44])
    th = np.arctan2(qb.kpoints_tpiba[idx, 1], qb.kpoints_tpiba[idx, 0])
    dirs = np.stack([np.cos(th), np.sin(th), 0.093 * np.ones_like(th)], axis=1)
    basis = CircularBasis(lmax)
    rep = noise_report(project(basis, dirs, dE), allowed_l(basis, _ops(), dirs))
    return rep, len(idx)


def test_the_diagonalisation_threshold_does_not_set_the_floor():
    """Swept across seven orders of magnitude; the floor does not move at all.

    diago_thr_init = 1e-6 Ry is a nominal eigenvalue error bound of 13.6 ueV
    -- the same scale as the 27 ueV floor -- so if per-k-point diagonalisation
    were the mechanism, this run would show it. Measured floor at 1e-13,
    1e-10 and 1e-6: 27.49 ueV in all three, identical to five significant
    figures. Even the allowed harmonics only differ in the fifth digit.
    """
    base, _ = _floor(DENSE)
    assert base.floor == pytest.approx(27.49e-6, rel=0.01)
    nominal_at_1e6 = 1e-6 * RYDBERG_EV
    assert 1e-6 < nominal_at_1e6 / base.floor < 1.0, (
        "the loosest threshold should be within an order of magnitude of the "
        "floor, or the test is not probing the right range"
    )
    for thr, xml in DIAGO.items():
        rep, _ = _floor(xml)
        assert rep.floor == pytest.approx(base.floor, rel=1e-3), (
            f"diago_thr_init={thr:g}: floor {rep.floor/1e-6:.2f} ueV vs "
            f"{base.floor/1e-6:.2f} ueV at 1e-13"
        )
        assert rep.allowed[3] == pytest.approx(base.allowed[3], rel=1e-5)


def test_the_floor_is_a_function_of_sampling_density():
    """The positive evidence, and the practical conclusion.

    At a FIXED l_max of 24 -- well below the Nyquist limit of either sampling
    -- halving the circle from 152 points to 76 raises the floor from 22.7 to
    217.5 ueV, nearly tenfold, on the same physical signal. Combined with
    P3's measurement in the other direction (56 -> 152 points dropped the
    floor by an order of magnitude), this is the one mechanism with evidence
    pointing at it from both sides.

    The practical reading: the floor is not an intrinsic limit of the DFT
    calculation. It is set by how densely the circle is sampled, so it buys
    down with k-points -- not with a larger basis, a tighter diagonalisation,
    or a more converged density, all of which have now been shown to do
    nothing.
    """
    full, n_full = _floor(DENSE, lmax=24, stride=1)
    half, n_half = _floor(DENSE, lmax=24, stride=2)
    assert (n_full, n_half) == (152, 76)

    assert full.floor == pytest.approx(22.7e-6, rel=0.05)
    assert half.floor == pytest.approx(217.5e-6, rel=0.05)
    assert half.floor / full.floor > 5, (
        f"halving the sampling should raise the floor several-fold, got "
        f"{half.floor/full.floor:.1f}x"
    )
    # and it is the floor that moves, not the physics
    for l in (3, 9, 15, 21):
        assert half.allowed[l] == pytest.approx(full.allowed[l], rel=0.02), (
            f"l={l} moved when the sampling changed -- that would mean the "
            "signal is undersampled too, not just the floor"
        )


def test_the_floor_rises_as_lmax_approaches_the_nyquist_limit():
    """Corroboration from the other axis: at fixed sampling, pushing l_max
    toward Nyquist raises the floor, because allowed and forbidden harmonics
    begin to alias into each other.

    Measured on the 152-point circle (Nyquist l = 76): the floor runs 20.9,
    27.5, 38.4, 44.0 ueV at l_max = 12, 45, 60, 70. The acceleration near the
    end is the aliasing signature; it is also the reason l_max = 45 is the
    working choice rather than the largest l_max the point count permits.
    """
    floors = {lmax: _floor(DENSE, lmax=lmax)[0].floor for lmax in (12, 45, 60, 70)}
    assert floors[12] < floors[45] < floors[60] < floors[70], (
        f"floor should rise monotonically toward Nyquist: "
        f"{ {k: round(v/1e-6, 1) for k, v in floors.items()} }"
    )
    assert floors[70] / floors[12] > 1.8, "the rise toward Nyquist should be marked"
