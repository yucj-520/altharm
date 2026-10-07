"""Band tracking applied to the real Cr2Se2As circle -- what it found there.

`tests/test_bandtrack.py` checks the tracker against hand-built two-level
models. This file is the other half: what happens when it is pointed at the
actual 152-point DFT circle the production analysis was run on, which nothing
had ever done. It turned up a real feature of the band structure, and it
justifies (rather than assumes) the fixed-band-index convention the harmonic
fits use.

Independence rule: the symmetry claim below is checked against the raw
eigenvalues' own gap structure -- a quantity the harmonic fit never touches
and `altharm.symmetry` never sees.
"""

from pathlib import Path

import numpy as np
import pytest

from altharm.qe.reader import read_bands
from altharm.bandtrack import track_by_continuity

XML = Path(__file__).parent / "data" / "qe_cr2seas_dense" / "data-file-schema.xml"

pytestmark = pytest.mark.skipif(
    not XML.exists(), reason="the dense-circle Cr2Se2As fixture is not present"
)

BAND = 44           # the band nearest E_F; the one every harmonic fit uses
NEIGHBOUR = 43


@pytest.fixture(scope="module")
def bands():
    return read_bands(XML)


def _gap(bands, spin):
    e = bands.energies_ev[spin]
    return e[:, BAND] - e[:, NEIGHBOUR]


def test_bands_43_and_44_never_actually_cross():
    """The fixed-index convention is legitimate here, and this is why.

    Measured: the gap between bands 43 and 44 runs 18.0-276.2 meV around the
    circle and is never zero. They undergo *avoided* crossings, not real ones,
    so the sorted index 44 is a continuous, single-valued adiabatic branch all
    the way round -- which is what makes `energies_ev[:, :, 44]` a well-defined
    function of theta rather than a splice of two different bands.
    """
    qb = read_bands(XML)
    for spin in (0, 1):
        gap = _gap(qb, spin)
        assert gap.min() > 0, f"spin {spin}: bands 43/44 touch -- the index convention breaks"
        assert np.isclose(gap.min() * 1e3, 18.0, atol=0.3), (
            f"spin {spin}: min gap moved to {gap.min()*1e3:.2f} meV (was 18.0)"
        )
        assert np.isclose(gap.max() * 1e3, 276.2, atol=1.0)


def test_six_avoided_crossings_per_loop():
    """Six local gap minima around the circle, in three depth classes.

    Measured (spin up): minima at k = 13, 38, 63, 89, 114, 139 -- spacing
    152/6 -- with depths 95.7, 18.0, 95.7, 25.7, 95.5, 25.7 meV. Sharp angular
    features like these are the plausible reason this band's Delta E(theta)
    needs harmonics past l=45 to converge (PLAN.md P3).
    """
    qb = read_bands(XML)
    gap = _gap(qb, 0)
    nk = len(gap)
    minima = [i for i in range(nk)
              if gap[i] <= gap[(i - 1) % nk] and gap[i] <= gap[(i + 1) % nk]]
    assert minima == [13, 38, 63, 89, 114, 139], f"minima moved to {minima}"
    depths_meV = sorted(round(gap[i] * 1e3, 1) for i in minima)
    assert depths_meV[0] == pytest.approx(18.0, abs=0.3)
    assert depths_meV[-1] == pytest.approx(95.7, abs=0.5)


def test_the_two_spin_channels_gap_curves_differ_by_exactly_half_a_turn():
    """An independent confirmation of the chi = -1 C_2z symmetry, from raw
    eigenvalues rather than from the fitted harmonics.

    altharm.symmetry derives, from the structure alone, that M_z combined with
    time reversal gives C_2z with character -1, i.e. E_up(theta) =
    E_down(theta + pi). If that is true, any band-structure observable of the
    up channel must equal the down channel's rotated by half a turn. Measured
    on the 43-44 gap curve: rolling the down channel by nk/2 = 76 points
    reproduces the up channel to 0.004 meV. Nothing in the harmonic fit or the
    symmetry module is involved in either side of that comparison.
    """
    qb = read_bands(XML)
    gu, gd = _gap(qb, 0), _gap(qb, 1)
    nk = len(gu)
    assert nk % 2 == 0
    dev_meV = np.abs(np.roll(gd, nk // 2) - gu).max() * 1e3
    assert dev_meV < 0.05, f"half-turn symmetry broken by {dev_meV:.4f} meV"
    # and it is specifically a HALF turn, not any roll
    worse = min(np.abs(np.roll(gd, s) - gu).max() * 1e3
                for s in (nk // 4, nk // 3, nk // 6))
    assert worse > 1.0, "some other rotation matches just as well -- the test proves nothing"


def test_the_continuity_tracker_switches_index_at_the_tight_crossings():
    """The tracker follows the diabatic branch, and says so -- documented limit.

    Started on band 44, `track_by_continuity` hands back a trajectory that
    swaps to 43 and back at the 18 meV avoided crossings, because its linear
    extrapolation overshoots a branch that is turning over. That is the
    eigenvalue-only limitation `bandtrack.py`'s own docstring states, not a
    defect -- but it does mean the tracker must NOT be used to pick the band
    for these harmonic fits, and this test pins that behaviour so the choice
    stays deliberate.
    """
    qb = read_bands(XML)
    tracked = track_by_continuity(qb.energies_ev[0], start_index=BAND, wrap=True)
    assert set(np.unique(tracked.index)) == {NEIGHBOUR, BAND}
    switches = np.where(np.diff(tracked.index) != 0)[0]
    assert len(switches) >= 2, "expected the tracker to swap index at the tight crossings"
    # it flags the ambiguity rather than hiding it
    assert tracked.suspect().any(), "the tracker should flag at least one suspect point here"
    # and the loop does not close on the band it started from -- the tell-tale
    assert tracked.closure_error > 0.05, (
        "an odd number of index swaps should leave a large closure error"
    )
