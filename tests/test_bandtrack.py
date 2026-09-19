"""Band-tracking tests.

Independence rule: the reference is a hand-built two-level model whose true
(diabatic) trajectories are known in closed form -- e1(k) = k, e2(k) = -k,
crossing at k = 0 -- not anything this module produced.
"""

import numpy as np
import pytest

from altharm.bandtrack import track_by_continuity


def two_level_crossing(k):
    """e1 = k, e2 = -k: an EXACT crossing at k=0, zero coupling.

    Sorted ascending, the lower column is e1 for k<0 and e2 for k>0 -- the
    sorted index swaps character exactly at the crossing. A tracker that
    just follows "column 0" would report a false discontinuity there; a
    correct diabatic tracker must not.
    """
    e1, e2 = k, -k
    return np.sort(np.stack([e1, e2], axis=1), axis=1)


def avoided_crossing(k, gap=0.3):
    """e = +-sqrt(k^2 + gap^2): adiabatic, no exact crossing, gap = 2*|k=0 value|."""
    e = np.sqrt(k ** 2 + gap ** 2)
    return np.sort(np.stack([-e, e], axis=1), axis=1)


def test_continuity_tracker_follows_the_diabatic_branch_through_a_crossing():
    """The whole point of the module: recover e1(k) = k despite the index swap.

    Checked against the closed form, not against the sorted array.
    """
    k = np.linspace(-1.0, 1.0, 41)
    energies = two_level_crossing(k)
    tracked = track_by_continuity(energies, start_energy=k[0], wrap=False)
    assert np.allclose(tracked.energy, k, atol=1e-9), (
        f"max deviation from e1(k)=k: {np.abs(tracked.energy - k).max():.3g}"
    )


def test_naive_sorted_index_would_have_jumped_here():
    """Sanity on the test fixture itself: column 0 alone is NOT band 1.

    min(k, -k) = -|k| is continuous in *value* -- it has a kink, not a jump --
    which is exactly why naively taking "the lower sorted band" looks fine at
    a glance. What breaks is *identity*: for k>0 the lower sorted column is
    e2 = -k, not e1 = k. If this failed, the crossing test above would not be
    testing anything, because there would be no identity swap to avoid.
    """
    k = np.linspace(-1.0, 1.0, 41)
    energies = two_level_crossing(k)
    naive = energies[:, 0]
    positive = k > 0
    assert np.allclose(naive[positive], -k[positive], atol=1e-9), (
        "expected the naive lower band to have switched to e2=-k for k>0"
    )
    assert not np.allclose(naive[positive], k[positive], atol=1e-9)


def test_avoided_crossing_never_swaps_and_the_gap_matches_the_formula():
    """No exact crossing here, so continuity tracking is trivial -- but the
    reported gap must match 2*sqrt(k^2+g^2) exactly, by construction."""
    k = np.linspace(-1.0, 1.0, 41)
    g = 0.3
    energies = avoided_crossing(k, gap=g)
    tracked = track_by_continuity(energies, start_energy=-1.0, wrap=False)
    expected = -np.sqrt(k ** 2 + g ** 2)
    assert np.allclose(tracked.energy, expected, atol=1e-9)
    assert np.allclose(tracked.gap, 2 * np.sqrt(k ** 2 + g ** 2), atol=1e-9)


def test_suspect_flags_a_close_pass_from_another_band():
    """A band with steady slope, passed close by another band at one point.

    Unlike a symmetric avoided crossing -- where both branches flatten out
    together at the closest approach, so the *local step* shrinks exactly
    where the gap does, and "gap < step" never fires -- here the tracked band
    keeps a constant slope throughout. The other band only dips near it at
    k=0.3, without ever crossing. That isolates what ``suspect`` is meant to
    catch: a point where the gap to a neighbour is not comfortably larger
    than how fast the tracked band is itself moving.
    """
    k = np.linspace(-1.0, 1.0, 81)
    step = k[1] - k[0]
    e_tracked = k                                            # slope 1 throughout
    closest_gap = 0.3 * step                                 # well under the local step
    e_other = e_tracked + closest_gap + 3.0 * (k - 0.3) ** 2  # always above, dips near k=0.3
    energies = np.sort(np.stack([e_tracked, e_other], axis=1), axis=1)

    tracked = track_by_continuity(energies, start_energy=e_tracked[0], wrap=False)
    assert np.allclose(tracked.energy, e_tracked, atol=1e-9)  # never actually ambiguous here

    flagged = tracked.suspect()
    near_dip = np.abs(k - 0.3) < 3 * step
    assert flagged[near_dip].any(), "the close pass near k=0.3 should be flagged"
    assert not flagged[np.abs(k - 0.3) > 0.3].any(), "far from the dip nothing should be flagged"


def test_closure_on_a_genuine_closed_loop():
    """A band with no winding must close on itself around a full circle."""
    theta = np.linspace(0.0, 2 * np.pi, 64, endpoint=False)
    e1 = 1.0 + 0.1 * np.cos(3 * theta)     # smooth, single-valued on the circle
    e2 = -1.0 + 0.1 * np.sin(2 * theta)
    energies = np.sort(np.stack([e1, e2], axis=1), axis=1)
    tracked = track_by_continuity(energies, start_energy=e1[0], wrap=True)
    assert tracked.closure_error is not None
    assert tracked.closure_error < 0.05        # small compared to the band's own variation
    assert np.allclose(tracked.energy, e1, atol=1e-9)


def test_gap_to_nearest_other_is_a_genuine_min_distance():
    energies = np.array([[0.0, 1.0, 5.0], [0.0, 1.0, 5.0]])
    tracked = track_by_continuity(energies, start_index=1, wrap=False)
    # column 1 (value 1.0): nearest other is 0.0, distance 1.0 -- not 5.0
    assert np.allclose(tracked.gap, 1.0)


def test_refuses_a_single_point_path():
    with pytest.raises(ValueError, match="at least two"):
        track_by_continuity(np.array([[0.0, 1.0]]))


def test_start_by_energy_picks_the_nearest_column():
    energies = np.array([[0.0, 1.0, 2.0], [0.1, 1.1, 2.1]])
    tracked = track_by_continuity(energies, start_energy=1.05, wrap=False)
    assert tracked.index[0] == 1
