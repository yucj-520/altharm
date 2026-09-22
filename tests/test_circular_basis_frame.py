"""The circle's plane and the basis frame have to be the same plane.

`CircularBasis` takes `e1`/`e2` and measures the angle in that frame.  Handing
it directions from a *different* plane is not an error the code could catch by
type or shape: the arctan2 still returns numbers, the design matrix still has
the right shape, and the fit still converges.  What comes out is a
parametrisation that is not monotone along the sample, and an allowed set that
looks plausible and is wrong.

Measured on Cr2SeAs, a 152-point circle in the y-z plane:

| basis frame | d.n spread | off-block max in S | allowed l |
|---|---|---|---|
| default (x, y) -- wrong | 1.96 | **165** | {3,6,8,9,10,11,12,13,14,15,16,19,21} |
| (y, z) -- correct | 5.6e-17 | 1.4e-07 | **{2,4,6,...,24}** |

The correct frame gives the even harmonics, which is what the four surviving
operations require: one of them is a two-fold axis sending theta -> theta+pi.
The wrong frame gives thirteen scattered values with no group-theoretic reading,
because the symmetrizer is no longer block-diagonal in l and
`allowed_harmonics`'s per-block diagonalisation has lost its premise.

Nothing raised.  The guard added to `angles()` is what makes it raise now, and
this file is the record of what it costs to not have one.
"""

import numpy as np
import pytest

from altharm.basis import CircularBasis
from altharm.kpath import circle


def _yz_circle(radius=0.5, offset=0.1, n=200):
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.stack([offset * np.ones(n), radius * np.cos(th), radius * np.sin(th)], axis=1)


def test_the_default_frame_refuses_a_yz_circle():
    with pytest.raises(ValueError, match="not on a cone"):
        CircularBasis(24).angles(_yz_circle())


def test_the_matching_frame_accepts_it():
    b = CircularBasis(24, e1=(0, 1, 0), e2=(0, 0, 1))
    th = b.angles(_yz_circle())
    assert th.shape == (200,)
    # and the angle really is the circle's own parameter, monotone up to the wrap
    d = np.diff(th)
    assert np.sum(d < 0) == 1, "exactly one wrap at +pi"


def test_a_great_circle_through_the_origin_is_still_fine():
    """d.n = 0 is the other valid case, and the commonest one."""
    b = CircularBasis(12)
    pts = circle(0.5, 200, center=(0, 0, 0))
    assert b.angles(pts).shape == (200,)


def test_a_cone_at_fixed_kz_is_still_fine():
    """d.n = const != 0 -- what every Cr2SeAs run in this repository uses."""
    b = CircularBasis(12)
    pts = circle(0.5, 200, center=(0, 0, 0.1))
    d = pts / np.linalg.norm(pts, axis=1, keepdims=True)
    assert np.ptp(d @ np.cross(b.e1, b.e2)) < 1e-12
    assert b.angles(pts).shape == (200,)


def test_the_guard_does_not_fire_on_two_directions():
    """Two points define no cone; the check needs a sample to have an opinion."""
    b = CircularBasis(4)
    assert b.angles(np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.3]])).shape == (2,)


def test_evaluate_inherits_the_guard():
    with pytest.raises(ValueError, match="not on a cone"):
        CircularBasis(8).evaluate(_yz_circle())
