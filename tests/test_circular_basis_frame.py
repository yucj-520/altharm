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

from pathlib import Path

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


# --------------------------------------------------------------------------
# why the cone guard is not duplicated upstream as an orthogonality check
# --------------------------------------------------------------------------
#
# The guard's message used to blame the frame (e1, e2) only.  That is one
# cause; the other is a cell that is not the lattice its symmetry operations
# belong to, which makes their Cartesian form slightly non-orthogonal and
# tilts a rotated circle off its own cone.  A mistyped lattice constant does
# it, and the frame is then innocent.
#
# Catching that at its source -- raising when |R R^T - I| exceeds a threshold
# in _cartesian_rotation -- was implemented and REVERTED, because no threshold
# separates the two populations.  Measured over all 28 QE fixtures plus the
# two hard-coded P8 structures:
#
#   cubic / tetragonal (CoF2, Fe)                      0
#   hexagonal, exact Wyckoff (MnTe)                    2.7e-11
#   hexagonal, relaxed, from QE's XML (Cr2SeAs)        2.2e-07   <- legitimate
#   one lattice constant wrong by 16% (MnTe c)         5.2e-06   <- broken
#
# A factor of 24 between the worst legitimate case and the broken one is too
# narrow for a hard error: the project's own main material already sits at
# 2.2e-07, so a slightly rougher cell would be rejected for no reason.  The
# consequence is caught downstream instead, where it is unambiguous.
#
# The spread has one cause.  A cubic or tetragonal cell carries no sqrt(3), so
# its Cartesian rotations are exact; a hexagonal cell's C3 is only as orthogonal
# as its stored b_y approximates (sqrt(3)/2)a.  MnTe's cell is good to 1.8e-11
# there and Cr2SeAs's only to 2.9e-08 -- that is the whole four orders.  So
# "try another material" cannot rescue a threshold: the spread between
# materials IS the obstacle, and this project's own is at the bad end.
#
# Orthogonalising rot_cart (polar decomposition) was also considered and
# rejected.  It would remove the symptom for a legitimate cell and equally for
# a wrong one, handing back clean operations of the wrong structure with
# nothing raised anywhere -- the silent-wrong-answer mode this package exists
# to prevent.

def test_cr2seas_cartesian_rotations_are_not_exactly_orthogonal():
    """The measurement that killed the upstream check -- keep it visible."""
    import spglib
    from altharm.qe.reader import read_structure

    st = read_structure(Path(__file__).parent / "data" / "qe_cr2seas_sm24"
                        / "data-file-schema.xml")
    ds = spglib.get_symmetry_dataset(
        (st.cell_ang, st.scaled_positions, st.numbers), symprec=1e-4)
    rots = ds.rotations if hasattr(ds, "rotations") else ds["rotations"]
    a = np.asarray(st.cell_ang, float)
    worst = max(
        float(np.abs((a.T @ np.asarray(r, float) @ np.linalg.inv(a.T))
                     @ (a.T @ np.asarray(r, float) @ np.linalg.inv(a.T)).T
                     - np.eye(3)).max())
        for r in rots)
    assert 1e-8 < worst < 1e-6, worst       # legitimate, yet above 1e-8


def test_spin_operations_accepts_a_cell_that_is_merely_imprecise():
    """It must not raise on Cr2SeAs, whatever a tidier material would give."""
    import spglib
    from altharm.qe.reader import read_structure
    from altharm.symmetry import spin_operations

    st = read_structure(Path(__file__).parent / "data" / "qe_cr2seas_sm24"
                        / "data-file-schema.xml")
    ops = spin_operations(st.cell_ang, st.scaled_positions, st.numbers,
                          st.magmoms, symprec=1e-3)
    assert len(ops) > 0


def test_the_guard_message_names_both_causes():
    b = CircularBasis(8, e1=(1, 0, 0), e2=(0, 1, 0))
    dirs = np.array([[0.0, np.cos(t), np.sin(t)] for t in
                     np.linspace(0, 2 * np.pi, 32, endpoint=False)])
    with pytest.raises(ValueError) as e:
        b.angles(dirs)
    msg = str(e.value)
    assert "e1, e2" in msg and "cell" in msg and "non-orthogonal" in msg
