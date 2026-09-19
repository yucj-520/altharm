"""Following one band's identity around a path of k-points.

Eigenvalues out of a plane-wave code come back sorted by energy at each
k-point, with no information carried over about which state is "the same
band" from one k-point to the next. Where two bands stay well separated this
does not matter -- sorted order and band identity coincide. Near a crossing or
a near-degeneracy it does: the splitting has to be
:math:`E_\\uparrow^{(n)}(k) - E_\\downarrow^{(n)}(k)` for the *same* orbital
character :math:`n` at every :math:`k`, not for whichever sorted index happens
to be smallest there. Following the sorted index instead makes
:math:`\\Delta E_n(\\theta)` jump wherever two bands change order, and a fit
mistakes that jump for a high harmonic.

This module cannot resolve character across a genuine, exact crossing -- that
needs the wavefunctions' overlap, which a plain eigenvalue list does not
carry. What it does instead:

* track a band by continuity, extrapolating from the last two points visited
  and matching to the nearest eigenvalue at the next one, which recovers the
  diabatic (character-preserving) branch through a crossing as long as the
  k-step is fine enough that no other band gets closer to the prediction, and
* report, at every step, how close the tracked value came to its nearest
  neighbour, so a fit's residual near a small gap can be attributed to
  band-tracking ambiguity rather than mistaken for physics.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class TrackedBand:
    """The result of following one band around a path.

    Attributes
    ----------
    index : (n,) int array
        Which column of the sorted eigenvalue array was picked at each point.
    energy : (n,) float array
        Its value there -- ``energies[i, index[i]]``.
    gap : (n,) float array
        Distance from ``energy[i]`` to the nearest *other* eigenvalue at the
        same point. Small compared to the step-to-step change in ``energy``
        is where tracking becomes unreliable: see :meth:`suspect`.
    closure_error : float or None
        For a closed path (``wrap=True``), the mismatch between the tracked
        energy extrapolated past the last point and the value actually found
        back at the first one. Large compared to ``gap`` signals that the
        loop does not return to the band it started from -- a genuine
        topological feature (a band inversion around the loop), not a
        tracking failure, but one that breaks the assumption that
        :math:`\\Delta E_n(\\theta)` is a single-valued function on the circle.
    """

    index: np.ndarray
    energy: np.ndarray
    gap: np.ndarray
    closure_error: float | None = None

    def suspect(self, factor: float = 1.0):
        """Points where the gap is not comfortably larger than the local step.

        ``factor`` scales the comparison: the default flags a point whenever
        its gap is smaller than the change in tracked energy from the
        previous point, which is the regime where a step to the wrong band
        costs less than staying on the right one.
        """
        step = np.abs(np.diff(self.energy, prepend=self.energy[0]))
        return self.gap < factor * np.maximum(step, np.finfo(float).eps)


def _nearest_column(row, target):
    return int(np.argmin(np.abs(row - target)))


def _gap_to_nearest_other(row, col):
    row = np.asarray(row, dtype=float)
    masked = np.delete(row, col)
    if masked.size == 0:
        return np.inf
    return float(np.min(np.abs(masked - row[col])))


def track_by_continuity(energies, start_index=0, start_energy=None, wrap=True):
    """Follow one band across a sequence of k-points by energy continuity.

    Parameters
    ----------
    energies : (n, nbnd) array
        Eigenvalues at each of ``n`` k-points, sorted ascending as a
        plane-wave code returns them (the sorting itself is not relied on
        past the first point).
    start_index : int
        Which column to start from at point 0, if ``start_energy`` is not
        given.
    start_energy : float, optional
        Start instead from whichever column at point 0 is nearest this value
        -- the usual choice is a target like the Fermi level.
    wrap : bool
        Treat the path as a closed loop (the last point's neighbour is the
        first). This is the case for a circle in k-space, which is what this
        package fits harmonics on. With ``wrap=False`` the path is treated as
        open and ``closure_error`` is left ``None``.

    Notes
    -----
    The prediction for point ``i`` is linear extrapolation from the tracked
    values at ``i-2`` and ``i-1`` (or just the value at ``i-1`` for the second
    point). This is what lets the tracker follow a *smoothly varying* band
    through an exact crossing with another one: near the crossing both bands
    have similar values, but only one of them continues the trend the tracker
    has already seen. Nearest-to-previous-value alone (no extrapolation)
    would tend to latch onto whichever band is closer at that instant, which
    is the wrong side of a crossing exactly where it matters.
    """
    energies = np.asarray(energies, dtype=float)
    n, nbnd = energies.shape
    if n < 2:
        raise ValueError("need at least two k-points to track a band")

    idx = np.empty(n, dtype=int)
    val = np.empty(n, dtype=float)
    gap = np.empty(n, dtype=float)

    idx[0] = (
        _nearest_column(energies[0], start_energy)
        if start_energy is not None
        else int(start_index)
    )
    val[0] = energies[0, idx[0]]
    gap[0] = _gap_to_nearest_other(energies[0], idx[0])

    idx[1] = _nearest_column(energies[1], val[0])
    val[1] = energies[1, idx[1]]
    gap[1] = _gap_to_nearest_other(energies[1], idx[1])

    for i in range(2, n):
        predicted = 2.0 * val[i - 1] - val[i - 2]     # linear extrapolation
        idx[i] = _nearest_column(energies[i], predicted)
        val[i] = energies[i, idx[i]]
        gap[i] = _gap_to_nearest_other(energies[i], idx[i])

    closure_error = None
    if wrap:
        predicted = 2.0 * val[-1] - val[-2]
        closure_error = float(np.abs(predicted - val[0]))

    return TrackedBand(index=idx, energy=val, gap=gap, closure_error=closure_error)
