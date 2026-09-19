"""Angular basis functions on the unit circle and the unit sphere.

The spin splitting of an altermagnet is a function on directions in k-space.
This module supplies the two bases we expand it in:

* :class:`CircularBasis` -- cos(l t), sin(l t) on a circle lying in a chosen
  plane of k-space.  This is the 2D case: a monolayer, or a fixed-k_z cut.
  These are orthogonal on a uniform circle but **not normalised** -- the Gram
  matrix over N sample points is diag(N, N/2, N/2, ...), not the identity.
  Nothing here depends on that, because the expansion is done by least squares
  rather than by a projection integral, but anything comparing a circular
  coefficient with a spherical one must not assume otherwise.
* :class:`SphericalBasis` -- real spherical harmonics Y_lm on a sphere,
  genuinely orthonormal.  This is the 3D case.

Both take directions as Cartesian unit vectors, so nothing in this module ever
touches fractional coordinates.  That is deliberate: expressing an angle in
fractional coordinates silently applies the lattice metric to it, which for a
hexagonal cell stretches a circle into an ellipse and mixes the harmonics.
Directions come in Cartesian, always.
"""

from __future__ import annotations

import numpy as np
from scipy.special import lpmv, gammaln


def _unit(v):
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    if np.any(n == 0.0):
        raise ValueError("zero-length direction vector")
    return v / n


def real_sph_harm(l: int, m: int, theta, phi):
    """Real spherical harmonic Y_lm(theta, phi), orthonormal on the unit sphere.

    ``theta`` is the polar angle from +z, ``phi`` the azimuth from +x.  The
    Condon-Shortley phase is inherited from ``scipy.special.lpmv``; it flips the
    sign of some components but does not affect orthonormality, which is what
    the expansion relies on and what ``tests/test_basis.py`` checks.
    """
    am = abs(m)
    if am > l:
        raise ValueError(f"|m| > l ({m}, {l})")
    # log-space normalisation: (l-|m|)!/(l+|m|)! underflows for large l
    log_norm = 0.5 * (
        np.log((2 * l + 1) / (4.0 * np.pi))
        + gammaln(l - am + 1)
        - gammaln(l + am + 1)
    )
    norm = np.exp(log_norm)
    P = lpmv(am, l, np.cos(theta))
    if m == 0:
        return norm * P
    if m > 0:
        return np.sqrt(2.0) * norm * P * np.cos(m * np.asarray(phi))
    return np.sqrt(2.0) * norm * P * np.sin(am * np.asarray(phi))


class CircularBasis:
    """cos/sin harmonics on a circle spanned by two orthonormal vectors.

    Parameters
    ----------
    lmax : int
        Highest harmonic index retained.
    e1, e2 : array_like, shape (3,)
        Orthonormal vectors spanning the plane of the circle.  The angle is
        measured from ``e1`` towards ``e2``.  Defaults to the x-y plane.
    """

    def __init__(self, lmax: int, e1=(1.0, 0.0, 0.0), e2=(0.0, 1.0, 0.0)):
        if lmax < 0:
            raise ValueError("lmax must be >= 0")
        self.lmax = int(lmax)
        self.e1 = _unit(e1)
        self.e2 = _unit(e2)
        if abs(float(self.e1 @ self.e2)) > 1e-10:
            raise ValueError("e1 and e2 must be orthogonal")
        self.labels = [(0, "c")]
        for l in range(1, self.lmax + 1):
            self.labels.append((l, "c"))
            self.labels.append((l, "s"))

    @property
    def size(self) -> int:
        return len(self.labels)

    def angles(self, dirs):
        """Cartesian directions -> angle in the plane of the circle."""
        d = _unit(np.atleast_2d(dirs))
        return np.arctan2(d @ self.e2, d @ self.e1)

    def evaluate(self, dirs):
        """Design matrix, shape (n_directions, size)."""
        t = self.angles(dirs)
        out = np.empty((t.size, self.size))
        for j, (l, kind) in enumerate(self.labels):
            out[:, j] = np.cos(l * t) if kind == "c" else np.sin(l * t)
        return out

    def l_of(self):
        return np.array([l for l, _ in self.labels], dtype=int)


class SphericalBasis:
    """Real spherical harmonics up to ``lmax`` on the unit sphere."""

    def __init__(self, lmax: int):
        if lmax < 0:
            raise ValueError("lmax must be >= 0")
        self.lmax = int(lmax)
        self.labels = [(l, m) for l in range(self.lmax + 1) for m in range(-l, l + 1)]

    @property
    def size(self) -> int:
        return len(self.labels)

    def angles(self, dirs):
        d = _unit(np.atleast_2d(dirs))
        theta = np.arccos(np.clip(d[:, 2], -1.0, 1.0))
        phi = np.arctan2(d[:, 1], d[:, 0])
        return theta, phi

    def evaluate(self, dirs):
        theta, phi = self.angles(dirs)
        out = np.empty((theta.size, self.size))
        for j, (l, m) in enumerate(self.labels):
            out[:, j] = real_sph_harm(l, m, theta, phi)
        return out

    def l_of(self):
        return np.array([l for l, _ in self.labels], dtype=int)
