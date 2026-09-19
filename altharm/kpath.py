"""Circles and spheres of k-points, built in Cartesian axes.

Every path here is constructed directly in Cartesian reciprocal space and only
converted to a code's input units at the last moment.  Building a circle by
stepping uniformly in fractional coordinates instead would apply the lattice
metric to it: in a hexagonal cell the result is an ellipse, and an ellipse
sampled as if it were a circle mixes the angular harmonics into each other.
"""

from __future__ import annotations

import numpy as np


def circle(radius, n, e1=(1.0, 0.0, 0.0), e2=(0.0, 1.0, 0.0), center=(0.0, 0.0, 0.0),
           endpoint=False):
    """``n`` points on a circle of given radius, in Cartesian reciprocal units.

    Parameters
    ----------
    radius : float
        Radius, in whatever Cartesian reciprocal unit the caller is working in
        (1/Angstrom throughout this package until a writer converts it).
    n : int
        Number of points.  For a fit up to ``lmax`` this must exceed
        ``2*lmax+1``; a comfortable factor of four is the usual choice, because
        the extra points are what let the residual mean anything.
    endpoint : bool
        Leave False for a fit.  A repeated endpoint double-weights one angle and
        biases the least-squares solution.

    Returns
    -------
    (n,3) array
        The k-points, and ``angles`` via :func:`angles_of` if wanted.
    """
    e1 = np.asarray(e1, dtype=float)
    e2 = np.asarray(e2, dtype=float)
    e1 = e1 / np.linalg.norm(e1)
    e2 = e2 / np.linalg.norm(e2)
    if abs(float(e1 @ e2)) > 1e-10:
        raise ValueError("e1 and e2 must be orthogonal")
    t = np.linspace(0.0, 2.0 * np.pi, int(n), endpoint=endpoint)
    pts = radius * (np.cos(t)[:, None] * e1 + np.sin(t)[:, None] * e2)
    return pts + np.asarray(center, dtype=float)


def sphere_grid(ntheta, nphi, radius=1.0):
    """Gauss-Legendre x uniform-phi product grid on a sphere, with weights.

    Gauss-Legendre in :math:`\\cos\\theta` integrates a polynomial of degree
    ``2*ntheta-1`` exactly, and the uniform azimuth is exact for
    :math:`|m| < nphi`.  So for a band-limited function of index ``lmax`` the
    quadrature is exact once ``ntheta > lmax`` and ``nphi > 2*lmax``: the
    projection carries no quadrature error to hide behind.

    Returns
    -------
    points : (ntheta*nphi, 3)
    weights : (ntheta*nphi,)
        Summing to :math:`4\\pi`.
    """
    x, wx = np.polynomial.legendre.leggauss(int(ntheta))
    theta = np.arccos(x)
    phi = 2.0 * np.pi * np.arange(int(nphi)) / int(nphi)
    TH, PH = np.meshgrid(theta, phi, indexing="ij")
    W = np.outer(wx, np.full(int(nphi), 2.0 * np.pi / int(nphi)))
    st = np.sin(TH)
    pts = np.stack([st * np.cos(PH), st * np.sin(PH), np.cos(TH)], axis=-1)
    return radius * pts.reshape(-1, 3), W.ravel()


def angles_of(points, e1=(1.0, 0.0, 0.0), e2=(0.0, 1.0, 0.0), center=(0.0, 0.0, 0.0)):
    """Angle of each point in the plane spanned by ``e1``, ``e2``."""
    d = np.atleast_2d(np.asarray(points, dtype=float)) - np.asarray(center, dtype=float)
    e1 = np.asarray(e1, dtype=float)
    e2 = np.asarray(e2, dtype=float)
    return np.arctan2(d @ e2, d @ e1)
