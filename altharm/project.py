"""Least-squares expansion of a sampled function onto an angular basis."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class HarmonicFit:
    """The result of expanding sampled values on an angular basis.

    Attributes
    ----------
    labels : list
        Basis labels, ``(l, 'c'|'s')`` for a circle or ``(l, m)`` for a sphere.
    coeffs : (size,) array
        Expansion coefficients, in the units of the input values.
    residual_rms : float
        Root-mean-square of ``values - model`` over the sample points, in the
        units of the input values.  This is the part of the signal the basis
        could not represent: if it is not small, ``lmax`` is too low or the
        sampled quantity is not smooth (a band crossing will do that).
    cond : float
        Condition number of the design matrix.  A large value means the sample
        points do not resolve the basis, and the coefficients are not to be
        trusted no matter how small the residual is.
    """

    labels: list
    coeffs: np.ndarray
    residual_rms: float
    cond: float
    values_rms: float = field(default=0.0)

    def amplitude_by_l(self):
        """Amplitude carried by each ``l``: the norm of its coefficient block.

        For a circle this is :math:`\\sqrt{a_l^2 + b_l^2}`, which is the physical
        amplitude of that harmonic and does not depend on where the angle was
        measured from.  Reporting ``a_l`` alone would depend on the choice of
        ``e1`` and is therefore not a physical number.
        """
        out = {}
        for (lab, c) in zip(self.labels, self.coeffs):
            l = int(lab[0])
            out[l] = out.get(l, 0.0) + float(c) ** 2
        return {l: np.sqrt(v) for l, v in sorted(out.items())}

    def phase_by_l(self):
        """Phase of each circular harmonic, in radians, or None where undefined."""
        a, b = {}, {}
        for (lab, c) in zip(self.labels, self.coeffs):
            if len(lab) != 2 or lab[1] not in ("c", "s"):
                return {}
            (a if lab[1] == "c" else b)[int(lab[0])] = float(c)
        out = {}
        for l in sorted(a):
            if l == 0:
                continue
            out[l] = np.arctan2(b.get(l, 0.0), a[l]) / l
        return out


def project(basis, dirs, values, weights=None):
    """Expand ``values`` sampled at ``dirs`` onto ``basis`` by least squares.

    Least squares rather than a quadrature sum, because the sample points are
    whatever the DFT run produced and need not be a quadrature rule.  Where they
    do form one, the two agree -- the tests check that.
    """
    values = np.asarray(values, dtype=float).ravel()
    design = basis.evaluate(dirs)
    if design.shape[0] != values.size:
        raise ValueError(
            f"{design.shape[0]} directions but {values.size} values"
        )
    if design.shape[0] < design.shape[1]:
        raise ValueError(
            f"{design.shape[0]} sample points cannot determine "
            f"{design.shape[1]} coefficients; use more k-points or a lower lmax"
        )
    if weights is not None:
        w = np.sqrt(np.asarray(weights, dtype=float).ravel())
        design = design * w[:, None]
        values = values * w

    coeffs, *_ = np.linalg.lstsq(design, values, rcond=None)
    model = design @ coeffs
    resid = values - model
    return HarmonicFit(
        labels=list(basis.labels),
        coeffs=coeffs,
        residual_rms=float(np.sqrt(np.mean(resid ** 2))),
        cond=float(np.linalg.cond(design)),
        values_rms=float(np.sqrt(np.mean(values ** 2))),
    )


def evaluate(basis, fit, dirs):
    """Reconstruct the fitted function at arbitrary directions."""
    return basis.evaluate(dirs) @ fit.coeffs
