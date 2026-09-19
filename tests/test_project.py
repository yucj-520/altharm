"""Projection tests.

Independence rule: the reference is a set of coefficients chosen by hand before
the fit ran, or an invariance that must hold on physical grounds regardless of
how the fit is implemented.
"""

import numpy as np
import pytest

from altharm.basis import CircularBasis, SphericalBasis
from altharm.kpath import circle, sphere_grid
from altharm.project import project, evaluate


def test_recovers_coefficients_chosen_by_hand():
    """Build dE(t) from literal numbers, fit it, get the same numbers back."""
    a2, b2, b6 = -0.4, 1.3, 0.75
    n = 64
    pts = circle(1.0, n)
    t = CircularBasis(1).angles(pts)
    values = a2 * np.cos(2 * t) + b2 * np.sin(2 * t) + b6 * np.sin(6 * t)

    fit = project(CircularBasis(8), pts, values)
    coeff = dict(zip(fit.labels, fit.coeffs))
    assert np.isclose(coeff[(2, "c")], a2, atol=1e-12)
    assert np.isclose(coeff[(2, "s")], b2, atol=1e-12)
    assert np.isclose(coeff[(6, "s")], b6, atol=1e-12)
    for lab, c in coeff.items():
        if lab not in {(2, "c"), (2, "s"), (6, "s")}:
            assert abs(c) < 1e-12, f"spurious weight {c:.2e} at {lab}"
    assert fit.residual_rms < 1e-12


def test_amplitude_is_independent_of_where_the_angle_starts():
    """Rotating e1 must not change the amplitude of a harmonic.

    a_l alone depends on the choice of axes and is therefore not a physical
    number; sqrt(a_l^2 + b_l^2) is.  This checks the reported quantity is the
    invariant one.
    """
    n, alpha = 96, 0.7
    pts = circle(1.0, n)
    t = CircularBasis(1).angles(pts)
    values = 0.9 * np.cos(6 * t + 0.3) + 0.2 * np.sin(2 * t)

    e1 = np.array([np.cos(alpha), np.sin(alpha), 0.0])
    e2 = np.array([-np.sin(alpha), np.cos(alpha), 0.0])

    amp_a = project(CircularBasis(8), pts, values).amplitude_by_l()
    amp_b = project(CircularBasis(8, e1=e1, e2=e2), pts, values).amplitude_by_l()
    for l in amp_a:
        assert np.isclose(amp_a[l], amp_b[l], atol=1e-10), (
            f"l={l}: {amp_a[l]:.6g} vs {amp_b[l]:.6g} after rotating the axes"
        )


def test_least_squares_agrees_with_quadrature_on_the_sphere():
    """Two routes that share no code: lstsq on samples vs a quadrature sum."""
    lmax = 4
    basis = SphericalBasis(lmax)
    pts, w = sphere_grid(lmax + 4, 2 * lmax + 4)
    rng = np.random.default_rng(3)
    true = rng.normal(size=basis.size)
    values = basis.evaluate(pts) @ true

    quad = basis.evaluate(pts).T @ (w * values)      # <Y|f>, exact for this grid
    fit = project(basis, pts, values)
    assert np.allclose(fit.coeffs, quad, atol=1e-10)
    assert np.allclose(fit.coeffs, true, atol=1e-10)


def test_round_trip_through_evaluate():
    pts = circle(1.0, 50)
    t = CircularBasis(1).angles(pts)
    values = np.cos(4 * t) - 0.5 * np.sin(t)
    basis = CircularBasis(6)
    fit = project(basis, pts, values)
    assert np.allclose(evaluate(basis, fit, pts), values, atol=1e-12)


def test_refuses_an_underdetermined_fit():
    """Fewer k-points than coefficients must raise, not quietly return garbage."""
    pts = circle(1.0, 5)
    with pytest.raises(ValueError, match="cannot determine"):
        project(CircularBasis(8), pts, np.zeros(5))


def test_refuses_mismatched_lengths():
    pts = circle(1.0, 20)
    with pytest.raises(ValueError, match="directions but"):
        project(CircularBasis(4), pts, np.zeros(19))


def test_condition_number_flags_bad_sampling():
    """Points bunched into a quarter of the circle cannot resolve lmax=8."""
    t = np.linspace(0.0, np.pi / 2, 40)
    pts = np.stack([np.cos(t), np.sin(t), np.zeros_like(t)], axis=1)
    fit = project(CircularBasis(8), pts, np.cos(6 * t))
    assert fit.cond > 1e3, f"cond = {fit.cond:.3g}, expected an obviously bad number"
