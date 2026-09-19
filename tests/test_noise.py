"""Noise-floor tests.

Independence rule: the noise is injected at a known amplitude before the fit,
and the report must recover that amplitude.  The reference is the injected
number, not anything the fit produced.
"""

import numpy as np

from altharm.basis import CircularBasis
from altharm.kpath import circle
from altharm.project import project
from altharm.noise import noise_report
from altharm.symmetry import spin_operations, allowed_l
from tests.test_symmetry import rutile


def test_floor_recovers_an_injected_forbidden_amplitude():
    """Put a known amount into a forbidden harmonic; the floor must report it."""
    n = 120
    pts = circle(1.0, n)
    t = CircularBasis(1).angles(pts)
    signal, leak = 1.0e-3, 4.0e-5          # 1 meV signal, 40 ueV of fake noise
    values = signal * np.sin(6 * t) + leak * np.cos(5 * t)

    fit = project(CircularBasis(10), pts, values)
    rep = noise_report(fit, allowed_l_set={2, 6, 10})

    assert np.isclose(rep.allowed[6], signal, atol=1e-12)
    assert np.isclose(rep.forbidden[5], leak, atol=1e-12)
    # one forbidden harmonic out of l = 1,3,4,5,7,8,9 carries the leak
    n_forbidden = len(rep.forbidden)
    assert np.isclose(rep.floor, leak / np.sqrt(n_forbidden), atol=1e-12)
    assert rep.floor_max == rep.forbidden[5]


def test_signal_to_noise_ranks_what_is_believable():
    n = 120
    pts = circle(1.0, n)
    t = CircularBasis(1).angles(pts)
    rng = np.random.default_rng(11)
    values = (
        1.0e-3 * np.sin(2 * t)            # strong, believable
        + 2.0e-6 * np.sin(6 * t)          # weak, below the noise we inject
        + 3.0e-5 * rng.normal(size=n)     # white noise, spread over all l
    )
    fit = project(CircularBasis(10), pts, values)
    rep = noise_report(fit, allowed_l_set={2, 6, 10})

    assert rep.snr()[2] > 20.0, "the strong harmonic should be far above the floor"
    assert rep.snr()[6] < 5.0, "the buried harmonic must not look credible"
    assert 2 in rep.credible(threshold=5.0)
    assert 6 not in rep.credible(threshold=5.0)


def test_l0_is_excluded_from_the_floor_by_default():
    """A rigid offset between the spin channels is not angular noise."""
    n = 60
    pts = circle(1.0, n)
    t = CircularBasis(1).angles(pts)
    values = 5.0e-3 + 1.0e-3 * np.sin(2 * t)      # big constant offset
    fit = project(CircularBasis(6), pts, values)
    rep = noise_report(fit, allowed_l_set={2, 6})
    assert 0 not in rep.forbidden
    assert rep.floor < 1e-12, f"the offset leaked into the floor: {rep.floor:.3g}"


def test_end_to_end_on_rutile_symmetry():
    """Full chain: structure -> allowed l -> circle -> fit -> noise report.

    The splitting is built to respect rutile's symmetry, so every forbidden
    harmonic the fit returns is pure numerical noise, and here that noise is a
    known injected amount.
    """
    lmax = 10
    n = 4 * lmax + 8
    basis = CircularBasis(lmax)
    pts = circle(1.0, n)
    ops = spin_operations(*rutile())
    allowed = allowed_l(basis, ops, pts)
    assert allowed == {2, 6, 10}

    t = basis.angles(pts)
    noise_amp = 1.0e-5
    rng = np.random.default_rng(7)
    values = 1.2e-3 * np.sin(2 * t) + 8.0e-5 * np.sin(6 * t)
    values = values + noise_amp * rng.normal(size=n)

    rep = noise_report(project(basis, pts, values), allowed)

    assert np.isclose(rep.allowed[2], 1.2e-3, rtol=0.05)
    assert np.isclose(rep.allowed[6], 8.0e-5, rtol=0.2)
    # white noise of amplitude s spread over n points lands as ~s*sqrt(2/n) per harmonic
    assert rep.floor < 5.0 * noise_amp * np.sqrt(2.0 / n)
    assert rep.snr()[2] > 50.0
    assert "noise floor" in rep.summary(unit_name="meV", unit_scale=1e-3)
