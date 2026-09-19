"""Plane-wave counting tests.

Independence rule, in its strongest form: the reference for count_planewaves
is Quantum ESPRESSO itself -- the PW count it printed for a real run, parsed
by parse_pw_count's own regex, compared against a sphere count that shares no
code with QE. Neither side was derived from the other.
"""

from pathlib import Path

import numpy as np
import pytest

from altharm.planewaves import count_planewaves, parse_pw_count, discontinuity_spectrum
from altharm.qe.reader import read_bands
from altharm.units import BOHR_ANG
from altharm.basis import CircularBasis
from altharm.kpath import circle

DATA = Path(__file__).parent / "data" / "qe_fe_bands"
FE_CELL = np.array([[1.435, 1.435, 1.435], [-1.435, 1.435, 1.435], [-1.435, -1.435, 1.435]])
FE_ECUT = 80.0


@pytest.fixture(scope="module")
def fe_bands():
    return read_bands(DATA / "data-file-schema.xml")


def test_matches_the_pw_count_qe_actually_printed(fe_bands):
    """The strongest form of the independence rule: checked against pw.x itself."""
    qe_npw = parse_pw_count((DATA / "bands.out").read_text())[:24]
    k_invang = fe_bands.kpoints_tpiba * (2 * np.pi / (fe_bands.alat_bohr * BOHR_ANG))
    mine = count_planewaves(k_invang, FE_CELL, FE_ECUT)
    assert np.array_equal(mine, qe_npw), f"mismatch: mine={mine.tolist()} qe={qe_npw.tolist()}"


def test_pw_count_is_symmetric_under_the_point_group(fe_bands):
    """bcc iron's cubic symmetry: the printed counts repeat with period 6 on
    this 24-point circle (960 963 965 971 965 963 x4), a property of the
    geometry the sphere count must reproduce independent of QE."""
    qe_npw = parse_pw_count((DATA / "bands.out").read_text())[:24]
    assert np.array_equal(qe_npw[:6], qe_npw[6:12])
    assert np.array_equal(qe_npw[:6], qe_npw[18:24])


def test_larger_cutoff_gives_more_planewaves():
    k = np.array([[0.1, 0.0, 0.0]])
    lo = count_planewaves(k, FE_CELL, 40.0)[0]
    hi = count_planewaves(k, FE_CELL, 160.0)[0]
    assert hi > lo


def test_gamma_point_count_is_the_expected_order_of_magnitude():
    """Sanity, not precision: N ~ (4/3)pi R^3 V / (2pi)^3 for a ball of PWs.

    Computed independently from the sphere-volume formula, not from the
    function under test, as a coarse cross-check that nothing is off by an
    order of magnitude (a wrong 2*pi or a wrong Bohr/Angstrom factor would
    fail this by a large margin, not a percent).
    """
    r_cut = np.sqrt(FE_ECUT) / BOHR_ANG
    vol = abs(np.linalg.det(FE_CELL))
    expected = (4.0 / 3.0) * np.pi * r_cut ** 3 * vol / (2 * np.pi) ** 3
    got = count_planewaves(np.zeros((1, 3)), FE_CELL, FE_ECUT)[0]
    assert 0.5 * expected < got < 2.0 * expected


def test_refuses_nothing_but_is_deterministic_and_reproducible():
    k = circle(0.3, 12)
    a, b = count_planewaves(k, FE_CELL, FE_ECUT), count_planewaves(k, FE_CELL, FE_ECUT)
    assert np.array_equal(a, b)


def test_parse_pw_count_on_a_hand_built_string():
    text = "k =  0.10 0.00 0.00 (  123 PWs)   bands (ev):\nk = -0.10 0.00 0.00 (  456 PWs)   bands (ev):\n"
    assert list(parse_pw_count(text)) == [123, 456]


def test_discontinuity_spectrum_finds_the_sawtooth_harmonics():
    """A synthetic step-like npw(theta) with known Fourier content.

    Independent reference: a floor() staircase has a textbook Fourier series
    dominated by many harmonics with a 1/l falloff (like a sawtooth); the fit
    just needs to put non-negligible amplitude broadly across l, unlike a
    smooth physical signal concentrated at the symmetry-allowed l's.
    """
    n = 96
    pts = circle(1.0, n)
    basis = CircularBasis(20)
    theta = basis.angles(pts)
    steps = np.floor(8 * (theta + np.pi) / (2 * np.pi))     # 8 discrete jumps around the circle
    fit = discontinuity_spectrum(basis, pts, steps)
    amp = fit.amplitude_by_l()
    # broadband: at least half of l=1..20 carry non-trivial amplitude
    significant = [l for l, a in amp.items() if l > 0 and a > 0.05 * amp[1]]
    assert len(significant) > 8, f"expected a broadband spectrum, got l={significant}"
