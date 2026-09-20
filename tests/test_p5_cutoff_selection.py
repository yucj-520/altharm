"""P5, the free half: the basis-set discontinuity as a function of ecutwfc.

`altharm.planewaves.count_planewaves` reproduces QE's own plane-wave count
exactly (tests/test_planewaves.py), and needs no DFT run -- so the size of
the basis-set discontinuity at any candidate cutoff can be measured before
spending a single core-hour. This file records what that scan found on the
real Cr2SeAs circle, and it is actionable rather than merely descriptive.

Independence rule: the reference is arithmetic on the lattice -- how many
reciprocal-lattice points fall inside a sphere of a given radius -- computed
by a routine already validated bit-for-bit against QE's own count. Nothing
about the DFT eigenvalues enters.
"""

from pathlib import Path

import numpy as np
import pytest

from altharm.qe.reader import read_bands
from altharm.planewaves import count_planewaves, discontinuity_spectrum
from altharm.basis import CircularBasis
from altharm.units import BOHR_ANG

XML = Path(__file__).parent / "data" / "qe_cr2seas_dense" / "data-file-schema.xml"

pytestmark = pytest.mark.skipif(
    not XML.exists(), reason="the dense-circle Cr2SeAs fixture is not present"
)


def _geometry():
    qb = read_bands(XML)
    k = qb.kpoints_tpiba * (2 * np.pi / (qb.alat_bohr * BOHR_ANG))
    th = np.arctan2(k[:, 1], k[:, 0])
    dirs = np.stack([np.cos(th), np.sin(th), np.zeros_like(th)], axis=1)
    return qb, k, dirs


def _l6_ppm(k, cell, dirs, ecut):
    """l=6 amplitude of the plane-wave count, per million plane waves.

    The relative figure is the one that matters: an absolute wobble of ten
    plane waves means something different in a basis of 8,000 than in one of
    30,000.
    """
    npw = count_planewaves(k, cell, ecutwfc_ry=float(ecut))
    amp = discontinuity_spectrum(CircularBasis(12), dirs, npw).amplitude_by_l()
    return amp[6] / npw.mean() * 1e6


def test_planewave_count_follows_the_sphere_volume_law():
    """N_pw should scale as ecut^1.5 -- the volume of the cutoff sphere.

    A power-law fit over 40-120 Ry recovers the exponent to within 0.02 of
    the analytic 1.5, which is the cheapest end-to-end check that the counter
    is counting a sphere and not something else.
    """
    qb, k, _ = _geometry()
    cuts = np.array([40.0, 50.0, 60.0, 70.0, 80.0, 100.0, 120.0])
    n = np.array([count_planewaves(k, qb.cell_ang, ecutwfc_ry=c).mean() for c in cuts])
    slope = np.polyfit(np.log(cuts), np.log(n), 1)[0]
    assert abs(slope - 1.5) < 0.02, f"N_pw ~ ecut^{slope:.3f}, expected 1.5"


def test_the_discontinuity_does_not_shrink_smoothly_with_cutoff():
    """Measured: it is an arithmetic accident of each cutoff, not a trend.

    Over 50-100 Ry in 2 Ry steps the l=6 relative amplitude swings between
    14.5 and 1167.5 ppm with no monotone behaviour -- neighbouring cutoffs
    differ by more than an order of magnitude. Which reciprocal-lattice
    points happen to sit in a thin shell at the cutoff radius is a
    number-theoretic property of that particular radius, so "raise ecutwfc
    until it converges" is the wrong instinct for this particular noise
    source.
    """
    qb, k, dirs = _geometry()
    cuts = np.arange(50, 101, 2.0)
    ppm = np.array([_l6_ppm(k, qb.cell_ang, dirs, c) for c in cuts])

    assert ppm.max() / ppm.min() > 50, (
        f"spread only {ppm.max()/ppm.min():.1f}x -- the scan looks smooth, "
        "which would contradict the recorded finding"
    )
    # not monotone: the largest value is not at the smallest cutoff
    assert cuts[ppm.argmax()] > cuts[0], "worst cutoff should not be the lowest one"
    # and a higher cutoff is not reliably cleaner
    assert ppm[cuts == 80.0][0] > ppm[cuts == 66.0][0], (
        "80 Ry should be dirtier than 66 Ry despite using ~33% more plane waves"
    )


def test_the_production_cutoff_is_the_worst_in_its_range():
    """The finding that motivated the P5 DFT runs.

    60 Ry -- the cutoff the group's production setup uses -- measures 1167 ppm,
    the largest value anywhere in 50-100 Ry, while 66 Ry measures 14.5 ppm for
    only ~15% more plane waves. If the plane-wave discontinuity is what sets
    the even-l noise floor, that is an 80x improvement available for almost
    nothing; if the 66 Ry DFT run does not show it, the mechanism is not what
    sets that floor. Either way the prediction is falsifiable, which is the
    point of recording it here.
    """
    qb, k, dirs = _geometry()
    at60, at66 = _l6_ppm(k, qb.cell_ang, dirs, 60.0), _l6_ppm(k, qb.cell_ang, dirs, 66.0)
    assert at60 == pytest.approx(1167.5, rel=0.02)
    assert at66 == pytest.approx(14.5, rel=0.05)
    assert at60 / at66 > 50

    n60 = count_planewaves(k, qb.cell_ang, ecutwfc_ry=60.0).mean()
    n66 = count_planewaves(k, qb.cell_ang, ecutwfc_ry=66.0).mean()
    assert n66 / n60 < 1.2, "the cleaner cutoff should cost well under 20% more plane waves"
