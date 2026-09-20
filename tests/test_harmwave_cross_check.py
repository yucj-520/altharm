"""Cross-check: the QE-native Fortran tool (harmwave.x) against this
package's independent Python implementation, on the same real Cr2SeAs run.

This is the strongest form of the independence rule the project follows.
Two completely separate pipelines computed these numbers:

  Fortran (qe_native/harmwave.f90): QE's own find_sym() for symmetry, QE's
  own read_file() for the eigenvalues (Rydberg, in memory, no XML parsing),
  LAPACK's DGELSS for the least-squares fit.

  Python (altharm.qe.reader + altharm.project): a hand-parsed
  data-file-schema.xml (Hartree), numpy's lstsq.

Different language, different symmetry engine (QE's own vs spglib),
different linear-algebra library, different unit-conversion path. Agreement
between them is not a tautology -- and it was not automatic: a first version
of harmwave.f90 passed compilation cleanly but classified the wrong set of
harmonics as allowed ({0,3,6,9,12} instead of {3,9,...}), because its
time-reversal closure only doubled operations that had already survived a
plane-preservation filter, silently dropping the sublattice-exchanging
mirror's time-reversed partner. This test fixes that number down so it
cannot silently regress.
"""

import re
from pathlib import Path

import numpy as np
import pytest

from altharm.qe.reader import read_bands
from altharm.basis import CircularBasis
from altharm.kpath import circle
from altharm.project import project
from altharm.symmetry import SpinOp, allowed_harmonics

XML = Path(__file__).parent / "data" / "qe_cr2seas_p3" / "data-file-schema.xml"
REPORT = Path(__file__).parent / "data" / "harmwave_cr2seas" / "harmwave_report.txt"

pytestmark = pytest.mark.skipif(
    not XML.exists() or not REPORT.exists(),
    reason="the real Cr2SeAs fixtures are not present",
)

BAND0 = 44           # 0-based Python index == Fortran's band_index=45 (1-based).
                     # This is the band nearest E_F in the 80-band DFT calculation.
                     # An earlier version of this fixture used index 22, copied
                     # unthinkingly from the Wannier route's 38-band numbering: the
                     # 80-band DFT run has 22 extra low-lying semicore/buffer bands
                     # below the projected p-d manifold (matching Cr2SeAs_wannier_v03's
                     # own "16 semicore + 6 Se/As-4s" comment exactly), so "band 22"
                     # in the two models is not the same physical band. Caught by
                     # checking the band's mean energy against E_F -- it was 4.4 eV
                     # below it, not the ~0.04 eV the Wannier route reports.
LMAX = 12


def parse_report():
    ops, pts, coeffs = [], [], {}
    for line in REPORT.read_text().splitlines():
        if line.startswith("OP"):
            parts = line.split()
            chi = int(parts[1])
            R = np.array([float(x) for x in parts[2:11]]).reshape(3, 3)
            ops.append(SpinOp(rot_frac=None, trans=None, rot_cart=R, chi=chi))
        elif line.startswith("PT"):
            theta, de = (float(x) for x in line.split()[1:3])
            pts.append((theta, de))
        elif line.startswith("CF"):
            parts = line.split()
            l, kind, val = int(parts[1]), parts[2], float(parts[3])
            coeffs[(l, kind)] = val
    return ops, np.array(pts), coeffs


def test_fortran_and_python_read_the_same_eigenvalue_difference():
    """dE(theta) from QE's own in-memory Ry->eV path vs the XML's Ha->eV path."""
    _, pts, _ = parse_report()
    theta_f, dE_f = pts[:, 0], pts[:, 1]

    qb = read_bands(XML)
    dE_py = qb.energies_ev[0, :, BAND0] - qb.energies_ev[1, :, BAND0]
    theta_py = np.arctan2(qb.kpoints_tpiba[:, 1], qb.kpoints_tpiba[:, 0])

    # Both pipelines read the same underlying k-point list off the same run
    # (Fortran via read_file(), Python via the XML) in the same order -- no
    # reordering should be needed, and needing one would itself be a red flag.
    assert np.allclose(theta_py, theta_f, atol=1e-8), (
        f"max angle mismatch = {np.abs(theta_py - theta_f).max():.2e} rad "
        "(the two pipelines are not walking the k-points in the same order)"
    )
    assert np.allclose(dE_py, dE_f, atol=1e-9), (
        f"max |Python - Fortran| = {np.abs(dE_py - dE_f).max():.2e} eV"
    )


def test_fortran_and_python_fits_agree():
    """DGELSS (Fortran, LAPACK) vs numpy lstsq (Python): same data, same answer."""
    _, pts, coeffs_f = parse_report()
    theta, dE = pts[:, 0], pts[:, 1]
    dirs = np.stack([np.cos(theta), np.sin(theta), np.zeros_like(theta)], axis=1)

    fit = project(CircularBasis(LMAX), dirs, dE)
    coeffs_py = dict(zip(fit.labels, fit.coeffs))

    for l in range(LMAX + 1):
        py_c = coeffs_py.get((l, "c"), 0.0)
        assert np.isclose(py_c, coeffs_f[(l, "c")], atol=1e-9), f"l={l} cos"
        if l > 0:
            py_s = coeffs_py.get((l, "s"), 0.0)
            assert np.isclose(py_s, coeffs_f[(l, "s")], atol=1e-9), f"l={l} sin"


def test_qe_native_symmetry_matches_the_physical_selection_rule():
    """QE's own find_sym()-derived operations, fed to altharm's projector,
    must reproduce the {3, 9, ...} selection rule -- not the over-inclusive
    {0, 3, 6, 9, 12} a buggy time-reversal closure produced before the fix."""
    ops, pts, _ = parse_report()
    theta = pts[:, 0]
    dirs = circle(1.0, len(theta))                     # any radius; only angle matters
    # re-angle them to match the report's actual theta ordering/values
    dirs = np.stack([np.cos(theta), np.sin(theta), 0.093 * np.ones_like(theta)], axis=1)

    allowed = allowed_harmonics(CircularBasis(LMAX), ops, dirs)
    assert sorted(allowed) == [3, 9], f"got {sorted(allowed)}"


def test_the_allowed_harmonics_are_where_the_real_signal_is():
    """Close the loop: the harmonics QE's own symmetry engine allows are
    exactly the ones the independently-fitted real data puts its signal in.

    Measured: l=3 at 59.1 meV, l=9 at 38.9 meV, against a largest forbidden
    harmonic of 1.22 meV (l=11) -- roughly 30-50x, not the 100x+ this package
    sees on a cleaner circle (P3's 56-point sample is not dense enough to
    separate genuine higher-l signal from the forbidden floor at this l_max;
    see PLAN.md P3's honest accounting of what changed once the band index
    was corrected).
    """
    _, pts, coeffs_f = parse_report()
    amp = {
        l: (coeffs_f.get((l, "c"), 0.0) ** 2 + coeffs_f.get((l, "s"), 0.0) ** 2) ** 0.5
        for l in range(LMAX + 1)
    }
    allowed_amp = [amp[l] for l in (3, 9)]
    forbidden_amp = [amp[l] for l in range(LMAX + 1) if l not in (3, 9)]
    assert min(allowed_amp) > 20 * max(forbidden_amp), (
        f"allowed {allowed_amp} not well above forbidden {max(forbidden_amp):.2e}"
    )
    assert np.isclose(amp[3], 59.13e-3, rtol=0.01)
    assert np.isclose(amp[9], 38.86e-3, rtol=0.01)
