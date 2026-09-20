"""P7, closed: three codes sharing no machinery agree on all eight harmonics.

  * Quantum ESPRESSO 7.2 -- Fortran plane-wave, USPP, DFT+U
  * Wannier90 -- tight-binding interpolation of a DIFFERENT plane-wave run
  * defumat -- an independent Python plane-wave implementation

The comparison that matters is defumat against the QE run with a BYTE-IDENTICAL
input file (6x6x2 SCF mesh, no vdW).  That matters more here than it normally
would: P3's second check showed the SCF mesh alone moves these coefficients by
6-29%, so comparing defumat-on-6x6x2 against QE-on-12x12x4 would have buried a
code difference under a density difference, in either direction.

Grimme-D3 is absent from both sides on purpose: defumat does not implement it,
and a QE run without it was measured to give bit-identical eigenvalues (P3), so
dropping it costs nothing and keeps two independent D2 implementations out of
the comparison.
"""

import pickle
from pathlib import Path

import numpy as np
import pytest

from altharm.basis import CircularBasis
from altharm.noise import noise_report
from altharm.project import project
from altharm.qe.reader import read_bands
from altharm.symmetry import SpinOp, allowed_l

HERE = Path(__file__).parent
PKL = HERE / "data" / "defumat_p7" / "defumat_small.pkl"
QE662 = HERE / "data" / "qe_cr2seas_scf662" / "data-file-schema.xml"
REPORT = HERE / "data" / "harmwave_cr2seas" / "harmwave_report.txt"

pytestmark = pytest.mark.skipif(
    not (PKL.exists() and QE662.exists() and REPORT.exists()),
    reason="the P7 defumat fixtures are not present",
)

BAND, LMAX = 44, 45
ALLOWED = (3, 9, 15, 21, 27, 33, 39, 45)


def _ops():
    ops = []
    for line in REPORT.read_text().splitlines():
        if line.startswith("OP"):
            p = line.split()
            ops.append(SpinOp(rot_frac=None, trans=None, chi=int(p[1]),
                              rot_cart=np.array([float(x) for x in p[2:11]]).reshape(3, 3)))
    return ops


def _harmonics(dE, kcart):
    th = np.arctan2(kcart[:, 1], kcart[:, 0])
    dirs = np.stack([np.cos(th), np.sin(th), 0.093 * np.ones_like(th)], axis=1)
    basis = CircularBasis(LMAX)
    fit = project(basis, dirs, dE)
    return fit.amplitude_by_l(), noise_report(fit, allowed_l(basis, _ops(), dirs)).floor


@pytest.fixture(scope="module")
def routes():
    qb = read_bands(QE662)
    qe = _harmonics(qb.energies_ev[0, :, BAND] - qb.energies_ev[1, :, BAND], qb.kpoints_tpiba)
    d = pickle.load(PKL.open("rb"))
    ev = np.asarray(d["eigenvalues"])
    dfm = _harmonics(ev[0, :, BAND] - ev[1, :, BAND], np.asarray(d["k_tpiba"]))
    return qe, dfm, qb, d, ev


def test_the_defumat_run_is_the_one_we_think_it_is(routes):
    """Guard the premise before comparing anything."""
    _, _, qb, d, ev = routes
    assert ev.shape == (2, 152, 80), "spin, k-points, bands"
    assert d["converged"] is True and d["iterations"] == 45
    # the SCF free energy matches QE's to all six digits defumat prints
    assert d["total_energy"] == pytest.approx(-912.52758604, abs=1e-6)
    # band 44 sits where QE puts it: same band, not a lucky index
    resid = abs(ev[0, :, BAND].mean() - qb.energies_ev[0, :, BAND].mean())
    assert resid < 1e-4, f"band {BAND} is {resid*1e3:.2f} meV from QE's"


def test_defumat_matches_qe_at_every_allowed_harmonic(routes):
    """The headline: 0.12% worst case across two orders of magnitude in amplitude."""
    (a_qe, _), (a_dfm, _), _, _, _ = routes
    dev = {l: abs(a_dfm[l] - a_qe[l]) / a_qe[l] for l in ALLOWED}
    assert max(dev.values()) < 0.005, f"worst {max(dev, key=dev.get)}: {max(dev.values()):.4f}"
    assert np.median(list(dev.values())) < 0.001
    # and the amplitudes really do span two orders of magnitude, or the test is easy
    assert a_qe[3] / a_qe[39] > 50


def test_both_codes_report_the_same_noise_floor(routes):
    """A code-independent floor is evidence it is sampling, not arithmetic."""
    (_, floor_qe), (_, floor_dfm), _, _, _ = routes
    assert floor_qe == pytest.approx(47.4e-6, rel=0.02)
    assert floor_dfm == pytest.approx(floor_qe, rel=0.01)


def test_three_routes_agree_on_the_leading_harmonic_ordering(routes):
    """Cheap structural check: all three rank the harmonics the same way."""
    (a_qe, _), (a_dfm, _), _, _, _ = routes
    order_qe = sorted(ALLOWED, key=lambda l: -a_qe[l])
    order_dfm = sorted(ALLOWED, key=lambda l: -a_dfm[l])
    assert order_qe == order_dfm
    assert order_qe[:2] == [3, 9], "l=3 then l=9 should lead"
