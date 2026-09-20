"""harmwave.x rerun on the group's best (152-point) data.

The original cross-check (test_harmwave_cross_check.py) used the earlier
56-point circle at lmax=12. This is the same independent-implementation
comparison -- QE's own find_sym + LAPACK DGELSS vs spglib + numpy lstsq --
run on the final 152-point, lmax=45 dataset, so the two routes are validated
against each other on the data everything else in this project's results
actually uses.
"""

from pathlib import Path

import numpy as np
import pytest

from altharm.qe.reader import read_bands
from altharm.basis import CircularBasis
from altharm.project import project

DATA = Path(__file__).parent / "data"
XML = DATA / "qe_cr2seas_dense" / "data-file-schema.xml"
REPORT_152 = DATA / "harmwave_cr2seas" / "harmwave_report_152.txt"

pytestmark = pytest.mark.skipif(
    not XML.exists() or not REPORT_152.exists(),
    reason="the 152-point harmwave fixtures are not present",
)

BAND0 = 44
LMAX = 45


def _parse(report):
    coeffs, pts = {}, []
    for line in report.read_text().splitlines():
        if line.startswith("CF"):
            p = line.split()
            coeffs[(int(p[1]), p[2])] = float(p[3])
        elif line.startswith("PT"):
            p = line.split()
            pts.append((float(p[1]), float(p[2])))
    return coeffs, np.array(pts)


def test_fortran_and_python_agree_to_machine_precision_on_the_final_dataset():
    """46 coefficients (l=0..45), Fortran DGELSS vs Python lstsq, on the
    152-point circle every other P3/P6 result in this project is built on."""
    coeffs_f, pts = _parse(REPORT_152)
    theta_f, dE_f = pts[:, 0], pts[:, 1]

    qb = read_bands(XML)
    dE_py = qb.energies_ev[0, :, BAND0] - qb.energies_ev[1, :, BAND0]
    theta_py = np.arctan2(qb.kpoints_tpiba[:, 1], qb.kpoints_tpiba[:, 0])

    assert np.allclose(theta_py, theta_f, atol=1e-8)
    assert np.allclose(dE_py, dE_f, atol=1e-9)

    dirs = np.stack([np.cos(theta_py), np.sin(theta_py), np.zeros_like(theta_py)], axis=1)
    fit = project(CircularBasis(LMAX), dirs, dE_py)
    coeffs_py = dict(zip(fit.labels, fit.coeffs))

    max_dev = 0.0
    for l in range(LMAX + 1):
        for kind in ("c", "s"):
            if l == 0 and kind == "s":
                continue
            max_dev = max(max_dev, abs(coeffs_py.get((l, kind), 0.0) - coeffs_f[(l, kind)]))
    assert max_dev < 1e-12, f"max coefficient deviation {max_dev:.2e} eV"
