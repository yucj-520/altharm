"""P3's second independent check: re-run the same circle from a DIFFERENT
converged density.

The phase asked for this from the beginning and it stayed open until now. The
answer is a caveat, not a reassurance: the SCF k-mesh propagates into the
harmonic coefficients at the 6-29% level, far above every other uncertainty
this project has measured (P5's ~20 ueV per-coefficient spread from ecutwfc,
P3's 27 ueV noise floor).

Two densities, both self-consistently converged to 1e-10 Ry, both then used for
an identical frozen-density 152-point bands run on the same circle:

  * 12x12x4 (69 irreducible k-points) -- the group's production mesh
  *  6x6x2  (14 irreducible k-points) -- deliberately coarse

These are NOT two samples of the same quantity with noise between them. The
coarse mesh gives a genuinely different self-consistent state: E_F moves by
21 meV and band 44 sits at +0.0242 eV instead of +0.0415 eV above it. The bands
have moved, so the harmonics moved with them.

What this closes: the floor measured in P3 is reproducible only at fixed SCF
mesh. What it opens: the production 12x12x4 mesh has never itself been shown
converged -- 6->12 changed the answer by up to 29%, so 12->24 might too. That is
now the largest known uncertainty in the whole pipeline and it is NOT covered by
the forbidden-harmonic floor, which measures the fit, not the density.
"""

from pathlib import Path

import numpy as np
import pytest

from altharm.basis import CircularBasis
from altharm.noise import noise_report
from altharm.project import project
from altharm.qe.reader import read_bands
from altharm.symmetry import SpinOp, allowed_l

HERE = Path(__file__).parent
DENSE = HERE / "data" / "qe_cr2seas_dense" / "data-file-schema.xml"     # 12x12x4
COARSE = HERE / "data" / "qe_cr2seas_scf662" / "data-file-schema.xml"   # 6x6x2
REPORT = HERE / "data" / "harmwave_cr2seas" / "harmwave_report.txt"

pytestmark = pytest.mark.skipif(
    not (DENSE.exists() and COARSE.exists() and REPORT.exists()),
    reason="the two-density Cr2SeAs fixtures are not present",
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


def _fit(xml):
    qb = read_bands(xml)
    assert qb.kpoints_tpiba.shape[0] == 152
    th = np.arctan2(qb.kpoints_tpiba[:, 1], qb.kpoints_tpiba[:, 0])
    dirs = np.stack([np.cos(th), np.sin(th), 0.093 * np.ones_like(th)], axis=1)
    dE = qb.energies_ev[0, :, BAND] - qb.energies_ev[1, :, BAND]
    fit = project(CircularBasis(LMAX), dirs, dE)
    floor = noise_report(fit, allowed_l(CircularBasis(LMAX), _ops(), dirs)).floor
    return fit.amplitude_by_l(), floor, qb


def test_the_two_runs_really_are_different_densities():
    """Guard the premise: if these were the same density the test proves nothing."""
    _, _, dense = _fit(DENSE)
    _, _, coarse = _fit(COARSE)
    assert dense.fermi_ev == pytest.approx(8.1169, abs=1e-3)
    assert coarse.fermi_ev == pytest.approx(8.1379, abs=1e-3)
    # 21 meV apart -- a different self-consistent state, not a repeat of one.
    assert abs(dense.fermi_ev - coarse.fermi_ev) > 0.015
    # and band 44 sits at a different height above it in each
    hi = dense.energies_ev[0, :, BAND].mean() - dense.fermi_ev
    lo = coarse.energies_ev[0, :, BAND].mean() - coarse.fermi_ev
    assert hi == pytest.approx(0.0415, abs=0.002)
    assert lo == pytest.approx(0.0242, abs=0.002)


def test_scf_kmesh_moves_the_harmonics_by_up_to_thirty_percent():
    """The number this phase owes: how much does the density matter?"""
    a_dense, _, _ = _fit(DENSE)
    a_coarse, _, _ = _fit(COARSE)
    dev = {l: abs(a_coarse[l] - a_dense[l]) / a_dense[l] for l in ALLOWED}

    # Not noise-level: every single allowed harmonic moves by more than 5%.
    assert min(dev.values()) > 0.05, "a coarser density should not be this harmless"
    assert max(dev.values()) == pytest.approx(0.289, abs=0.03)
    assert max(dev, key=dev.get) == 45, "the weakest harmonic should be the most fragile"

    # The leading harmonics move least -- the effect grows toward high l.
    assert dev[3] < 0.10
    assert dev[45] > 0.25


def test_the_density_effect_dwarfs_every_other_uncertainty():
    """Put the three measured uncertainties on one scale, at l = 3."""
    a_dense, floor_dense, _ = _fit(DENSE)
    a_coarse, floor_coarse, _ = _fit(COARSE)

    from_density = abs(a_coarse[3] - a_dense[3])          # this test
    from_cutoff = 20e-6                                   # P5, eV per coefficient
    from_fit = floor_dense                                # P3 noise floor

    assert from_fit == pytest.approx(27.49e-6, rel=0.02)
    assert from_density > 100 * from_cutoff
    assert from_density > 100 * from_fit
    # The floor itself is density-dependent too, though far less so.
    assert floor_coarse == pytest.approx(47.4e-6, rel=0.05)
