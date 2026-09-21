"""P13: the pseudopotential and the functional -- the last two untested inputs.

P12 established that the Hubbard U moves the coefficients 21% per eV and that
the forbidden-harmonic floor cannot see it.  Two inputs of the same kind were
still untested, and this closes them.  Both are compared against the corrected
reference (12x12x4, degauss 0.0150, U = 3.0, band 44, 152-point circle) with
only the one knob changed.

| input | median | worst | floor sees it? |
|---|---|---|---|
| Hubbard U, per eV | 21% | -- | no |
| **PAW vs ultrasoft** | **9.7%** | 17.8% (l=45) | no |
| **PBEsol vs PBE** | **7.3%** | 33.2% (l=39) | no |
| smearing set too narrow | 12.9% | -- | no |
| geometry at 1e-3 A | <4% | -- | partly |
| *the floor's own number* | *0.04%* | -- | -- |

**They eat the high harmonics specifically.**  l = 3 moves -1.3% (PBEsol) and
+0.4% (PAW); l = 39 moves -33.2% and -12.1%, l = 45 moves -21.7% and -17.8%.
So neither is a rescaling: both leave the leading harmonic essentially intact
and attack the tail -- which is exactly the part this project exists to measure.

**And the floor goes DOWN while the coefficients go wrong.**  24.30 -> 8.15
(PBEsol) and 14.28 ueV (PAW), because the floor tracks the residual above
l_max (0.795 -> 0.495, 0.647 meV), not the error.  That is P12's mechanism
correction seen at its sharpest: a run whose coefficients are 20% off reports a
*better* noise floor than the reference.

**A constraint found on the way, which is general rather than incidental.**
The first attempt used SG15 ONCV norm-conserving pseudopotentials -- a bigger
contrast with ultrasoft than PAW is.  It cannot be done: their `<PP_PSWFC>`
block is empty, and DFT+U needs atomic wavefunctions to build the Hubbard
projectors, so `pw.x` stops with `determine_hubbard_occ (1): no atomic
wavefunctions in pseudopotential file`.  **A Hubbard calculation restricts which
pseudopotential families are available at all**, which is a real limitation on
this kind of sensitivity test rather than a misconfiguration.

A second trap on the same path: pslibrary's default Se PAW (`-dn`) carries 3d in
valence, 16 electrons rather than 6.  That takes the cell from 90 electrons to
130 and band 44 stops being the same band.  The `-n` variant restores 6, so all
three runs have 90 electrons and 54 KS states -- checked below, because the
CoF2 result in `test_p8_cof2.py` is what happens when that check is skipped.
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
RUNS = {"base": "qe_cr2seas_sm12", "pbesol": "qe_cr2seas_xc_pbesol",
        "paw": "qe_cr2seas_ps_paw"}
REPORT = HERE / "data" / "harmwave_cr2seas" / "harmwave_report.txt"

pytestmark = pytest.mark.skipif(
    not (all((HERE / "data" / d / "data-file-schema.xml").exists() for d in RUNS.values())
         and REPORT.exists()),
    reason="the P13 fixtures are not present",
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


@pytest.fixture(scope="module")
def runs():
    out = {}
    for tag, folder in RUNS.items():
        q = read_bands(HERE / "data" / folder / "data-file-schema.xml")
        th = np.arctan2(q.kpoints_tpiba[:, 1], q.kpoints_tpiba[:, 0])
        dirs = np.stack([np.cos(th), np.sin(th), 0.093 * np.ones_like(th)], axis=1)
        basis = CircularBasis(LMAX)
        f = project(basis, dirs, q.energies_ev[0, :, BAND] - q.energies_ev[1, :, BAND])
        floor = noise_report(f, allowed_l(basis, _ops(), dirs)).floor
        out[tag] = (f.amplitude_by_l(), floor, f, q)
    return out


def test_all_three_runs_describe_the_same_band():
    """The premise, checked rather than assumed -- see P8's CoF2 failure."""
    for tag, folder in RUNS.items():
        text = (HERE / "data" / folder / "scf.in").read_text()
        assert "degauss     = 0.0150" in text and "12 12 4" in text
        assert "U CrU-3d 3.0" in text
    for tag in RUNS:
        q = read_bands(HERE / "data" / RUNS[tag] / "data-file-schema.xml")
        assert q.energies_ev.shape[1] == 152
        # band 44 within 5 meV of the same distance above E_F in all three
        gap = q.energies_ev[0, :, BAND].mean() - q.fermi_ev
        assert 0.030 < gap < 0.040, f"{tag}: band {BAND} at {gap:+.4f} eV"


def test_only_one_knob_differs_in_each_comparison():
    base = (HERE / "data" / RUNS["base"] / "scf.in").read_text()
    pbesol = (HERE / "data" / RUNS["pbesol"] / "scf.in").read_text()
    paw = (HERE / "data" / RUNS["paw"] / "scf.in").read_text()
    # PBEsol: functional only, same pseudopotentials
    assert "PBESOL" in pbesol and "uspp" in pbesol
    assert base.replace("'PBE'", "'PBESOL'") == pbesol
    # PAW: pseudopotentials (and the ecutrho they require), same functional
    assert "kjpaw" in paw and "'PBE'" in paw and "PBESOL" not in paw
    assert "ecutrho     = 480" in paw and "ecutwfc     = 60" in paw


def test_neither_input_is_a_rescaling(runs):
    """Both leave l=3 alone and eat the tail -- the opposite of harmless."""
    base, pbesol, paw = (runs[t][0] for t in ("base", "pbesol", "paw"))
    for other in (pbesol, paw):
        assert abs(other[3] - base[3]) / base[3] < 0.02, "l=3 barely moves"
    assert abs(pbesol[39] - base[39]) / base[39] == pytest.approx(0.332, abs=0.03)
    assert abs(paw[45] - base[45]) / base[45] == pytest.approx(0.178, abs=0.03)
    # the ratio of the worst shift to the l=3 shift is large for both
    for other in (pbesol, paw):
        worst = max(abs(other[l] - base[l]) / base[l] for l in ALLOWED)
        assert worst / (abs(other[3] - base[3]) / base[3]) > 10


def test_the_sizes_against_the_other_known_blind_spots(runs):
    base, pbesol, paw = (runs[t][0] for t in ("base", "pbesol", "paw"))
    med = lambda o: float(np.median([abs(o[l] - base[l]) / base[l] for l in ALLOWED]))
    assert med(pbesol) == pytest.approx(0.073, abs=0.02)
    assert med(paw) == pytest.approx(0.097, abs=0.02)
    # both are far above the floor's own claim (0.04%) and comparable to smearing
    for o in (pbesol, paw):
        assert med(o) > 100 * 4e-4


def test_the_floor_improves_while_the_answer_gets_worse(runs):
    """The sharpest single demonstration that the floor is not an error bar."""
    f_base, f_pbesol, f_paw = (runs[t][1] for t in ("base", "pbesol", "paw"))
    assert f_base == pytest.approx(24.30e-6, rel=0.03)
    assert f_pbesol == pytest.approx(8.15e-6, rel=0.05)
    assert f_paw == pytest.approx(14.28e-6, rel=0.05)
    # both report a BETTER floor than the reference while being 7-10% off
    assert f_pbesol < f_base and f_paw < f_base
    # and it tracks the residual, as P12 established
    r_base, r_pbesol, r_paw = (runs[t][2].residual_rms for t in ("base", "pbesol", "paw"))
    assert r_pbesol < r_base and r_paw < r_base
    assert (f_pbesol / f_base) == pytest.approx(r_pbesol / r_base, rel=0.5)
