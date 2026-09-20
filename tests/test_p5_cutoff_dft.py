"""P5's expensive half: the basis-set-discontinuity hypothesis, falsified.

`tests/test_p5_cutoff_selection.py` measured, with no DFT at all, that the
plane-wave count's l=6 angular amplitude swings 80x between 60, 66 and 80 Ry.
That made a prediction: if that discontinuity is what sets the even-l noise
floor in the eigenvalues, the floor should swing with it.

It does not. Three real scf+bands runs at those three cutoffs put the
measured floor within 6% of itself while the predicted mechanism varies by a
factor of 80. The mechanism is ruled out -- for even l as well as the odd l
P4 had already excluded it for.

Independence rule: the prediction comes from lattice arithmetic (a G-vector
count), the measurement from DFT eigenvalues. Neither computes the other.
"""

from pathlib import Path

import numpy as np
import pytest

from altharm.qe.reader import read_bands
from altharm.basis import CircularBasis
from altharm.project import project
from altharm.noise import noise_report
from altharm.symmetry import SpinOp, allowed_l

DATA = Path(__file__).parent / "data"
XMLS = {
    60: DATA / "qe_cr2seas_dense" / "data-file-schema.xml",
    66: DATA / "qe_cr2seas_p5" / "ecut66_data-file-schema.xml",
    80: DATA / "qe_cr2seas_p5" / "ecut80_data-file-schema.xml",
}
REPORT = DATA / "harmwave_cr2seas" / "harmwave_report.txt"
LMAX = 45

# what count_planewaves predicted for each cutoff, in ppm of the mean count
PREDICTED_PPM = {60: 1167.5, 66: 14.5, 80: 975.2}

pytestmark = pytest.mark.skipif(
    not all(p.exists() for p in list(XMLS.values()) + [REPORT]),
    reason="the P5 cutoff fixtures are not present",
)


def _ops():
    ops = []
    for line in REPORT.read_text().splitlines():
        if line.startswith("OP"):
            p = line.split()
            ops.append(SpinOp(rot_frac=None, trans=None, chi=int(p[1]),
                               rot_cart=np.array([float(x) for x in p[2:11]]).reshape(3, 3)))
    return ops


def _report(ecut):
    qb = read_bands(XMLS[ecut])
    band = int(np.argmin(np.abs(qb.energies_ev[0].mean(axis=0) - qb.fermi_ev)))
    assert band == 44, f"ecut={ecut}: band nearest E_F moved to {band}"
    dE = qb.energies_ev[0, :, band] - qb.energies_ev[1, :, band]
    th = np.arctan2(qb.kpoints_tpiba[:, 1], qb.kpoints_tpiba[:, 0])
    dirs = np.stack([np.cos(th), np.sin(th), 0.093 * np.ones_like(th)], axis=1)
    basis = CircularBasis(LMAX)
    return noise_report(project(basis, dirs, dE), allowed_l(basis, _ops(), dirs))


def test_all_three_cutoffs_converged_to_the_same_compensated_state():
    for ecut in XMLS:
        qb = read_bands(XMLS[ecut])
        assert qb.kpoints_tpiba.shape[0] == 152
        assert qb.fermi_ev == pytest.approx(8.1168, abs=0.001), f"ecut={ecut}"


def test_the_allowed_harmonics_are_converged_in_ecutwfc():
    """P5's actual convergence question, answered: yes, 60 Ry is enough.

    Measured across 60/66/80 Ry, the useful form is the **absolute** spread,
    because that is what turns out to be constant: every harmonic from l=9 to
    l=45 varies by 0.012-0.025 meV regardless of whether the coefficient
    itself is 39 meV or 1 meV. l=3, the largest, moves 0.21 meV -- still only
    0.35% of itself.

    So the cutoff contributes an uncertainty of roughly **20 ueV per
    coefficient**, which is the same scale as P3's 27 ueV noise floor. Two
    independent numerical uncertainties landing on the same scale is a
    consistent picture, and it sets the honest error bar on any coefficient
    this project quotes: not the 1e-9 eV the fit reproduces internally, but
    a few tens of ueV.
    """
    amps = {e: _report(e).allowed for e in XMLS}
    for l in (9, 15, 21, 27, 33, 39, 45):
        vals = np.array([amps[e][l] for e in (60, 66, 80)])
        spread_meV = (vals.max() - vals.min()) / 1e-3
        assert spread_meV < 0.04, (
            f"l={l} spans {spread_meV:.4f} meV across cutoffs (was <=0.026)"
        )
    l3 = np.array([amps[e][3] for e in (60, 66, 80)])
    assert (l3.max() - l3.min()) / l3.mean() < 0.006, (
        f"l=3 spans {(l3.max()-l3.min())/l3.mean()*100:.2f}% (was 0.35%)"
    )


def test_the_noise_floor_ignores_an_eighty_fold_swing_in_the_predicted_mechanism():
    """The falsification, stated as a number.

    Predicted basis-set discontinuity: 1167.5 / 14.5 / 975.2 ppm at 60/66/80
    Ry -- an 80x range. Measured noise floor: 27.5 / 27.0 / 28.7 ueV -- a 6%
    range, with the *cleanest-predicted* cutoff (66 Ry) not the cleanest
    measured one. The even-l harmonics specifically, which are the only ones
    the discontinuity could touch since it is exactly even in l, sit at
    0.42 / 0.66 / 0.64 ueV -- again flat, and again with 66 Ry not winning.
    """
    reps = {e: _report(e) for e in XMLS}

    predicted = np.array([PREDICTED_PPM[e] for e in (60, 66, 80)])
    assert predicted.max() / predicted.min() > 50, "the prediction should span ~80x"

    floors = np.array([reps[e].floor for e in (60, 66, 80)])
    assert floors.max() / floors.min() < 1.15, (
        f"floor spans {floors.max()/floors.min():.2f}x -- it should be flat: {floors/1e-6} ueV"
    )

    even_max = np.array([
        max(reps[e].forbidden[l] for l in range(2, LMAX + 1, 2)) for e in (60, 66, 80)
    ])
    assert even_max.max() < 1e-6, "even-l harmonics should all be sub-ueV"
    # the decisive point: 66 Ry was predicted 80x cleaner and is not cleaner at all
    assert even_max[1] >= even_max[0], (
        "66 Ry was predicted to be the cleanest by far; if it now measures cleaner "
        "than 60 Ry, the falsification recorded in PLAN.md P4/P5 needs revisiting"
    )


def test_raising_the_cutoff_does_not_reduce_the_floor_either():
    """80 Ry uses 53% more plane waves than 60 Ry and buys nothing."""
    f60, f80 = _report(60).floor, _report(80).floor
    assert f80 >= f60 * 0.95, (
        f"80 Ry floor {f80/1e-6:.2f} ueV vs 60 Ry {f60/1e-6:.2f} ueV -- if the higher "
        "cutoff now helps materially, the recorded conclusion changes"
    )
