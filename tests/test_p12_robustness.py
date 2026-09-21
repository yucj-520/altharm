"""Which conclusions survive the choice of U, and which do not.

Standard practice is to pick a reasonable U and check that it does not change
the conclusion.  P12's scan is what makes the second half of that sentence a
measurement rather than an assertion, so this file separates the two groups
explicitly.  It exists because the first pass at this separation (in
conversation, not in code) overstated three of the entries, and the overstated
versions are recorded here as the assertions that would fail.

**Robust over U = 2..5 eV**

- the allowed set (pure symmetry, no DFT input at all)
- every allowed harmonic being pure sin(l theta) -- worst cos/sin ratio 5e-4,
  and the cos amplitude below that run's own noise floor at every U
- the leading two harmonics being l = 3 then l = 9
- l = 45 measurable: 15x to 383x the floor
- lmax = 45 sufficing, and the fit staying well-conditioned: residual 0.6-1.9%
  of the signal, condition number 1.4 at every U

**Not robust**

- any individual coefficient's value: l=3 spans 43% over the range, ~21% per eV
- the full ordering beyond the leading two: l=15 moves from last (U=2) to third
- l=15 existing at all: 0.3 ueV at U=2, which is BELOW that run's floor
- and therefore "every allowed harmonic is above the floor", which is true at
  U = 3,4,5 and false at U = 2

**Not tested, and so claimed for neither group:** whether the QE/Wannier/defumat
agreement holds away from U = 3. All three routes were only ever compared there.
Asserting robustness for it would be an extrapolation; the test below records
the gap instead of hiding it.
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
RUNS = {2.0: "qe_cr2seas_u2", 3.0: "qe_cr2seas_sm12",
        4.0: "qe_cr2seas_u4", 5.0: "qe_cr2seas_u5"}
REPORT = HERE / "data" / "harmwave_cr2seas" / "harmwave_report.txt"

pytestmark = pytest.mark.skipif(
    not (all((HERE / "data" / d / "data-file-schema.xml").exists() for d in RUNS.values())
         and REPORT.exists()),
    reason="the Hubbard-U scan fixtures are not present",
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
def scan():
    out = {}
    for U, folder in RUNS.items():
        q = read_bands(HERE / "data" / folder / "data-file-schema.xml")
        th = np.arctan2(q.kpoints_tpiba[:, 1], q.kpoints_tpiba[:, 0])
        dirs = np.stack([np.cos(th), np.sin(th), 0.093 * np.ones_like(th)], axis=1)
        basis = CircularBasis(LMAX)
        f = project(basis, dirs, q.energies_ev[0, :, BAND] - q.energies_ev[1, :, BAND])
        floor = noise_report(f, allowed_l(basis, _ops(), dirs)).floor
        out[U] = (f.amplitude_by_l(), floor, dict(zip(f.labels, f.coeffs)), f)
    return out


# ---------------------------------------------------------------- robust

def test_pure_sine_holds_at_every_u(scan):
    """Checked per U this time.

    The first pass at this claim ran a loop with an empty body and printed the
    conclusion without computing it.  The conclusion happened to be right; the
    output was not evidence.
    """
    for U in RUNS:
        a, floor, c, _ = scan[U]
        for l in ALLOWED:
            cos, sin = abs(c[(l, "c")]), abs(c[(l, "s")])
            assert cos < 3 * floor, f"U={U}, l={l}: cos {cos*1e6:.2f} ueV vs floor {floor*1e6:.2f}"
            if sin > 10 * floor:
                assert cos / sin < 1e-3, f"U={U}, l={l}: cos/sin {cos/sin:.2e}"


def test_the_leading_two_harmonics_are_u_independent(scan):
    for U in RUNS:
        a = scan[U][0]
        assert sorted(ALLOWED, key=lambda l: -a[l])[:2] == [3, 9], f"U={U}"


def test_l45_is_measurable_at_every_u(scan):
    """15x to 383x the floor.  This is the honest form of 'high harmonics exist'."""
    margins = {}
    for U in RUNS:
        a, floor, _, _ = scan[U]
        margins[U] = a[45] / floor
        assert margins[U] > 10, f"U={U}: only {margins[U]:.1f}x the floor"
    assert margins[2.0] == pytest.approx(15, abs=3)
    assert margins[5.0] == pytest.approx(383, rel=0.15)


def test_the_fit_stays_well_posed_at_every_u(scan):
    """lmax=45 sufficing is a premise, and it does not depend on U."""
    for U in RUNS:
        _, _, _, f = scan[U]
        assert f.cond == pytest.approx(1.4, abs=0.1), f"U={U}"
        assert f.residual_rms / f.values_rms < 0.025, f"U={U}"


# ------------------------------------------------------------ NOT robust

def test_individual_coefficients_are_not_u_robust(scan):
    v = [scan[U][0][3] * 1e3 for U in (2.0, 3.0, 4.0, 5.0)]
    assert (max(v) - min(v)) / np.mean(v) == pytest.approx(0.43, abs=0.05)


def test_the_full_ordering_is_not_u_robust(scan):
    """l=15 goes from last at U=2 to third at U=5."""
    order2 = sorted(ALLOWED, key=lambda l: -scan[2.0][0][l])
    order5 = sorted(ALLOWED, key=lambda l: -scan[5.0][0][l])
    assert order2 != order5
    assert order2[-1] == 15 and order5[2] == 15


def test_not_every_allowed_harmonic_is_above_the_floor_at_every_u(scan):
    """The claim that had to be weakened.

    'l=45 is measurable at every U' is true.  'every allowed harmonic is
    measurable at every U' is false: at U=2 the weakest is l=15, which sits
    below that run's own floor.
    """
    a2, floor2, _, _ = scan[2.0]
    weakest = min(ALLOWED, key=lambda l: a2[l])
    assert weakest == 15
    assert a2[weakest] < floor2, "at U=2 the weakest allowed harmonic is below the floor"
    # while at U >= 3 the statement does hold
    for U in (3.0, 4.0, 5.0):
        a, floor, _, _ = scan[U]
        assert min(a[l] for l in ALLOWED) > 10 * floor, f"U={U}"


def test_the_floor_is_not_a_single_number(scan):
    """'micro-eV level' was loose: it is 1.3 to 82 ueV across the range."""
    floors = [scan[U][1] for U in RUNS]
    assert min(floors) == pytest.approx(1.3e-6, rel=0.15)
    assert max(floors) == pytest.approx(81.8e-6, rel=0.05)
    assert max(floors) / min(floors) > 40


def test_the_three_route_agreement_was_never_measured_away_from_u3():
    """Records a gap rather than asserting robustness for it.

    P7 compared QE, Wannier and defumat only at U = 3.0 eV.  The Wannier model
    is a fit to a DFT run at that U and cannot be re-evaluated at another one
    without refitting, so this cannot be closed inside this repository.
    """
    u3_fixtures = [HERE / "data" / "defumat_p7" / "defumat_full.pkl",
                   HERE / "data" / "qe_cr2seas_dense" / "data-file-schema.xml"]
    assert all(p.exists() for p in u3_fixtures)
    for U in (2.0, 4.0, 5.0):
        assert not (HERE / "data" / f"defumat_u{U:.0f}.pkl").exists(), (
            "if a defumat run at another U appears, extend P7 rather than "
            "assuming the agreement carries over"
        )
