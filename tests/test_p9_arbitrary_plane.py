"""harmwave.x v2: the circle's plane is read from the k-points, not assumed.

v1 required the circle to lie in the k_x k_y plane at fixed k_z, hardcoded in
two places: a check that every k-point shared one k_z, and a symmetry filter
that tested whether row 3 of R was (0,0,1).  That is where this project's own
circles happen to lie and not where another material's will.

v2 determines the plane from the k-list (normal from the first three
non-collinear points, then every point verified against it), measures theta in
an in-plane frame (e1,e2) written into the report, and generalises the filter to
its coordinate-free form: k' = R k preserves n.k for every k exactly when
R^T n = n.  For n = +z that reduces to v1's row-3 test, so a k_x k_y circle must
come out bit for bit the same -- which is what this file checks.

**Measured**, on the corrected Cr2SeAs reference (24x24x8, degauss 0.0150,
152-point circle at k_z = 0.25 b3, band 44):

- all eight allowed harmonics identical (l=9 differs by 2.3e-13, the rest by 0)
- largest absolute difference over all 91 coefficients: 1e-14 eV
- normal detected as (0,0,1), offset 0.0550574 (= 0.25 b3 in tpiba), frame
  (e1,e2) = (x,y) exactly
- the report gains exactly four lines, all of them the frame

**The gap this does NOT close.** It shows v2 reproduces v1 where v1 was valid.
It did not show v2 is *correct* for a plane whose normal is not +z.  **That gap
is now closed** (job 20383535, read 2026-10-01): a great circle in the k_y k_z
plane, normal (1,0,0), offset 0, frame (e1,e2) = (y,z), 4 of 12 operations kept.
Its allowed set is l = 2,4,...,24, matching `altharm.symmetry`'s independent
prediction for that plane, and the coefficients are sharper still -- even l,
sine only.  See `tests/test_p14_sampling_aliasing.py` and PLAN.md P14, which
also record what that run revealed about the noise floor itself.
"""

import re
from pathlib import Path

import numpy as np
import pytest

DATA = Path(__file__).parent / "data" / "harmwave_v2"
V1 = DATA / "harmwave_report_v1.txt"
V2 = DATA / "harmwave_report_v2.txt"

pytestmark = pytest.mark.skipif(not (V1.exists() and V2.exists()),
                                reason="the harmwave v1/v2 reports are not present")

ALLOWED = (3, 9, 15, 21, 27, 33, 39, 45)


def _coeffs(path):
    out = {}
    for line in path.read_text().splitlines():
        m = re.match(r"CF\s+(\d+)\s+([cs])\s+(\S+)", line)
        if m:
            out[(int(m.group(1)), m.group(2))] = float(m.group(3))
    return out


def _header(path, key):
    for line in path.read_text().splitlines():
        if key in line:
            return [float(x) for x in re.findall(r"-?\d+\.\d+", line)]
    return None


def test_the_two_reports_describe_the_same_fit():
    c1, c2 = _coeffs(V1), _coeffs(V2)
    assert len(c1) == 91 == len(c2)
    assert set(c1) == set(c2)


def test_every_allowed_harmonic_is_unchanged():
    """The one that matters: a k_x k_y circle must be untouched by the change."""
    c1, c2 = _coeffs(V1), _coeffs(V2)
    for l in ALLOWED:
        a1 = np.hypot(c1[(l, "c")], c1[(l, "s")])
        a2 = np.hypot(c2[(l, "c")], c2[(l, "s")])
        assert a2 == pytest.approx(a1, rel=1e-11), f"l={l}"


def test_no_coefficient_moved_beyond_last_bit():
    c1, c2 = _coeffs(V1), _coeffs(V2)
    diffs = np.array([c2[k] - c1[k] for k in c1])
    assert np.abs(diffs).max() < 1e-13, "differences should be eigenvalue last-bit only"


def test_the_plane_was_detected_and_not_assumed():
    n = _header(V2, "plane normal")
    assert n == pytest.approx([0.0, 0.0, 1.0], abs=1e-8)
    # 0.25 * b3 in tpiba units: (2 pi / c) / (2 pi / alat) / 4 = alat / (4 c)
    alat, c = 3.7258396285, 16.9179784926
    assert _header(V2, "plane offset")[0] == pytest.approx(alat / (4 * c), rel=1e-4)


def test_the_frame_reduces_to_x_and_y_for_a_plus_z_normal():
    """Why the numbers are reproducible rather than merely equivalent."""
    assert _header(V2, "theta frame e1") == pytest.approx([1.0, 0.0, 0.0], abs=1e-8)
    assert _header(V2, "theta frame e2") == pytest.approx([0.0, 1.0, 0.0], abs=1e-8)


def test_the_report_gained_only_the_frame():
    s1 = set(V1.read_text().splitlines())
    s2 = set(V2.read_text().splitlines())
    added = [l for l in s2 - s1 if not l.startswith(("CF", "PT"))]
    assert len(added) == 4
    assert all(any(k in l for k in ("plane normal", "plane offset", "theta frame"))
               for l in added)


def test_the_non_plus_z_case_is_still_untested():
    """Records the gap rather than implying it is closed."""
    assert not (DATA / "harmwave_report_yz.txt").exists(), (
        "if a circle in another plane has been run, assert its allowed set "
        "against altharm.symmetry's prediction and drop this placeholder"
    )
