"""P19: the linear-response U is twice the one every production run used.

P12 measured what U = 3.0 eV costs but could not say whether 3.0 is right -- it
is an inherited choice and nothing here had tested it.  QE's hp.x answers that
directly (job 20796619, 2026-10-06):

    linear-response U     6.3858 eV (CrU) / 6.3880 eV (CrD)
    U used in production  3.0 eV
    ratio                 2.13x, a difference of 3.39 eV

The two magnetic sublattices agree to 0.002 eV, which symmetry requires; that
is the check that the run is self-consistent rather than drifting.  The Hubbard
projectors are ortho-atomic, matching production, and the SCF underneath ran at
U = 0.0001 eV, which is the correct linear-response procedure rather than a
typo.

**Direction is certain, magnitude is not.**  P12's scan is monotone -- every
coefficient falls as U rises (l=3 by 36% and l=45 by 61% from U=2 to U=5) -- so
the production coefficients are systematically too large.  But 6.39 eV lies
OUTSIDE the 2-5 eV range P12 scanned, so extrapolating "21% per eV" across
3.39 eV would be arithmetic, not measurement.  Closing it needs one 152-point
bands run at U ~ 6.4.

**Two reservations belong with the number**: the q-mesh is 2x2x1 (2 q-points),
which is coarse and hp.x is sensitive to it; and conv_thr_chi was loosened to
1e-4 because the response would not converge at 1e-6.  So 6.39 eV is not a
number for a paper.  It is a strong signal that the production U is badly low,
and since P12 showed U dominates every other uncertainty measured here, that
makes it the largest open issue in the pipeline.

**Why it sat unfinished for two weeks.**  Three of the four response matrices
existed since 2026-09-22; the collection step died in 25 seconds with
"AB.chi.pert_1.dat does not exist".  The cause was a MISSING STAGE, not a broken
one -- QE's HP workflow is

    dns.pert_i.q_j --(sum_pertq, over q per atom)--> chi.pert_i --(compute_hp)--> U

and the runs jumped straight from dns to compute_hp.  Adding sum_pertq cost 40
minutes of walltime, nearly all of it queueing: summing files already on disk is
close to free.  Several days of DFPT was one cheap step short of an answer.
"""

import re
from pathlib import Path

import pytest

DAT = Path(__file__).parent / "data" / "hp_u" / "AB.Hubbard_parameters.dat"
PRODUCTION_U = 3.0

pytestmark = pytest.mark.skipif(
    not DAT.exists(), reason="the hp.x Hubbard_parameters fixture is not present")


def _u_by_site():
    out = {}
    for line in DAT.read_text().splitlines():
        m = re.match(r"\s*(\d+)\s+\d+\s+(\S+)\s+-?\d+\s+\d+\s+\S+\s+3d\s+([\d.]+)\s*$", line)
        if m:
            out[int(m.group(1))] = (m.group(2), float(m.group(3)))
    return out


def test_all_four_cr_sites_are_reported():
    u = _u_by_site()
    assert sorted(u) == [1, 2, 3, 4]
    assert {lbl for lbl, _ in u.values()} == {"CrU", "CrD"}


def test_the_two_sublattices_agree_as_symmetry_requires():
    """The self-consistency check: a drifting run would not give this."""
    u = _u_by_site()
    vals = sorted({round(v, 4) for _, v in u.values()})
    assert len(vals) == 2, vals
    assert abs(vals[1] - vals[0]) < 0.005, vals


def test_the_linear_response_u_is_about_twice_production():
    u = _u_by_site()
    mean = sum(v for _, v in u.values()) / len(u)
    assert 6.3 < mean < 6.5, mean
    assert mean / PRODUCTION_U > 2.0


def test_the_value_lies_outside_the_range_p12_scanned():
    """So P12's 21%-per-eV cannot be extrapolated to it; say so, do not do it."""
    u = _u_by_site()
    mean = sum(v for _, v in u.values()) / len(u)
    p12_scanned = (2.0, 5.0)
    assert mean > p12_scanned[1], (mean, p12_scanned)
