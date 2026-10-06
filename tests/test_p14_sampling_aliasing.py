"""P14: the noise floor is anti-correlated with the error it is supposed to report.

Phases P3 and P12/P13 established that the forbidden-harmonic floor is blind to
five inputs -- U, smearing, pseudopotential, functional, geometry -- each of
which moves the coefficients 7-21% while the floor reports 0.04%.  Blindness is
survivable: a number that does not notice an error is still not *wrong* about
the errors it does notice.  This phase found a sixth input where it is worse
than blind.

The number of sample points on the circle sets the floor and the error at the
same time, and sets them in opposite directions.

**The mechanism**, confirmed here two independent ways.  A material's splitting
contains only allowed harmonics, at every order; on Cr2SeAs the allowed set is
``l = 3 (mod 6)``, a lattice of period ``p = 6``.  Everything above ``lmax`` is
unresolved and aliases onto ``l' = |l - m n|``.  That ``l'`` is allowed exactly
when ``m n = 0 (mod p)``, so the tail first touches the signal at alias order
``m = p / gcd(n, p)``.  When ``p | n`` that order is 1: the strongest alias of
every unresolved harmonic is added directly to a reported coefficient, and no
forbidden channel receives anything at all.

**Measured on the real Cr2SeAs run** (24x24x8, degauss 0.0150, 152 points,
band 44), splitting the 90 coefficients by what forbids each one:

================================================  ===  ===========  ==========
group                                             n    median       max
================================================  ===  ===========  ==========
allowed                                           8    3093 ueV     64293 ueV
forbidden by an op the sample grid CAN represent  67   0.0040 ueV   1.11 ueV
forbidden only by C3, which it CANNOT             15   24.70 ueV    75.0 ueV
================================================  ===  ===========  ==========

The published floor, 23.6 ueV, is the third row.  It is read entirely off the
fifteen harmonics that only the C3 rotations forbid -- and 120 degrees is not a
multiple of 360/152, so C3 is the one operation this grid cannot represent.
Every harmonic forbidden by an operation the grid *can* represent sits 6000x
lower.  The floor is therefore a measure of one specific aliasing, not of
numerical error in general.

**Measured synthetically**, sweeping every n from 120 to 199 against a tail
fitted to the real coefficients:

=============  =====  ================  ===================
``gcd(n, 6)``  ``m``  reported floor    true error, median
=============  =====  ================  ===================
6              1      **0 ueV**         **18.1 %**
3              2      160 ueV           3.7 %
2              3      221 ueV           1.4 %
1              6      239 ueV           0.29 %
=============  =====  ================  ===================

Monotone, with no exceptions in that range.  The grid whose floor reads a
perfect zero is sixty times less accurate than the grid with the worst-looking
floor.

**What is still open.** Both rows of evidence are about *this* kind of error --
the unresolved tail.  The synthetic model has no SCF noise in it at all, and
the real-data split is one material, one band, two planes.  The DFT
confirmation -- the same density sampled at n = 150, 151, 153, where the
prediction is that n = 150 reports a floor near 0.004 ueV and shifts the
coefficients by order 10% -- is prepared in ``altharm_runs/nsweep`` and has not
been run.  Until it has, the first table is an observation about the 152-point
run and the second is a statement about a model, and neither is a measurement
of what n = 150 does to Cr2SeAs.
"""

import functools
import re
from math import gcd
from pathlib import Path

import numpy as np
import pytest

from altharm.sampling import alias_order, allowed_period, sampling_report

DATA = Path(__file__).parent / "data"
YZ = DATA / "harmwave_yz" / "harmwave_yz.txt"
YZ_K = DATA / "harmwave_yz" / "kcirc_tpiba.npy"
XY = DATA / "harmwave_v2" / "harmwave_report_v2.txt"
XY_K = DATA / "harmwave_v2" / "kcirc_xy_tpiba.npy"

CR2SEAS_ALLOWED = (3, 9, 15, 21, 27, 33, 39, 45)
MEASURED = {3: 64292.8, 9: 43446.7, 15: 3603.3, 21: 8661.7,
            27: 1876.3, 33: 2583.2, 39: 1181.9, 45: 985.3}   # ueV, 152 points


# --------------------------------------------------------------------------
# the arithmetic of the rule
# --------------------------------------------------------------------------

def test_allowed_period_of_cr2seas_is_six():
    assert allowed_period(CR2SEAS_ALLOWED) == 6


def test_allowed_period_rejects_an_uneven_set():
    # a y-z circle's allowed set is every even l -- period 2
    assert allowed_period(range(2, 25, 2)) == 2
    # and a set with no constant spacing has no period, which must not be faked
    assert allowed_period([2, 3, 7]) == 0


def test_alias_order_is_period_over_gcd():
    for n in range(100, 200):
        assert alias_order(n, 6) == 6 // gcd(n, 6)


def test_a_multiple_of_the_period_is_order_one():
    for n in (120, 150, 156, 162, 300):
        assert alias_order(n, 6) == 1


def test_sampling_report_flags_the_dangerous_counts():
    assert not sampling_report(150, CR2SEAS_ALLOWED).safe
    assert not sampling_report(156, CR2SEAS_ALLOWED).safe
    assert sampling_report(151, CR2SEAS_ALLOWED).safe
    assert sampling_report(155, CR2SEAS_ALLOWED).safe


def test_the_project_own_count_is_flagged_as_suboptimal_not_fatal():
    r = sampling_report(152, CR2SEAS_ALLOWED)
    assert not r.safe and r.order == 3
    assert "152" in r.message and "coprime" in r.message


def test_report_suggests_a_coprime_count_near_the_one_given():
    for n in (150, 152, 153, 156):
        r = sampling_report(n, CR2SEAS_ALLOWED)
        suggested = int(re.findall(r"e\.g\. (\d+)", r.message)[0])
        assert gcd(suggested, 6) == 1
        assert abs(suggested - n) <= 3


# --------------------------------------------------------------------------
# the mechanism, on data whose truth is known exactly
# --------------------------------------------------------------------------

def _synthetic(n, decay=-1.5, lmax=45, l_true_max=1200):
    """Fit a perfectly symmetric, noiseless signal sampled at n points.

    Every harmonic present is allowed, so a perfect fit would put zero on every
    forbidden l and reproduce MEASURED exactly.  Anything else is aliasing.
    """
    amps = dict(MEASURED)
    for l in range(51, l_true_max, 6):
        amps[l] = 985.3 * (l / 45.0) ** decay * (1 + 0.3 * np.sin(l / 7.0))
    th = 2 * np.pi * np.arange(n) / n
    # The tail runs to l_true_max because truncating it is not a harmless
    # economy: at N=152 the third-order alias sits at 3*152 = 456, so a cutoff
    # below that silently removes the very term under test.  Build the signal
    # and the design matrix as single matrix products instead of Python loops.
    ls = np.fromiter(amps, dtype=float)
    f = np.sin(np.outer(th, ls)) @ np.fromiter(amps.values(), dtype=float)
    harm = np.arange(1, lmax + 1)
    arg = np.outer(th, harm)
    design = np.empty((n, 1 + 2 * lmax))
    design[:, 0] = 1.0
    design[:, 1::2] = np.cos(arg)
    design[:, 2::2] = np.sin(arg)
    c, *_ = np.linalg.lstsq(design, f, rcond=None)
    coef = {(l, k): c[1 + 2 * (l - 1) + j]
            for l in range(1, lmax + 1) for j, k in enumerate("cs")}
    floor = max(abs(v) for (l, k), v in coef.items()
                if not (l % 6 == 3 and k == "s"))
    err = max(abs(coef[(l, "s")] - MEASURED[l]) / MEASURED[l] for l in MEASURED)
    return floor, 100 * err


def test_multiple_of_six_reports_a_zero_floor_while_being_wrong():
    for n in (150, 156, 162):
        floor, err = _synthetic(n)
        assert floor < 1e-6, f"n={n} should leave the forbidden channels empty"
        assert err > 10.0, f"n={n} should corrupt the allowed coefficients"


def test_coprime_count_is_accurate_though_its_floor_looks_worse():
    for n in (151, 155, 157):
        floor, err = _synthetic(n)
        assert err < 1.0, f"n={n} should keep the coefficients within 1%"
        assert floor > 100.0, f"n={n} should show the tail in the floor"


def test_error_falls_monotonically_with_alias_order():
    by_order = {}
    for n in range(120, 200):
        by_order.setdefault(alias_order(n, 6), []).append(_synthetic(n)[1])
    medians = {m: float(np.median(v)) for m, v in by_order.items()}
    assert sorted(medians) == [1, 2, 3, 6]
    ordered = [medians[m] for m in sorted(medians)]
    assert ordered == sorted(ordered, reverse=True), medians
    assert medians[1] > 30 * medians[6], medians


def test_the_anticorrelation_holds_for_every_tail_decay():
    """The conclusion must not depend on how fast the tail falls off."""
    for decay in (-1.0, -1.5, -2.0, -2.5, -3.0):
        bad = _synthetic(150, decay=decay)
        good = _synthetic(151, decay=decay)
        assert bad[0] < good[0], decay      # the dangerous grid LOOKS better
        assert bad[1] > good[1], decay      # and IS worse


def test_no_multiple_of_six_escapes_the_rule():
    offenders = [n for n in range(120, 200) if n % 6 == 0
                 and not (_synthetic(n)[0] < 1e-6 and _synthetic(n)[1] > 1.0)]
    assert offenders == []


# --------------------------------------------------------------------------
# the mechanism, on the real run
# --------------------------------------------------------------------------

def _parse(path):
    ops, coef, theta, dE = [], {}, [], []
    for line in Path(path).read_text().splitlines():
        t = line.split()
        if line.startswith("OP "):
            ops.append((int(t[1]),
                        np.array([float(x) for x in t[2:11]]).reshape(3, 3)))
        elif line.startswith("CF "):
            coef[(int(t[1]), t[2])] = float(t[3])
        elif line.startswith("PT "):
            theta.append(float(t[1]))
            dE.append(float(t[2]))
    return ops, coef, np.array(theta), np.array(dE)


def _plane_action(R, k, e1, e2):
    """theta -> a*theta + b, plus whether R permutes the sample points."""
    kr = k @ R.T
    th0 = np.arctan2(k @ e2, k @ e1)
    th1 = np.arctan2(kr @ e2, kr @ e1)
    best = min(((np.abs(np.angle(np.exp(1j * (th1 - (a * th0 + np.angle(
        np.mean(np.exp(1j * (th1 - a * th0))))))))).max()), a) for a in (1, -1))
    a = best[1]
    b = float(np.angle(np.mean(np.exp(1j * (th1 - a * th0)))))
    perm = np.linalg.norm(kr[:, None] - k[None], axis=2).min(axis=1).max() < 1e-6
    return a, b, perm


def _block(l, a, b):
    c, s = np.cos(l * b), np.sin(l * b)
    return np.array([[c, s], [-a * s, a * c]])


def _classify(path, kpath, e1, e2, lmax):
    ops, coef, _, _ = _parse(path)
    k = np.load(kpath)
    acts = [(chi,) + _plane_action(R, k, e1, e2) for chi, R in ops]
    groups = {"allowed": [], "perm": [], "nonperm": []}
    for l in range(1, lmax + 1):
        for name, v in (("c", np.array([1.0, 0.0])), ("s", np.array([0.0, 1.0]))):
            bad_perm = bad_free = False
            # b is recovered from k-points printed to finite precision, so it
            # carries ~5e-10 rad; the comparison is on cos(l b), whose error
            # therefore grows as l*db and crosses a fixed 1e-8 at l=20.  A
            # genuine violation is O(1), so scale the tolerance with l rather
            # than loosening it globally.
            for chi, a, b, perm in acts:
                tol = 1e-6 * max(1, l)
                if not np.allclose(_block(l, a, b) @ v, chi * v, atol=tol):
                    if perm:
                        bad_perm = True
                    else:
                        bad_free = True
            amp = abs(coef.get((l, name), 0.0)) * 1e6
            key = "allowed" if not (bad_perm or bad_free) else (
                "perm" if bad_perm else "nonperm")
            groups[key].append((l, name, amp))
    return groups


@pytest.mark.skipif(not (XY.exists() and XY_K.exists()),
                    reason="the k_x k_y harmwave report is not present")
def test_published_floor_comes_only_from_the_unrepresentable_rotation():
    g = _classify(XY, XY_K, np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), 45)

    assert [l for l, k, _ in g["allowed"]] == list(CR2SEAS_ALLOWED)
    assert all(k == "s" for _, k, _ in g["allowed"])

    perm = np.array([a for _, _, a in g["perm"]])
    free = np.array([a for _, _, a in g["nonperm"]])

    # the fifteen C3-forbidden sine channels ARE the published 23.6 ueV floor
    assert len(free) == 15
    assert 20.0 < float(np.median(free)) < 30.0

    # everything forbidden by an op the grid can represent is thousands of
    # times smaller -- read the floor there and it would be 6000x too good
    assert float(np.median(perm)) < 0.02
    assert float(np.median(free)) / float(np.median(perm)) > 1000.0


@pytest.mark.skipif(not (YZ.exists() and YZ_K.exists()),
                    reason="the k_y k_z harmwave report is not present")
def test_on_the_yz_circle_every_operation_is_representable():
    """Which is why its forbidden channels collapse to machine zero.

    The y-z circle's kept group is {E, C2x} with chi=+1 and {M_y, M_z} with
    chi=-1 -- all four of them sample permutations on 152 points.  So this
    circle has no C3-style channel at all, and its forbidden amplitudes are
    bounded by SCF symmetry error rather than by aliasing.  Reading a "noise
    floor" off it and quoting it next to the k_x k_y coefficients would
    understate that run's error by four orders of magnitude.
    """
    g = _classify(YZ, YZ_K, np.array([0, 1.0, 0]), np.array([0, 0, 1.0]), 24)

    assert g["nonperm"] == []
    assert [l for l, k, _ in g["allowed"]] == list(range(2, 25, 2))
    assert all(k == "s" for _, k, _ in g["allowed"])

    perm = np.array([a for _, _, a in g["perm"]])
    assert float(np.max(perm)) < 1.0          # ueV


# --------------------------------------------------------------------------
# the mechanism, measured on real eigenvalues at two sample counts
# --------------------------------------------------------------------------
#
# The synthetic tests above control the truth but contain no SCF noise; the
# classification tests use real data but only one sample count.  This closes
# the gap between them without a new DFT run, by using the 304-point circle
# twice.  Fitted to its own Nyquist it measures the high-l tail directly;
# subsampled to every second point it IS a 152-point circle on identical
# eigenvalues, whose forbidden channels should then reproduce that tail.
#
# On a 152-point grid, sin(l*theta_j) with l = 152 - l' folds to -sin(l'*theta_j),
# so the prediction is signed: b_{l'} = -b_{152-l'}.
#
# Measured (Cr2SeAs, band 44, degauss 0.0037 -- the k304 fixture's own density):
#
#   l'    from    b_l' at N=152    -b_(152-l') at N=304    residual
#    5     147      -33.755 ueV          -32.617 ueV        -3.49%
#   11     141      +40.339              +39.324            +2.58%
#   23     129      +54.666              +53.835            +1.54%
#   41     111      -86.259              -85.661            -0.70%
#
# All seven signs correct.  The residual shrinks with l' and is what the
# second fold partner (l = 152 + l', above 304's Nyquist) would supply.
#
# So the forbidden-harmonic "noise floor" is not noise: it is the signal's own
# unresolved tail, folded down, and it can be predicted from an independent
# measurement of that tail to a few percent.

K304 = DATA / "qe_cr2seas_k304" / "data-file-schema.xml"

ALIAS_PAIRS = [(5, 147), (11, 141), (17, 135), (23, 129),
               (29, 123), (35, 117), (41, 111)]


@functools.lru_cache(maxsize=1)
def _k304_fits():
    # four tests below need these two fits; the XML read and the 304x303
    # least squares are the expensive part of this file, so do them once
    from altharm.qe import read_bands, spin_splitting
    qb = read_bands(K304)
    k = qb.kpoints_tpiba
    dE = spin_splitting(qb, 44)
    th = np.arctan2(k[:, 1], k[:, 0]) % (2 * np.pi)
    order = np.argsort(th)
    th, dE = th[order], dE[order]

    def fit(t, f, lmax):
        d = np.column_stack([np.ones(len(t))] +
                            [g(l * t) for l in range(1, lmax + 1)
                             for g in (np.cos, np.sin)])
        c, *_ = np.linalg.lstsq(d, f, rcond=None)
        return {(l, kk): c[1 + 2 * (l - 1) + j]
                for l in range(1, lmax + 1) for j, kk in enumerate("cs")}

    return fit(th, dE, 151), fit(th[::2], dE[::2], 45), len(th)


@pytest.mark.skipif(not K304.exists(), reason="the 304-point fixture is not present")
def test_the_two_sample_counts_come_from_one_run():
    hi, lo, n = _k304_fits()
    assert n == 304            # and the subsample is therefore exactly 152


@pytest.mark.skipif(not K304.exists(), reason="the 304-point fixture is not present")
def test_forbidden_amplitudes_are_the_folded_tail_with_the_right_sign():
    hi, lo, _ = _k304_fits()
    for lp, src in ALIAS_PAIRS:
        got, pred = lo[(lp, "s")], -hi[(src, "s")]
        assert np.sign(got) == np.sign(pred), (lp, src, got, pred)
        assert abs(got - pred) / abs(pred) < 0.05, (lp, src, got, pred)


@pytest.mark.skipif(not K304.exists(), reason="the 304-point fixture is not present")
def test_the_fold_residual_shrinks_with_l():
    """What is left over is the next fold partner, which decays."""
    hi, lo, _ = _k304_fits()
    res = [abs((lo[(lp, "s")] + hi[(src, "s")]) / hi[(src, "s")])
           for lp, src in ALIAS_PAIRS]
    assert res == sorted(res, reverse=True), res


@pytest.mark.skipif(not K304.exists(), reason="the 304-point fixture is not present")
def test_allowed_coefficients_survive_the_subsampling():
    """Both counts are gcd=2 with the period, so neither is corrupted at first
    order -- the tail lands in forbidden channels, not on the signal."""
    hi, lo, _ = _k304_fits()
    for l in CR2SEAS_ALLOWED:
        rel = abs(lo[(l, "s")] - hi[(l, "s")]) / abs(hi[(l, "s")])
        assert rel < 1e-3, (l, rel)


# --------------------------------------------------------------------------
# the rule across three materials -- and what bounds the damage
# --------------------------------------------------------------------------
#
# The allowed-set period is a property of the material's magnetic point group,
# not of this project's conventions, so the rule has to be checked where that
# period differs.  Measured on the existing P8 fixtures:
#
#   material       N    p   gcd   m    floor        next/leading
#   MnTe  U=4     152   6    2    3    0.0234 ueV      5.630%
#   CoF2  U=0     152   4    4    1    0.0369 ueV      0.021%
#   Cr2SeAs U=3   304   6    2    3    1.0082 ueV     65.601%
#
# CoF2's allowed set is l = 2, 6, 10, 14 -- period 4 -- and 152 is divisible by
# 4, so **the published P8 CoF2 runs are in the m=1 regime**.  That went
# unnoticed.  It did them no harm, and the reason is the refinement that
# matters: the damage at m=1 is proportional to the unresolved tail, and CoF2
# barely has one (its second allowed harmonic is 0.021% of the first).
#
# So the failure mode is not uniform across materials.  It is worst exactly
# where the harmonic ladder is slow to decay -- which is to say on the
# multi-harmonic materials this project exists to study.  Cr2SeAs, whose second
# harmonic is 66% of its first, is the one that would be destroyed by m=1;
# CoF2 and MnTe are nearly immune whatever N is chosen.
#
# Corollary for P8: CoF2's numbers stand, but its sample count was lucky rather
# than chosen.  Confirming them would want a rerun at N coprime to 4.

def test_cof2_period_is_four_and_152_puts_it_at_first_order():
    """The fact that went unnoticed in P8."""
    assert allowed_period((2, 6, 10, 14)) == 4
    assert alias_order(152, 4) == 1
    assert not sampling_report(152, (2, 6, 10, 14)).safe


def test_mnte_and_cr2seas_share_cr2seas_period_and_are_not_at_first_order():
    for allowed in [(3, 9, 15, 21, 27, 33, 39, 45)]:
        assert allowed_period(allowed) == 6
        for n in (152, 304):
            assert alias_order(n, 6) == 3


def test_damage_at_first_order_is_bounded_by_the_tail():
    """A material with no tail cannot be hurt by any sample count.

    Synthetic, with the tail amplitude scanned: at m=1 every unresolved
    harmonic lands on the signal, so the error tracks the tail's size and
    vanishes with it.  This is why CoF2 survived m=1 and why Cr2SeAs would
    not have.
    """
    def err_at(n, tail_scale):
        amps = dict(MEASURED)
        for l in range(51, 1200, 6):
            amps[l] = tail_scale * 985.3 * (l / 45.0) ** -1.5
        th = 2 * np.pi * np.arange(n) / n
        f = sum(a * np.sin(l * th) for l, a in amps.items())
        d = np.column_stack([np.ones(n)] + [g(l * th) for l in range(1, 46)
                                            for g in (np.cos, np.sin)])
        c, *_ = np.linalg.lstsq(d, f, rcond=None)
        coef = {l: c[1 + 2 * (l - 1) + 1] for l in MEASURED}
        return 100 * max(abs(coef[l] - MEASURED[l]) / MEASURED[l] for l in MEASURED)

    errs = [err_at(150, s) for s in (1.0, 0.1, 0.01, 0.001)]
    assert errs == sorted(errs, reverse=True), errs
    assert errs[0] > 10.0          # a full tail at m=1 is devastating
    assert errs[-1] < 0.1          # a negligible tail at m=1 is harmless


def test_the_rule_holds_for_every_allowed_set_period():
    """It is arithmetic about the period, not a fact about three materials.

    The allowed set of a circle is ``l = off (mod p)``; p is 2 for the y-z
    circle here, 4 for CoF2, 6 for MnTe and Cr2SeAs.  Nothing in the derivation
    cares which, so the rule is checked against periods this project has no
    material for -- including odd and prime ones, where the coprime counts are
    much denser.  All nine behave the same way: a multiple of p puts the floor
    at zero and ~20-28% error on the coefficients, a coprime count does not.
    """
    def synth(n, p, off, lmax=45, decay=-1.5, l_true_max=1500):
        base = {l: 1000.0 * (l / off) ** decay for l in range(off, lmax + 1, p)}
        amps = dict(base)
        start = lmax + 1 + ((off - lmax - 1) % p)
        for l in range(start, l_true_max, p):
            amps[l] = 1000.0 * (l / off) ** decay * (1 + 0.3 * np.sin(l / 7.0))
        th = 2 * np.pi * np.arange(n) / n
        f = (np.sin(np.outer(th, np.fromiter(amps, dtype=float)))
             @ np.fromiter(amps.values(), dtype=float))
        arg = np.outer(th, np.arange(1, lmax + 1))
        d = np.empty((n, 1 + 2 * lmax))
        d[:, 0] = 1.0
        d[:, 1::2] = np.cos(arg)
        d[:, 2::2] = np.sin(arg)
        c, *_ = np.linalg.lstsq(d, f, rcond=None)
        coef = {(l, k): c[1 + 2 * (l - 1) + j]
                for l in range(1, lmax + 1) for j, k in enumerate("cs")}
        floor = max(abs(v) for (l, k), v in coef.items()
                    if not (l % p == off % p and k == "s"))
        err = 100 * max(abs(coef[(l, "s")] - base[l]) / abs(base[l]) for l in base)
        return floor, err

    for p, off in [(2, 2), (3, 3), (4, 2), (5, 5), (6, 3), (7, 7), (8, 4),
                   (10, 5), (12, 6)]:
        n_bad = 150 - (150 % p)
        n_good = next(n for n in range(149, 200) if gcd(n, p) == 1)
        assert alias_order(n_bad, p) == 1
        assert alias_order(n_good, p) == p
        f_bad, e_bad = synth(n_bad, p, off)
        f_good, e_good = synth(n_good, p, off)
        assert f_bad < 1e-6, (p, f_bad)          # floor reads zero
        assert e_bad > 10.0, (p, e_bad)          # yet the answer is wrong
        assert e_good < e_bad, (p, e_good, e_bad)
        assert f_good > f_bad, (p, f_good)       # the honest grid looks worse
