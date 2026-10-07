"""How many points to put on the circle -- and why the obvious answer is wrong.

Every other knob in this project is checked by the forbidden-harmonic floor:
fit all ``l`` up to ``lmax``, and the amplitude landing on an ``l`` that
symmetry forbids is the numerical error, measured on the same footing as the
signal.  The number of sample points is the one knob the floor cannot check,
because the floor is *computed from* the sampling.  Worse than uninformative,
it is anti-correlated: the sample count that makes the floor read exactly zero
is the one that corrupts the reported coefficients most.

The mechanism is aliasing with a parity.  A material's spin splitting contains
only symmetry-allowed harmonics, at every order -- on Cr2Se2As the allowed set
is ``l = 3, 9, 15, ...``, i.e. ``l = 3 (mod 6)``, so the allowed set is a
lattice of period ``p = 6``.  The part of that ladder above ``lmax`` is not
resolved; on ``n`` sample points it aliases down onto ``l' = |l - m n|`` for
integer ``m``.  Whether ``l'`` is itself allowed decides everything:

* ``l' = l - m n`` is allowed exactly when ``m n = 0 (mod p)``, so the lowest
  alias order that reaches the allowed set is ``m = p / gcd(n, p)``.
* At ``m = 1`` -- which is to say whenever ``p`` divides ``n`` -- the *first*
  and strongest alias of every unresolved harmonic lands on an allowed ``l'``.
  It is added straight to a reported coefficient, and **not one forbidden
  channel receives anything**.  The floor reads zero and the answer is wrong.
* When ``gcd(n, p)`` is small the first alias lands in a forbidden channel,
  where it is visible as the floor, and the allowed coefficients are only
  reached at order ``m``, by which point the tail has decayed.

Measured on a tail fitted to Cr2Se2As (``tests/test_p14_sampling_aliasing.py``),
sweeping every ``n`` from 120 to 199:

===============  =====  ==================  ====================
``gcd(n, 6)``    ``m``  reported floor      true error, median
===============  =====  ==================  ====================
6                1      **0 ueV**           **18 %**
3                2      160 ueV             3.7 %
2                3      221 ueV             1.4 %
1                6      239 ueV             0.29 %
===============  =====  ==================  ====================

So the rule is to pick ``n`` **coprime to the period of the allowed set**, and
never a multiple of it.  This project's own runs used ``n = 152``, which is
``gcd = 2``: not the catastrophic case, but a factor of five worse than a
coprime count would have been, and chosen by accident rather than on purpose.

**How much a bad ``n`` costs depends on the material.**  At ``m = 1`` the error
is the size of the unresolved tail, so a nearly single-wave altermagnet is
almost immune and a multi-harmonic one is not.  The P8 CoF2 runs are in fact at
``m = 1`` -- its allowed set is ``l = 2, 6, 10, 14``, period 4, and 152 is
divisible by 4 -- and came to no harm, because its second allowed harmonic is
0.021% of its first.  Cr2Se2As, at 66%, would have been destroyed.  The failure
mode therefore bites hardest on exactly the materials this project is for.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd
from typing import Iterable, Sequence


def allowed_period(allowed_l: Iterable[int]) -> int:
    """The period of the lattice of allowed ``l``.

    The allowed set produced by :func:`altharm.symmetry.allowed_l` is, for a
    circle with a rotation axis through it, an arithmetic progression -- the
    spacing of that progression is what aliasing has to be measured against.
    Returns ``0`` when the set is too short or not evenly spaced, which is the
    caller's signal that the cheap rule below does not apply and the sampling
    has to be checked by brute force (run two coprime ``n`` and compare).
    """
    ls = sorted({int(l) for l in allowed_l if int(l) > 0})
    if len(ls) < 2:
        return 0
    steps = {b - a for a, b in zip(ls, ls[1:])}
    return steps.pop() if len(steps) == 1 else 0


def alias_order(n_points: int, period: int) -> int:
    """Lowest alias order at which the unresolved tail reaches an allowed ``l``.

    ``1`` is the dangerous value: the first alias lands on the signal itself.
    Larger is better, and ``period`` (the maximum) means ``n_points`` is
    coprime to the allowed-set spacing.
    """
    if n_points <= 0 or period <= 0:
        raise ValueError("n_points and period must be positive")
    return period // gcd(n_points, period)


@dataclass(frozen=True)
class SamplingReport:
    """What a given sample count does to the tail, and whether it is safe."""

    n_points: int
    period: int
    order: int
    safe: bool
    message: str

    def __str__(self) -> str:          # pragma: no cover - presentation only
        return self.message


def sampling_report(n_points: int, allowed_l: Sequence[int]) -> SamplingReport:
    """Judge a circle's sample count against the allowed set it will be fitted to.

    This is deliberately not folded into the fit.  A fit that silently repaired
    the sample count would hide exactly the failure this module exists to make
    visible, and the only honest repair is to recompute the eigenvalues on a
    different circle -- which costs a DFT run and is the user's call.
    """
    period = allowed_period(allowed_l)
    if period == 0:
        return SamplingReport(
            n_points, 0, 0, False,
            f"n = {n_points}: the allowed set {sorted(allowed_l)} is not evenly "
            f"spaced, so the alias-order rule does not apply; check the sampling "
            f"by running two counts whose difference is not a factor of the "
            f"splitting's angular period and comparing the coefficients.",
        )
    order = alias_order(n_points, period)
    if order == 1:
        return SamplingReport(
            n_points, period, order, False,
            f"n = {n_points} is a multiple of the allowed-set period p = {period}. "
            f"Every unresolved harmonic above lmax aliases straight onto an "
            f"ALLOWED l at first order, so the forbidden-harmonic floor will read "
            f"essentially zero while the reported coefficients carry the whole "
            f"tail as error. How bad that is scales with the tail: on a material "
            f"whose harmonic ladder decays slowly (Cr2Se2As, 2nd/1st = 66%) it is "
            f"~18%; on one that is nearly single-wave (CoF2, 0.021%) it is "
            f"negligible. Use a count coprime to {period}, "
            f"e.g. {_nearest_coprime(n_points, period)}.",
        )
    if order == period:
        return SamplingReport(
            n_points, period, order, True,
            f"n = {n_points} is coprime to the allowed-set period p = {period}; "
            f"the tail reaches the allowed coefficients only at alias order "
            f"{order}, and its first alias is visible in the forbidden channels.",
        )
    return SamplingReport(
        n_points, period, order, False,
        f"n = {n_points} shares the factor {gcd(n_points, period)} with the "
        f"allowed-set period p = {period}: the tail reaches the allowed "
        f"coefficients at alias order {order}. Usable, but a count coprime to "
        f"{period} (e.g. {_nearest_coprime(n_points, period)}) is strictly "
        f"better at the same cost.",
    )


def _nearest_coprime(n: int, period: int) -> int:
    """The closest sample count to ``n`` that is coprime to ``period``."""
    for d in range(1, period + 1):
        for cand in (n + d, n - d):
            if cand > 0 and gcd(cand, period) == 1:
                return cand
    return n + 1
