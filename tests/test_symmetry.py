"""Symmetry tests.

Independence rule: the reference is textbook physics about a named structure,
worked out by hand in the docstring, not anything this code produced.  The three
structures below were chosen because their answers are unambiguous:

* rutile with the two metal sites antialigned is the canonical d-wave
  altermagnet, and its splitting goes as k_x k_y,
* two sublattices related by a pure lattice translation is a conventional
  antiferromagnet and has no splitting at all,
* likewise two sublattices related by inversion.

The last two are the definition that separates an altermagnet from an ordinary
antiferromagnet, so a code that got them wrong would be worthless.
"""

import numpy as np
import pytest

from altharm.basis import CircularBasis
from altharm.kpath import circle
from altharm.symmetry import spin_operations, symmetrizer, allowed_harmonics, allowed_l


def rutile(u=0.3058, a=4.49, c=3.11):
    """RuO2-type. Ru at the corner and the body centre, antialigned."""
    cell = np.diag([a, a, c])
    pos = np.array([
        [0.0, 0.0, 0.0],            # Ru up
        [0.5, 0.5, 0.5],            # Ru down
        [u, u, 0.0],                # O
        [1 - u, 1 - u, 0.0],
        [0.5 + u, 0.5 - u, 0.5],
        [0.5 - u, 0.5 + u, 0.5],
    ])
    numbers = np.array([44, 44, 8, 8, 8, 8])
    magmoms = np.array([1.0, -1.0, 0.0, 0.0, 0.0, 0.0])
    return cell, pos, numbers, magmoms


def translation_afm():
    """Sublattices related by the body-centring translation: conventional AFM."""
    cell = np.diag([3.0, 3.0, 3.0])
    pos = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    return cell, pos, np.array([26, 26]), np.array([1.0, -1.0])


def inversion_afm():
    """Sublattices related by inversion through the origin: conventional AFM."""
    cell = np.diag([4.0, 4.0, 4.0])
    pos = np.array([[0.2, 0.0, 0.0], [0.8, 0.0, 0.0]])
    return cell, pos, np.array([26, 26]), np.array([1.0, -1.0])


def hexagonal():
    a, c = 4.0, 6.0
    cell = np.array([[a, 0.0, 0.0], [-a / 2, a * np.sqrt(3) / 2, 0.0], [0.0, 0.0, c]])
    pos = np.array([[1 / 3, 2 / 3, 0.25], [2 / 3, 1 / 3, 0.75]])
    return cell, pos, np.array([25, 25]), np.array([1.0, -1.0])


def xy_setup(lmax=8, n=None):
    basis = CircularBasis(lmax)
    n = n or (4 * lmax + 8)
    return basis, circle(1.0, n)


# ---------------------------------------------------------------- conventions

def test_cartesian_rotations_are_orthogonal_in_a_hexagonal_cell():
    """A wrong fractional->Cartesian convention survives a cubic cell and dies here.

    The lattice is non-orthogonal, so A^T R A^-T is orthogonal only if the
    convention matches spglib's.  This is the cheapest guard against the metric
    bug the package exists to avoid.
    """
    ops = spin_operations(*hexagonal())
    assert len(ops) > 1
    for op in ops:
        R = op.rot_cart
        assert np.allclose(R @ R.T, np.eye(3), atol=1e-10), (
            f"rot_cart not orthogonal:\n{R}"
        )
        assert np.isclose(abs(np.linalg.det(R)), 1.0, atol=1e-10)


def test_symmetrizer_is_a_projector():
    """S**2 = S.  Independent of any physics: it follows from S being a group average.

    It fails if the characters are inconsistent, if the set is not closed, or if
    the Cartesian rotations are wrong -- so it catches most ways of being wrong.
    """
    basis, dirs = xy_setup()
    ops = spin_operations(*rutile())
    S, kept = symmetrizer(basis, ops, dirs)
    assert np.allclose(S @ S, S, atol=1e-8), (
        f"max |S^2 - S| = {np.abs(S @ S - S).max():.2e} over {len(kept)} operations"
    )
    w = np.linalg.eigvals(S)
    assert np.allclose(np.sort_complex(w).real, np.round(np.sort_complex(w).real), atol=1e-8)


# --------------------------------------------------------------- the physics

def test_rutile_is_d_wave():
    """Rutile: l = 2 allowed, l = 4 forbidden, and the allowed function is sin(2t).

    Worked out by hand.  The fourfold rotation about z exchanges the two Ru
    sites, so chi(C4z) = -1, and cos/sin(4t) are invariant under a 90 degree
    turn -- they are averaged against -1 and die.  cos/sin(2t) change sign under
    the same turn, so they survive.  The {110} mirror is sublattice-preserving,
    chi = +1, and sends t -> pi/2 - t, which kills cos(2t) and keeps sin(2t).
    What is left is sin(2t), i.e. the splitting goes as k_x k_y: d-wave.
    """
    basis, dirs = xy_setup(lmax=8)
    ops = spin_operations(*rutile())
    allowed = allowed_harmonics(basis, ops, dirs)

    assert set(allowed) == {2, 6}, f"allowed l = {sorted(allowed)}, expected {{2, 6}}"
    assert 4 not in allowed, "g-wave harmonic must be forbidden in a d-wave altermagnet"

    # inside the l = 2 block the basis order is (cos, sin)
    vec = allowed[2]
    assert vec.shape[0] == 1, f"l=2 subspace has dimension {vec.shape[0]}, expected 1"
    cos_weight, sin_weight = abs(vec[0, 0]), abs(vec[0, 1])
    assert cos_weight < 1e-8 and abs(sin_weight - 1.0) < 1e-8, (
        f"l=2 allowed function is not sin(2t): weights cos={cos_weight:.3g} sin={sin_weight:.3g}"
    )


def test_odd_harmonics_are_forbidden_without_spin_orbit():
    """Without spin-orbit each spin channel is time-reversal symmetric, so
    dE(-k) = dE(k) and every odd harmonic must vanish."""
    basis, dirs = xy_setup(lmax=9)
    ops = spin_operations(*rutile())
    got = allowed_l(basis, ops, dirs)
    assert not any(l % 2 for l in got), f"odd harmonics allowed: {sorted(got)}"


def test_translation_related_antiferromagnet_has_no_splitting():
    """The textbook case: sublattices related by a lattice translation.

    The translation acts as the identity on k but carries chi = -1, so it is
    averaged against the identity's +1 and annihilates everything.  A finite
    splitting here would mean the code cannot tell an altermagnet from an
    ordinary antiferromagnet.
    """
    basis, dirs = xy_setup(lmax=8)
    ops = spin_operations(*translation_afm())
    assert allowed_l(basis, ops, dirs) == set(), "conventional AFM must allow nothing"


def test_inversion_related_antiferromagnet_has_no_splitting():
    basis, dirs = xy_setup(lmax=8)
    ops = spin_operations(*inversion_afm())
    assert allowed_l(basis, ops, dirs) == set(), "conventional AFM must allow nothing"


def test_characters_are_consistent():
    """Every returned operation is tagged +1 or -1, and both occur for rutile."""
    ops = spin_operations(*rutile())
    chis = {op.chi for op in ops}
    assert chis == {1, -1}
    assert all(op.chi in (1, -1) for op in ops)


def test_no_magnetic_atoms_is_refused():
    cell, pos, numbers, _ = rutile()
    with pytest.raises(ValueError, match="no magnetic atoms"):
        spin_operations(cell, pos, numbers, np.zeros(len(numbers)))


# ------------------------------------------------- the symprec cliff

def almost_rutile(jitter=3e-5):
    """Rutile with the mirror broken by a jitter smaller than a loose symprec.

    A relaxed structure is only approximately symmetric.  This reproduces that
    situation on purpose, because the operation that gets lost first is the one
    that exchanges the sublattices -- and losing it turns an altermagnet into
    something the code would describe as a ferromagnet.
    """
    cell, pos, numbers, magmoms = rutile()
    pos = pos.copy()
    pos[1, 2] += jitter          # nudge the body-centre Ru off its site
    return cell, pos, numbers, magmoms


def test_losing_the_exchanging_operation_is_refused_not_ignored():
    """With no chi = -1 operation the code must refuse, naming the cause.

    Silence here is the worst possible failure: every harmonic including l = 0
    would come out allowed, so a compensated magnet would be reported as if it
    carried a net splitting.
    """
    from altharm.symmetry import SymmetryToleranceError

    with pytest.raises(SymmetryToleranceError, match="exchanges the"):
        spin_operations(*almost_rutile(), symprec=1e-6)


def test_the_escape_hatch_is_explicit():
    """require_exchange=False returns the degraded set rather than refusing."""
    ops = spin_operations(*almost_rutile(), symprec=1e-6, require_exchange=False)
    assert all(op.chi == +1 for op in ops)


def test_a_looser_tolerance_recovers_the_operation():
    ops = spin_operations(*almost_rutile(), symprec=1e-3)
    assert any(op.chi == -1 for op in ops)


def test_symprec_scan_shows_where_the_group_changes():
    """The scan must bracket the jitter: tight tolerances lose it, loose ones keep it."""
    from altharm.symmetry import symprec_scan, format_symprec_scan

    rows = symprec_scan(*almost_rutile(jitter=3e-5))
    tight = [r for r in rows if r["symprec"] <= 1e-5]
    loose = [r for r in rows if r["symprec"] >= 1e-3]
    assert all(r["n_exchanging"] == 0 for r in tight), "a tight symprec should lose it"
    assert all(r["n_exchanging"] > 0 for r in loose), "a loose symprec should keep it"
    assert "symprec" in format_symprec_scan(rows)


def test_a_ferromagnet_is_not_refused():
    """The guard fires only when the moments cancel; a ferromagnet is legitimate."""
    cell, pos, numbers, _ = rutile()
    magmoms = np.array([1.0, 1.0, 0.0, 0.0, 0.0, 0.0])
    ops = spin_operations(cell, pos, numbers, magmoms)
    assert all(op.chi == +1 for op in ops)
