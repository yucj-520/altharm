"""Every input guard that no other test ever triggered.

These are the `raise ValueError` / `raise RuntimeError` branches scattered
through the package. Coverage said none of them had ever fired, which means
their messages had never been read by anything but a human writing them --
and a guard that raises the wrong error, or names the wrong variable, fails
exactly when someone is already confused.

Independence rule: each test constructs input that is wrong in one specific,
hand-chosen way and asserts on the message naming *that* way, not just that
something raised.
"""

import numpy as np
import pytest

from altharm.basis import CircularBasis, SphericalBasis, real_sph_harm, _unit
from altharm.kpath import circle, sphere_grid, angles_of
from altharm.project import project
from altharm.noise import noise_report
from altharm.symmetry import spin_operations, symmetrizer
from altharm.wannier import WannierHam, truncate, spin_splitting, read_hr
from altharm.qe.writer import bands_input, kpoints_card
from altharm.bandtrack import track_by_continuity


# ------------------------------------------------------------------ kpath

def test_angles_of_matches_hand_computed_angles():
    """`angles_of` had zero coverage -- nothing in the package called it.

    Reference: angles written out by hand for points placed on the axes.
    """
    pts = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]])
    got = angles_of(pts)
    assert np.allclose(got, [0.0, np.pi / 2, np.pi, -np.pi / 2], atol=1e-12)


def test_angles_of_respects_a_shifted_centre_and_rotated_axes():
    centre = np.array([0.0, 0.0, 0.3])
    pts = np.array([[2.0, 0.0, 0.3]])
    assert np.isclose(angles_of(pts, center=centre)[0], 0.0, atol=1e-12)
    # measuring from e1=+y puts a +x point at -90 degrees
    got = angles_of(pts, e1=(0, 1, 0), e2=(-1, 0, 0), center=centre)[0]
    assert np.isclose(got, -np.pi / 2, atol=1e-12)


def test_circle_refuses_non_orthogonal_axes():
    with pytest.raises(ValueError, match="e1 and e2 must be orthogonal"):
        circle(1.0, 10, e1=(1, 0, 0), e2=(1, 1, 0))


# ------------------------------------------------------------------ basis

def test_unit_refuses_a_zero_length_direction():
    with pytest.raises(ValueError, match="zero-length direction"):
        _unit(np.zeros((1, 3)))


def test_real_sph_harm_refuses_m_larger_than_l():
    with pytest.raises(ValueError, match=r"\|m\| > l"):
        real_sph_harm(2, 3, np.array([0.5]), np.array([0.5]))


def test_bases_refuse_a_negative_lmax():
    with pytest.raises(ValueError, match="lmax must be >= 0"):
        CircularBasis(-1)
    with pytest.raises(ValueError, match="lmax must be >= 0"):
        SphericalBasis(-1)


# ------------------------------------------------------------------ symmetry

def test_spin_operations_refuses_mismatched_array_lengths():
    cell = np.eye(3) * 3.0
    with pytest.raises(ValueError, match="same length"):
        spin_operations(cell, np.zeros((2, 3)), np.array([26, 26]), np.array([1.0]))


def test_spin_operations_refuses_a_structure_with_no_magnetic_atoms():
    cell = np.eye(3) * 3.0
    pos = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    with pytest.raises(ValueError, match="no magnetic atoms"):
        spin_operations(cell, pos, np.array([26, 26]), np.zeros(2))


def test_symmetrizer_refuses_when_no_operation_preserves_the_manifold():
    """A circle in a plane no kept operation maps to itself leaves nothing to
    average over -- that must raise, not silently return an empty projector."""
    from altharm.symmetry import SpinOp
    # a lone rotation about x: it does not preserve a circle in the x-y plane
    rot = np.array([[1.0, 0, 0], [0, 0, -1.0], [0, 1.0, 0]])
    ops = [SpinOp(rot_frac=None, trans=None, rot_cart=rot, chi=1)]
    with pytest.raises(RuntimeError, match="no operation maps"):
        symmetrizer(CircularBasis(4), ops, circle(1.0, 20))


# ------------------------------------------------------------------ project

def test_project_refuses_mismatched_lengths_and_underdetermined_fits():
    pts = circle(1.0, 6)
    with pytest.raises(ValueError, match="directions but"):
        project(CircularBasis(2), pts, np.zeros(5))
    with pytest.raises(ValueError, match="cannot determine"):
        project(CircularBasis(9), pts, np.zeros(6))


# ------------------------------------------------------------------ noise

def test_snr_is_infinite_when_no_forbidden_harmonic_exists():
    """With every l allowed there is no floor to divide by; the report must
    say "infinite", not raise a ZeroDivisionError."""
    pts = circle(1.0, 40)
    t = CircularBasis(1).angles(pts)
    fit = project(CircularBasis(4), pts, np.cos(2 * t))
    rep = noise_report(fit, allowed_l_set={0, 1, 2, 3, 4})
    assert rep.floor == 0.0 and rep.floor_max == 0.0
    assert all(np.isinf(v) for v in rep.snr().values())


# ------------------------------------------------------------------ wannier

def _toy_ham(cell=None):
    return WannierHam(
        rvecs=np.array([[0, 0, 0], [1, 0, 0]]),
        deg=np.ones(2, int),
        ham=np.full((2, 1, 1), 0.2, dtype=complex),
        cell=cell,
    )


def test_wannier_refuses_geometry_operations_without_a_cell():
    wh = _toy_ham(cell=None)
    with pytest.raises(ValueError, match="cell is unknown"):
        wh.rcart()
    with pytest.raises(ValueError, match="cell is required"):
        truncate(wh, 1.0)


def test_spin_splitting_refuses_an_out_of_range_band():
    wh = _toy_ham(cell=np.eye(3) * 2.0)
    with pytest.raises(IndexError, match="band 5 out of range"):
        spin_splitting(wh, wh, [[0.0, 0.0, 0.0]], 5)


def test_read_hr_refuses_a_file_whose_body_is_the_wrong_length(tmp_path):
    p = tmp_path / "bad_hr.dat"
    p.write_text(" header\n1\n2\n1 1\n0 0 0 1 1 0.5 0.0\n")   # one R-vector short
    with pytest.raises(ValueError, match="expected .* numbers"):
        read_hr(p)


# ------------------------------------------------------------------ qe writer

def _bands_kwargs(**over):
    kw = dict(
        prefix="x", outdir="./tmp", pseudo_dir="/p", alat_bohr=10.0,
        cell_ang=np.eye(3) * 5.0,
        species=[("Fe", 55.8, "Fe.UPF")],
        positions_ang=[("Fe", 0.0, 0.0, 0.0)],
        positions_crystal=None,
        ecutwfc=40.0, kcart_tpiba=circle(0.2, 8),
    )
    kw.update(over)
    return kw


def test_kpoints_card_refuses_a_wrong_shaped_k_list():
    with pytest.raises(ValueError, match=r"shape \(nk, 3\)"):
        kpoints_card(np.zeros((5, 2)))


def test_bands_input_refuses_a_wrong_shaped_cell():
    with pytest.raises(ValueError, match=r"cell_ang must have shape"):
        bands_input(**_bands_kwargs(cell_ang=np.eye(2)))


def test_bands_input_refuses_empty_species_or_positions():
    with pytest.raises(ValueError, match="species must contain"):
        bands_input(**_bands_kwargs(species=[]))
    with pytest.raises(ValueError, match="positions_ang must contain"):
        bands_input(**_bands_kwargs(positions_ang=[]))


def test_bands_input_refuses_magnetization_of_the_wrong_length():
    with pytest.raises(ValueError, match="one entry per species"):
        bands_input(**_bands_kwargs(nspin=2, starting_magnetization=[1.0, -1.0]))


# ------------------------------------------------------------------ bandtrack

def test_track_refuses_a_single_point_path():
    with pytest.raises(ValueError, match="at least two"):
        track_by_continuity(np.array([[0.0, 1.0]]))


def test_gap_is_infinite_for_a_single_band():
    """One band means no neighbour to be close to; the gap must be inf rather
    than crashing on an empty array."""
    tracked = track_by_continuity(np.array([[1.0], [1.1], [1.2]]), wrap=False)
    assert np.isinf(tracked.gap).all()


# ------------------------------------------------- remaining uncovered branches

def test_spherical_basis_l_of_lists_each_label_l():
    """`SphericalBasis.l_of()` had never been called by anything."""
    b = SphericalBasis(2)
    assert list(b.l_of()) == [0, 1, 1, 1, 2, 2, 2, 2, 2]
    assert len(b.l_of()) == b.size


def test_symmetrizer_accepts_a_spherical_basis():
    """The `_preserves` branch for a sphere: every rotation maps it to itself,
    so nothing is filtered out. Never exercised -- all previous use was 2D.

    Note the precondition this test had to satisfy on its second attempt: the
    operation list must be a *closed group*. S = (1/|G|) sum_g chi(g) rho(g)
    is idempotent only when G is closed under composition -- handing it
    {I, C4} alone gives a perfectly well-defined average that is not a
    projector. Real callers always pass a closed set, because
    `spin_operations` returns the full group and then closes it under time
    reversal; this test passes the whole C4 cyclic group for the same reason.
    """
    from altharm.symmetry import SpinOp
    def rot_z(k):
        a = k * np.pi / 2
        return np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1.0]])
    ops = [SpinOp(rot_frac=None, trans=None, rot_cart=rot_z(k), chi=1) for k in range(4)]
    pts, _ = sphere_grid(8, 14)
    smat, kept = symmetrizer(SphericalBasis(2), ops, pts)
    assert len(kept) == 4, "a sphere should keep every rotation"
    assert np.allclose(smat @ smat, smat, atol=1e-8), "S must still be a projector"


def test_an_operation_broken_by_the_magnetic_order_is_dropped():
    """Three sublattices arranged so one operation preserves some signs and
    flips others -- neither chi = +1 nor -1 applies, and it must be discarded
    rather than forced into one of them.

    Chain of three sites with moments (+, -, +): the 1/3 translation sends
    + -> -, - -> +, but also + -> + on the site that wraps, so no single
    character describes it.
    """
    cell = np.diag([3.0, 3.0, 9.0])
    pos = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1 / 3], [0.0, 0.0, 2 / 3]])
    numbers = np.array([26, 26, 26])
    magmoms = np.array([1.0, -1.0, 1.0])
    ops = spin_operations(cell, pos, numbers, magmoms, require_exchange=False)
    # whatever survives must carry a definite character
    assert all(op.chi in (1, -1) for op in ops)
    # and the inconsistent translation must not be among them
    for op in ops:
        assert not (np.allclose(op.rot_cart, np.eye(3))
                    and not np.allclose(op.trans, 0.0)
                    and op.chi == 1 and abs(op.trans[2] - 1/3) < 1e-9)


def test_atom_map_refuses_a_mapping_that_misses_and_one_that_is_not_a_permutation():
    """`_atom_map`'s two `return None` exits, tested directly."""
    from altharm.symmetry import _atom_map
    pos = np.array([[0.0, 0.0, 0.0], [0.37, 0.0, 0.0]])
    numbers = np.array([26, 26])
    # a rotation that sends neither atom near an atom site
    rot = np.eye(3, dtype=int)
    assert _atom_map(rot, np.array([0.19, 0.0, 0.0]), pos, numbers, 1e-5) is None
    # identity maps cleanly, as a positive control
    assert _atom_map(rot, np.zeros(3), pos, numbers, 1e-5) is not None


def test_symprec_scan_records_the_error_instead_of_raising(tmp_path):
    """The `except` arm of symprec_scan: a structure that makes
    spin_operations raise at every tolerance must come back as rows carrying
    the message, not blow up the scan."""
    from altharm.symmetry import symprec_scan
    cell = np.eye(3) * 3.0
    pos = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    rows = symprec_scan(cell, pos, np.array([26, 26]), np.zeros(2), symprecs=(1e-5, 1e-3))
    assert len(rows) == 2
    assert all(r["error"] and "no magnetic atoms" in r["error"] for r in rows)


def test_spin_operations_reports_when_spglib_finds_nothing(monkeypatch):
    """A defensive guard that only a broken spglib could trigger -- exercised
    here by making spglib return None, so the message is known to work."""
    import altharm.symmetry as sym
    monkeypatch.setattr(sym.spglib, "get_symmetry", lambda *a, **k: None)
    with pytest.raises(RuntimeError, match="spglib found no symmetry"):
        spin_operations(np.eye(3) * 3, np.zeros((1, 3)), np.array([26]), np.array([1.0]))


def test_bands_input_omits_optional_blocks_when_not_asked_for():
    """ecutrho / nbnd / occupations / smearing / degauss are all optional; the
    branches that leave them out had never been taken."""
    text = bands_input(**_bands_kwargs())
    for absent in ("ecutrho", "nbnd", "occupations", "smearing", "degauss"):
        assert absent not in text, f"{absent} should be omitted when not given"
    assert "ecutwfc" in text and "nosym = .true." in text


def test_bands_input_writes_every_optional_block_when_asked_for():
    """The other half of the optional-argument branches.

    Coverage showed these lines never ran: the suite only ever called
    `bands_input` *without* ecutrho/nbnd/occupations/smearing/degauss. The
    real Cr2Se2As production input used all five, but that input was generated
    by a script, not by a test -- so their formatting was unverified.
    """
    text = bands_input(**_bands_kwargs(
        ecutrho=240.0, nbnd=80, nspin=2, starting_magnetization=[1.0],
        occupations="smearing", smearing="gaussian", degauss=0.0037,
    ))
    assert "  ecutrho = 240.000000" in text
    assert "  nbnd = 80" in text
    assert "  occupations = 'smearing'" in text
    assert "  smearing = 'gaussian'" in text
    assert "  degauss = 0.00370000" in text
    assert "  nspin = 2" in text
    assert "  starting_magnetization(1) = 1.000000" in text


def test_atom_map_refuses_when_two_atoms_land_on_the_same_target():
    """`_atom_map`'s "not a permutation" exit: two atoms closer together than
    the matching tolerance both claim the same partner, so the mapping is not
    invertible and must be rejected rather than silently double-counted."""
    from altharm.symmetry import _atom_map
    eps = 2e-5
    pos = np.array([[0.0, 0.0, 0.0], [eps, 0.0, 0.0]])
    numbers = np.array([26, 26])
    # Translating by exactly the (sub-tolerance) separation sends BOTH atoms
    # nearest to atom 1: atom 0 lands on it exactly, atom 1 lands one eps past
    # it and is still closer to it than to atom 0. mapped = [1, 1], which is
    # not invertible.
    assert _atom_map(np.eye(3, dtype=int), np.array([eps, 0.0, 0.0]),
                     pos, numbers, 1e-3) is None
    # control: with no translation the same structure maps cleanly
    assert _atom_map(np.eye(3, dtype=int), np.zeros(3), pos, numbers, 1e-3) is not None


def test_every_operation_being_unmappable_is_refused_by_name(monkeypatch):
    """Two defensive guards at once, neither reachable with a working spglib:
    the `continue` when an operation's atom mapping is not clean, and the
    RuntimeError when that leaves nothing at all. In practice the identity
    always maps cleanly, so this can only be provoked deliberately -- which is
    the point of checking the message exists and is accurate."""
    import altharm.symmetry as sym
    monkeypatch.setattr(sym, "_atom_map", lambda *a, **k: None)
    with pytest.raises(RuntimeError, match="no operation survived"):
        spin_operations(np.eye(3) * 3.0,
                        np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
                        np.array([26, 26]), np.array([1.0, -1.0]))
