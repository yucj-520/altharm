"""Wannier-route tests.

Independence rule: the reference is a closed-form band structure written out by
hand, or a hand-built file whose contents the test itself knows.  Nothing is
checked against a Hamiltonian this code produced.
"""

import numpy as np
import pytest

from altharm.wannier import (
    read_hr, hk, eigenvalues, cart_to_frac_k, truncate, hopping_range, spin_splitting,
)


def write_hr(path, num_wann, entries, deg=None):
    """Write a hr.dat.  ``entries`` maps R -> (nw,nw) complex matrix H[i,j].

    Wannier90's column order is ``Rx Ry Rz i j Re Im`` with ``i`` running
    fastest, and the test writes it that way explicitly so that the reader's
    index convention is pinned by this file rather than by the reader itself.
    """
    rs = list(entries)
    deg = deg or [1] * len(rs)
    lines = [" written by the test suite", f"{num_wann}", f"{len(rs)}"]
    for start in range(0, len(rs), 15):
        lines.append(" ".join(str(d) for d in deg[start:start + 15]))
    for r, H in entries.items():
        for j in range(num_wann):
            for i in range(num_wann):
                lines.append(
                    f"{r[0]:5d}{r[1]:5d}{r[2]:5d}{i+1:5d}{j+1:5d}"
                    f"{H[i, j].real:15.8f}{H[i, j].imag:15.8f}"
                )
    path.write_text("\n".join(lines) + "\n")
    return path


def test_one_band_chain_against_its_closed_form(tmp_path):
    """eps(k) = e0 + 2 t cos(2 pi k), written out by hand."""
    e0, t = 0.2, -0.5
    entries = {
        (-1, 0, 0): np.array([[t]], dtype=complex),
        (0, 0, 0): np.array([[e0]], dtype=complex),
        (1, 0, 0): np.array([[t]], dtype=complex),
    }
    wh = read_hr(write_hr(tmp_path / "c_hr.dat", 1, entries))
    assert wh.num_wann == 1 and wh.num_rpts == 3

    kf = np.stack([np.linspace(0, 0.5, 11), np.zeros(11), np.zeros(11)], axis=1)
    got = eigenvalues(wh, kf)[:, 0]
    ref = e0 + 2 * t * np.cos(2 * np.pi * kf[:, 0])
    assert np.allclose(got, ref, atol=1e-12), f"max dev {np.abs(got - ref).max():.2e}"


def test_index_convention_is_not_transposed(tmp_path):
    """An asymmetric on-site block: H[0,1] and H[1,0] must not be swapped."""
    H = np.array([[0.0, 1.0 + 2.0j], [1.0 - 2.0j, 0.0]], dtype=complex)
    wh = read_hr(write_hr(tmp_path / "c_hr.dat", 2, {(0, 0, 0): H}))
    assert np.allclose(wh.ham[0], H, atol=1e-12), f"read back\n{wh.ham[0]}\nwanted\n{H}"
    assert np.allclose(hk(wh, [[0.0, 0.0, 0.0]])[0], H, atol=1e-12)


def test_degeneracy_divides_the_hopping(tmp_path):
    """A degeneracy of 2 must halve that R's contribution."""
    entries = {(0, 0, 0): np.array([[1.0]], dtype=complex)}
    plain = read_hr(write_hr(tmp_path / "a_hr.dat", 1, entries, deg=[1]))
    halved = read_hr(write_hr(tmp_path / "b_hr.dat", 1, entries, deg=[2]))
    k = [[0.1, 0.2, 0.3]]
    assert np.isclose(hk(plain, k)[0, 0, 0].real, 1.0)
    assert np.isclose(hk(halved, k)[0, 0, 0].real, 0.5)


def test_hk_is_hermitian_for_a_hermitian_model(tmp_path):
    rng = np.random.default_rng(5)
    A = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    entries = {(0, 0, 0): A + A.conj().T, (1, 0, 0): A, (-1, 0, 0): A.conj().T}
    wh = read_hr(write_hr(tmp_path / "c_hr.dat", 3, entries))
    H = hk(wh, [[0.13, 0.0, 0.0]])[0]
    assert np.allclose(H, H.conj().T, atol=1e-12)


def test_cart_to_frac_k_against_a_hand_built_reciprocal_lattice():
    """Hexagonal cell; the reciprocal vectors are formed independently here."""
    a, c = 3.7258396285, 16.9179784926
    cell = np.array([[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]])
    a1, a2, a3 = cell
    vol = float(a1 @ np.cross(a2, a3))
    b1 = 2 * np.pi * np.cross(a2, a3) / vol
    b2 = 2 * np.pi * np.cross(a3, a1) / vol
    b3 = 2 * np.pi * np.cross(a1, a2) / vol

    frac = np.array([[0.3, -0.2, 0.15]])
    kcart = frac[0, 0] * b1 + frac[0, 1] * b2 + frac[0, 2] * b3
    assert np.allclose(cart_to_frac_k(kcart, cell)[0], frac[0], atol=1e-12)


def test_truncate_drops_the_far_hoppings(tmp_path):
    cell = np.eye(3) * 2.0
    entries = {(n, 0, 0): np.array([[0.1 / (abs(n) + 1)]], dtype=complex)
               for n in range(-3, 4)}
    wh = read_hr(write_hr(tmp_path / "c_hr.dat", 1, entries), cell=cell)
    assert np.isclose(hopping_range(wh), 6.0)
    cut = truncate(wh, 4.1)
    assert cut.num_rpts == 5                       # R = -2..2, i.e. |R| <= 4 A
    assert np.isclose(hopping_range(cut), 4.0)

    # A radius below the first neighbour keeps the on-site term and nothing
    # else.  |R| = 0 always passes, so this is not an error case.
    onsite = truncate(wh, 0.5)
    assert onsite.num_rpts == 1
    assert tuple(onsite.rvecs[0]) == (0, 0, 0)


def test_truncate_refuses_when_nothing_survives():
    """With no on-site term there is a radius that keeps nothing, and that is an error."""
    from altharm.wannier import WannierHam
    wh = WannierHam(
        rvecs=np.array([[1, 0, 0], [-1, 0, 0]]),
        deg=np.ones(2, int),
        ham=np.full((2, 1, 1), 0.1, dtype=complex),
        cell=np.eye(3) * 2.0,
    )
    with pytest.raises(ValueError, match="keeps no hoppings"):
        truncate(wh, 0.5)


def test_truncation_damps_the_high_harmonic_first():
    """The claim the package rests on, on a model whose harmonics are known.

    Build a one-band 2D model whose dispersion is exactly cos(6 theta) at fixed
    |k| by putting equal hoppings on a ring of six lattice vectors, plus a
    nearer ring carrying a cos(2 theta) part.  Cutting the range must remove the
    l = 6 part and leave the l = 2 part, which is what "a truncated Wannier
    model cannot carry high harmonics" means quantitatively.
    """
    from altharm.basis import CircularBasis
    from altharm.kpath import circle
    from altharm.project import project
    from altharm.wannier import WannierHam

    cell = np.eye(3) * 1.0
    near = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]
    far = [(3, 0, 0), (-3, 0, 0), (0, 3, 0), (0, -3, 0)]
    rvecs, ham = [], []
    for r in near:
        rvecs.append(r)
        ham.append([[0.10 if r[0] else -0.10]])
    for r in far:
        rvecs.append(r)
        ham.append([[0.02 if r[0] else -0.02]])
    wh = WannierHam(
        rvecs=np.array(rvecs), deg=np.ones(len(rvecs), int),
        ham=np.array(ham, dtype=complex), cell=cell,
    )

    basis = CircularBasis(8)
    pts = circle(0.9, 64)
    frac = cart_to_frac_k(pts, cell)

    full = project(basis, pts, eigenvalues(wh, frac)[:, 0]).amplitude_by_l()
    cut = project(basis, pts, eigenvalues(truncate(wh, 1.5), frac)[:, 0]).amplitude_by_l()

    assert full[6] > 1e-4, "the full model must carry an l=6 amplitude to begin with"
    assert cut[6] < 0.02 * full[6], (
        f"truncation left {cut[6]:.3g} of {full[6]:.3g} at l=6"
    )
    assert cut[2] > 0.5 * full[2], "the low harmonic must survive the same cut"


def test_spin_splitting_refuses_mismatched_models(tmp_path):
    a = read_hr(write_hr(tmp_path / "a_hr.dat", 1, {(0, 0, 0): np.ones((1, 1), complex)}))
    b = read_hr(write_hr(tmp_path / "b_hr.dat", 2, {(0, 0, 0): np.eye(2, dtype=complex)}))
    with pytest.raises(ValueError, match="different sizes"):
        spin_splitting(a, b, [[0.0, 0.0, 0.0]], 0)


def test_read_hr_refuses_a_truncated_file(tmp_path):
    p = write_hr(tmp_path / "c_hr.dat", 2, {(0, 0, 0): np.eye(2, dtype=complex)})
    lines = p.read_text().splitlines()
    p.write_text("\n".join(lines[:-1]) + "\n")
    with pytest.raises(ValueError, match="expected .* numbers"):
        read_hr(p)
