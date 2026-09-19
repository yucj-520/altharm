"""Which angular harmonics the magnetic symmetry allows.

The selection rule is not looked up in a table here.  It is derived numerically
from the structure and its magnetic configuration, so that it holds for whatever
cell the user actually hands in -- including a strained or distorted one, where
the tabulated answer for the ideal lattice no longer applies.

The argument runs as follows.  Take a collinear magnet with two magnetic
sublattices carrying opposite moments, and ignore spin-orbit coupling, so the
two spin channels are separate spinless problems.  Let ``g`` be an operation of
the *non-magnetic* parent crystal.  Two things can happen:

* ``g`` maps each magnetic sublattice to itself.  Then it is a symmetry of the
  magnetic structure outright, so :math:`E_\\uparrow(Rk) = E_\\uparrow(k)`, and the
  splitting obeys :math:`\\Delta E(Rk) = +\\Delta E(k)`.
* ``g`` exchanges the two sublattices.  Then ``g`` maps the structure onto the
  one with every moment reversed, so ``g`` combined with a spin flip is a
  symmetry, :math:`E_\\uparrow(Rk) = E_\\downarrow(k)`, and
  :math:`\\Delta E(Rk) = -\\Delta E(k)`.

That sign is the character :math:`\\chi(g) = \\pm 1`.  A harmonic survives only
if it is invariant under the operator

.. math:: (S f)(\\hat k) = \\frac{1}{|G|} \\sum_g \\chi(g)\\, f(R_g^{-1}\\hat k)

which is a projector.  Its eigenvectors of eigenvalue 1 span the allowed
subspace; everything orthogonal to them is forbidden, and is therefore a direct
handle on the numerical noise (see :mod:`altharm.noise`).

The whole point of building the projector numerically is that we never have to
be right about a character table.  ``S**2 == S`` is checked in the tests, and a
wrong Cartesian convention would break it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import spglib


@dataclass(frozen=True)
class SpinOp:
    """One operation of the non-magnetic parent, tagged by its effect on the splitting.

    Attributes
    ----------
    rot_frac : (3,3) int array
        Rotation in the basis of the lattice vectors, as spglib returns it.
    trans : (3,) array
        Fractional translation.
    rot_cart : (3,3) array
        The same rotation in Cartesian axes.  Orthogonal by construction; the
        tests assert it.
    chi : int
        ``+1`` if the operation preserves the magnetic sublattices, ``-1`` if it
        exchanges them.
    label : str
        ``'identity'``, ``'time-reversal'`` or ``'spglib'``, recording where the
        operation came from.
    """

    rot_frac: np.ndarray
    trans: np.ndarray
    rot_cart: np.ndarray
    chi: int
    label: str = "spglib"


def _cartesian_rotation(rot_frac, cell):
    """Fractional rotation -> Cartesian, for a cell whose ROWS are lattice vectors.

    With row-vector positions ``r = x A`` and spglib's convention
    ``x'_i = R_ij x_j``, the Cartesian operation on column vectors is
    ``A^T R A^{-T}``.  It comes out orthogonal whenever ``R`` is a genuine
    operation of that lattice, which is the check the tests lean on.
    """
    A = np.asarray(cell, dtype=float)
    return A.T @ np.asarray(rot_frac, dtype=float) @ np.linalg.inv(A.T)


def spin_operations(
    cell,
    scaled_positions,
    numbers,
    magmoms,
    symprec: float = 1e-5,
    include_time_reversal: bool = True,
    moment_tol: float = 1e-3,
    require_exchange: bool = True,
):
    """Operations of the non-magnetic parent, each tagged with its character.

    Parameters
    ----------
    cell : (3,3) array
        Lattice vectors as ROWS, in Angstrom (ASE's convention).
    scaled_positions : (n,3) array
        Fractional coordinates.
    numbers : (n,) int array
        Atomic numbers.
    magmoms : (n,) array
        Collinear moments.  Only their signs matter here.
    include_time_reversal : bool
        Add the inversion :math:`k \\to -k` with character ``+1``.  Without
        spin-orbit coupling each spin channel is a time-reversal-symmetric
        spinless problem, so :math:`E_\\sigma(k) = E_\\sigma(-k)` and the
        splitting is even in ``k`` whether or not the crystal has inversion.
        This is what forbids all odd harmonics.  **Set it False if the
        calculation included spin-orbit coupling**, where it does not hold.
    require_exchange : bool
        Refuse to return a set that has no sublattice-exchanging operation when
        the moments cancel.  Such a set cannot distinguish an altermagnet from a
        ferromagnet, and the usual cause is a ``symprec`` slightly too tight for
        a relaxed structure.  See :func:`symprec_scan`.

    Returns
    -------
    list[SpinOp]
        Operations that are symmetries of the magnetic structure either
        directly or combined with a spin flip.  Operations that are neither are
        dropped: they are broken by the magnetic order and constrain nothing.
    """
    cell = np.asarray(cell, dtype=float)
    scaled_positions = np.asarray(scaled_positions, dtype=float)
    numbers = np.asarray(numbers, dtype=int)
    magmoms = np.asarray(magmoms, dtype=float)

    if not (len(scaled_positions) == len(numbers) == len(magmoms)):
        raise ValueError("positions, numbers and magmoms must have the same length")

    # The NON-magnetic parent: chemical species only, magnetic order ignored.
    sym = spglib.get_symmetry((cell, scaled_positions, numbers), symprec=symprec)
    if sym is None:
        raise RuntimeError("spglib found no symmetry for this cell")

    signs = np.sign(magmoms)
    signs[np.abs(magmoms) < moment_tol] = 0.0

    ops: list[SpinOp] = []
    for rot, trans in zip(sym["rotations"], sym["translations"]):
        mapped = _atom_map(rot, trans, scaled_positions, numbers, symprec)
        if mapped is None:
            continue  # numerically not a clean mapping; skip rather than guess
        s_new = signs[mapped]
        mag = np.abs(magmoms) >= moment_tol
        if not np.any(mag):
            raise ValueError("no magnetic atoms: nothing to classify")
        if np.all(s_new[mag] == signs[mag]):
            chi = +1
        elif np.all(s_new[mag] == -signs[mag]):
            chi = -1
        else:
            continue  # broken by the magnetic order in either sense
        ops.append(
            SpinOp(
                rot_frac=np.array(rot, dtype=int),
                trans=np.array(trans, dtype=float),
                rot_cart=_cartesian_rotation(rot, cell),
                chi=chi,
            )
        )

    if not ops:
        raise RuntimeError("no operation survived: check magmoms and symprec")

    if require_exchange and not any(op.chi == -1 for op in ops):
        net = float(np.sum(magmoms))
        if abs(net) < moment_tol:
            raise SymmetryToleranceError(
                "the moments cancel, but no operation of the parent exchanges the\n"
                "magnetic sublattices.  Without one the splitting is unconstrained and\n"
                "l = 0 comes out allowed, i.e. this would be reported as if it were a\n"
                "ferromagnet.  The usual cause is symprec: a relaxed structure carries\n"
                "numerical noise, and the sublattice-exchanging operation is often the\n"
                "first one spglib drops.  Run altharm.symmetry.symprec_scan() to see\n"
                f"where the operation count changes (current symprec = {symprec:g}).\n"
                "Pass require_exchange=False if the sublattices really are unrelated."
            )

    if include_time_reversal:
        ops = _close_with_time_reversal(ops)
    return ops


class SymmetryToleranceError(RuntimeError):
    """The symmetry found depends on symprec in a way that changes the physics."""


def symprec_scan(cell, scaled_positions, numbers, magmoms,
                 symprecs=(1e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 1e-2),
                 **kwargs):
    """How the operation count and its characters depend on ``symprec``.

    A relaxed structure is only approximately symmetric, so the group spglib
    reports is a function of the tolerance.  That is usually harmless, but for
    an altermagnet it is not: the sublattice-exchanging operation is the one
    that makes the material an altermagnet at all, and it is often the first
    casualty of a tight tolerance.

    Measured on the group's Cr2SeAs model, the parent is reported as P3m1 with
    six operations and *no* exchanging operation at symprec <= 1e-5, and as
    P-6m2 with six of each at symprec >= 5e-5.  The horizontal mirror maps every
    atom onto its partner to within 1.0e-5 Angstrom, which is exactly the
    default tolerance -- so the default silently returns the wrong physics.

    Returns
    -------
    list[dict]
        One row per tolerance: ``symprec``, ``spacegroup``, ``number``,
        ``n_preserving``, ``n_exchanging``, ``error``.
    """
    import spglib as _spglib

    rows = []
    for sp in symprecs:
        row = {"symprec": sp, "spacegroup": None, "number": None,
               "n_preserving": 0, "n_exchanging": 0, "error": None}
        try:
            ds = _spglib.get_symmetry_dataset(
                (np.asarray(cell, float), np.asarray(scaled_positions, float),
                 np.asarray(numbers, int)), symprec=sp)
            row["spacegroup"] = ds.international
            row["number"] = int(ds.number)
        except Exception as exc:                       # pragma: no cover
            row["error"] = str(exc)
        try:
            ops = spin_operations(cell, scaled_positions, numbers, magmoms,
                                  symprec=sp, include_time_reversal=False,
                                  require_exchange=False, **kwargs)
            row["n_preserving"] = sum(o.chi == +1 for o in ops)
            row["n_exchanging"] = sum(o.chi == -1 for o in ops)
        except Exception as exc:
            row["error"] = str(exc)
        rows.append(row)
    return rows


def format_symprec_scan(rows):
    """The scan as a table, for putting in a log or a record."""
    out = [f"{'symprec':>10} {'spacegroup':>12} {'no.':>5} {'chi=+1':>7} {'chi=-1':>7}"]
    for r in rows:
        out.append(
            f"{r['symprec']:10.0e} {str(r['spacegroup']):>12} "
            f"{str(r['number']):>5} {r['n_preserving']:7d} {r['n_exchanging']:7d}"
            + (f"   {r['error'].splitlines()[0]}" if r["error"] else "")
        )
    return "\n".join(out)


def _atom_map(rot, trans, scaled_positions, numbers, symprec):
    """Permutation induced on the atoms, or None if it is not clean."""
    new = scaled_positions @ np.asarray(rot, dtype=float).T + np.asarray(trans)
    n = len(scaled_positions)
    mapped = np.full(n, -1, dtype=int)
    for i in range(n):
        d = new[i] - scaled_positions
        d -= np.rint(d)                       # minimum image in fractional space
        dist = np.linalg.norm(d, axis=1)
        same = numbers == numbers[i]
        dist = np.where(same, dist, np.inf)
        j = int(np.argmin(dist))
        if dist[j] > max(symprec * 10, 1e-4):
            return None
        mapped[i] = j
    if len(set(mapped.tolist())) != n:
        return None                            # not a permutation
    return mapped


def _close_with_time_reversal(ops):
    """Append ``-R`` for every ``R``, with the character unchanged.

    Adding :math:`k \\to -k` with character ``+1`` and closing the set under it
    is what removes the odd harmonics.  Duplicates are dropped by comparing the
    Cartesian matrices.
    """
    out = list(ops)
    for op in ops:
        cart = -op.rot_cart
        if any(np.allclose(cart, o.rot_cart, atol=1e-8) and o.chi == op.chi for o in out):
            continue
        out.append(
            SpinOp(
                rot_frac=-op.rot_frac,
                trans=op.trans,
                rot_cart=cart,
                chi=op.chi,
                label="time-reversal",
            )
        )
    return out


def symmetrizer(basis, ops, dirs):
    """Matrix of the projector ``S`` in the given basis.

    ``dirs`` must sample the basis well enough to invert its design matrix --
    for a circular basis of ``lmax``, comfortably more than ``2*lmax+1`` points.
    Operations that do not map the sampled manifold to itself are dropped, since
    a function on that manifold cannot be evaluated at a point off it.
    """
    dirs = np.atleast_2d(np.asarray(dirs, dtype=float))
    design = basis.evaluate(dirs)

    kept = [op for op in ops if _preserves(basis, op, dirs)]
    if not kept:
        raise RuntimeError("no operation maps the sampled directions to themselves")

    acc = np.zeros_like(design)
    for op in kept:
        rinv = np.linalg.inv(op.rot_cart)
        acc += op.chi * basis.evaluate(dirs @ rinv.T)
    acc /= len(kept)

    smat, *_ = np.linalg.lstsq(design, acc, rcond=None)
    return smat, kept


def _preserves(basis, op, dirs):
    """Does the operation map the sampled manifold onto itself?

    For a circular basis the sampled directions all share one component along
    the plane normal.  A great circle through the origin has that component
    zero; a circle at fixed ``k_z`` is a cone and has it finite.  Both are
    manifolds of the form ``{d : d.n = c}``, and an operation preserves such a
    set exactly when it leaves that component alone.

    The distinction matters.  In Cr2SeAs the operation that exchanges the
    magnetic sublattices is the horizontal mirror ``M_z``, which maps the
    ``+k_z`` cone to the ``-k_z`` one and so does not constrain a single cone by
    itself -- but combined with time reversal it gives ``C_2z``, which does.
    Treating the cone as if it were a plane would throw both away.
    """
    moved = dirs @ np.linalg.inv(op.rot_cart).T
    moved /= np.linalg.norm(moved, axis=1, keepdims=True)
    if hasattr(basis, "e1"):  # circular
        normal = np.cross(basis.e1, basis.e2)
        d = np.asarray(dirs, dtype=float)
        d = d / np.linalg.norm(d, axis=1, keepdims=True)
        return np.allclose(moved @ normal, d @ normal, atol=1e-8)
    return True  # spherical: every rotation preserves the sphere


def allowed_harmonics(basis, ops, dirs, tol: float = 1e-6):
    """Allowed subspace, resolved by harmonic index ``l``.

    Returns
    -------
    dict
        ``l -> (k, m)`` array whose ``k`` rows are an orthonormal basis of the
        allowed combinations inside that ``l`` block, written in the basis
        functions of that block.  An ``l`` absent from the dict is forbidden
        outright.
    """
    smat, _ = symmetrizer(basis, ops, dirs)
    ls = basis.l_of()
    out = {}
    for l in sorted(set(ls.tolist())):
        idx = np.where(ls == l)[0]
        block = smat[np.ix_(idx, idx)]
        # S is a projector, so it is diagonalisable with eigenvalues 0 and 1.
        w, v = np.linalg.eigh(0.5 * (block + block.T))
        keep = v[:, w > 1.0 - tol].T
        if keep.size:
            out[int(l)] = keep
    return out


def allowed_l(basis, ops, dirs, tol: float = 1e-6):
    """The set of harmonic indices that carry at least one allowed function."""
    return set(allowed_harmonics(basis, ops, dirs, tol=tol).keys())
