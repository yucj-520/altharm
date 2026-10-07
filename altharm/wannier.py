"""The Wannier route: harmonics from a tight-binding Hamiltonian.

This is the baseline the DFT route is checked against, and it is also the object
of study in its own right.  A Wannier Hamiltonian is a finite sum

.. math:: H(k) = \\sum_R e^{i k \\cdot R} H(R) / \\deg(R)

over the lattice vectors kept in ``*_hr.dat``.  A harmonic of angular index
``l`` needs hoppings reaching far enough in ``R`` to produce that many
oscillations around the circle, so truncating the sum damps high ``l`` first.
That is the concrete sense in which a Wannier model can miss the high-harmonic
physics, and :func:`truncate` exists to measure how much: shrink the hopping
range, refit the harmonics, and watch which ones decay.

The Hamiltonian is real-space and finite, so nothing here has a k-mesh error --
unlike the DFT route, whose limitation is noise rather than truncation.  The two
routes therefore fail in different ways, which is what makes agreement between
them meaningful.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class WannierHam:
    """A ``seedname_hr.dat`` Hamiltonian.

    Attributes
    ----------
    rvecs : (nr,3) int array
        Lattice vectors, in units of the lattice vectors themselves.
    deg : (nr,) int array
        Wigner-Seitz degeneracies.  Wannier90 writes them and they must divide
        the corresponding ``H(R)``; forgetting that misweights the long-range
        hoppings, which is exactly where the high harmonics live.
    ham : (nr, nw, nw) complex array
    cell : (3,3) array or None
        Lattice vectors as ROWS, in Angstrom.  Needed to turn a Cartesian k into
        the fractional one the phase factor wants, and to measure ``|R|``.
    """

    rvecs: np.ndarray
    deg: np.ndarray
    ham: np.ndarray
    cell: np.ndarray | None = None

    @property
    def num_wann(self) -> int:
        return self.ham.shape[1]

    @property
    def num_rpts(self) -> int:
        return self.ham.shape[0]

    def rcart(self):
        """Lattice vectors in Cartesian Angstrom."""
        if self.cell is None:
            raise ValueError("cell is unknown; pass it to read_hr or set it")
        return self.rvecs @ np.asarray(self.cell, dtype=float)


def read_hr(path, cell=None) -> WannierHam:
    """Read a Wannier90 ``seedname_hr.dat``."""
    path = Path(path)
    with open(path) as fh:
        fh.readline()                                  # date stamp
        num_wann = int(fh.readline())
        nrpts = int(fh.readline())
        deg = []
        while len(deg) < nrpts:
            deg.extend(int(x) for x in fh.readline().split())
        body = np.fromstring(fh.read(), sep=" ")

    expected = nrpts * num_wann * num_wann * 7
    if body.size != expected:
        raise ValueError(
            f"{path}: expected {expected} numbers in the body "
            f"({nrpts} R-vectors x {num_wann}^2 x 7 columns), found {body.size}"
        )
    body = body.reshape(nrpts, num_wann * num_wann, 7)

    rvecs = body[:, 0, 0:3].astype(int)
    # column order is Rx Ry Rz i j Re Im, with i the fastest index
    ham = (body[:, :, 5] + 1j * body[:, :, 6]).reshape(nrpts, num_wann, num_wann)
    ham = np.transpose(ham, (0, 2, 1))                 # -> H[R, i, j]

    return WannierHam(
        rvecs=rvecs,
        deg=np.asarray(deg[:nrpts], dtype=int),
        ham=ham,
        cell=None if cell is None else np.asarray(cell, dtype=float),
    )


def hk(wh: WannierHam, kfrac):
    """``H(k)`` at fractional k-points, shape ``(nk, nw, nw)``."""
    kfrac = np.atleast_2d(np.asarray(kfrac, dtype=float))
    phase = np.exp(2j * np.pi * (kfrac @ wh.rvecs.T))          # (nk, nr)
    weighted = wh.ham / wh.deg[:, None, None]
    return np.einsum("kr,rij->kij", phase, weighted)


def cart_to_frac_k(kcart, cell):
    """Cartesian k in 1/Angstrom -> fractional, for a cell whose ROWS are a_i.

    With ``B = 2 pi (A^-1)^T`` holding the reciprocal vectors as rows,
    ``k_cart = k_frac B``, so ``k_frac = k_cart B^-1``.  Doing this conversion
    at the last possible moment is deliberate: the circle was built in Cartesian
    axes precisely so the lattice metric never reached it.
    """
    A = np.asarray(cell, dtype=float)
    B = 2.0 * np.pi * np.linalg.inv(A).T
    return np.atleast_2d(np.asarray(kcart, dtype=float)) @ np.linalg.inv(B)


def eigenvalues(wh: WannierHam, kfrac):
    """Eigenvalues in eV, sorted ascending, shape ``(nk, nw)``."""
    H = hk(wh, kfrac)
    H = 0.5 * (H + np.conjugate(np.transpose(H, (0, 2, 1))))   # enforce hermiticity
    return np.linalg.eigvalsh(H)


def spin_splitting(up: WannierHam, dn: WannierHam, kfrac, band):
    """``E_up - E_dn`` for one band index, in eV.

    Band index, not band identity: eigenvalues come out sorted by energy, so
    where two bands cross this difference jumps.  Check the result against
    :func:`altharm.project.HarmonicFit.residual_rms` -- a crossing shows up
    there as a residual the smooth basis cannot absorb.
    """
    if up.num_wann != dn.num_wann:
        raise ValueError(
            f"the two spin channels have different sizes ({up.num_wann}, {dn.num_wann})"
        )
    e_up = eigenvalues(up, kfrac)
    e_dn = eigenvalues(dn, kfrac)
    if not (0 <= band < e_up.shape[1]):
        raise IndexError(f"band {band} out of range 0..{e_up.shape[1] - 1}")
    return e_up[:, band] - e_dn[:, band]


def truncate(wh: WannierHam, rmax_ang: float, plane_normal=None) -> WannierHam:
    """Drop every hopping beyond ``rmax_ang``, in Angstrom.

    The tool for the truncation study: the high harmonics should be the first to
    go, and how fast they go is the measure of what a Wannier model of that
    range can and cannot represent.

    Parameters
    ----------
    plane_normal : array_like, shape (3,), optional
        Measure ``|R|`` only in the plane perpendicular to this vector, keeping
        every out-of-plane reach intact.

        **Use this whenever the circle lies in a plane.**  An angular harmonic
        on a circle in the ``k_x k_y`` plane is resolved by in-plane hoppings;
        the out-of-plane ones set the ``k_z`` dispersion and have nothing to do
        with it.  A plain isotropic cut confounds the two, and for a cell with a
        long c axis it does so violently: in Cr2Se2As, ``c = 16.9 A``, so any
        isotropic cut below that removes every ``R_z != 0`` hopping at once, the
        model loses its ``k_z`` dependence entirely, and a splitting that lives
        only at ``k_z != 0`` collapses to zero.  That looks like the harmonics
        dying of truncation when it is really the c axis being amputated.
    """
    if wh.cell is None:
        raise ValueError("cell is required to measure |R|; pass it to read_hr")
    rcart = wh.rcart()
    if plane_normal is not None:
        n = np.asarray(plane_normal, dtype=float)
        n = n / np.linalg.norm(n)
        rcart = rcart - np.outer(rcart @ n, n)
    keep = np.linalg.norm(rcart, axis=1) <= float(rmax_ang)
    if not np.any(keep):
        raise ValueError(f"rmax = {rmax_ang} A keeps no hoppings at all")
    return WannierHam(
        rvecs=wh.rvecs[keep],
        deg=wh.deg[keep],
        ham=wh.ham[keep],
        cell=wh.cell,
    )


def hopping_range(wh: WannierHam):
    """Largest ``|R|`` present, in Angstrom -- the model's reach."""
    return float(np.linalg.norm(wh.rcart(), axis=1).max())
