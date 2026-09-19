"""How many plane waves a k-point gets, and where that leaks into the fit.

At fixed ``ecutwfc``, the number of plane waves satisfying
:math:`|\\mathbf k+\\mathbf G|^2 \\le E_\\mathrm{cut}` is a step function of
:math:`\\mathbf k`: it changes by one exactly when a lattice point crosses the
cutoff sphere. Walking a circle in k-space therefore does not sample a smooth
basis -- the basis itself changes size, discretely, at a set of angles set by
the lattice rather than by the physics. That discontinuity is broadband in
:math:`\\theta` (a step has content at every harmonic), and it is a genuine
source of the numerical noise this package's whole strategy is built to catch
-- see :mod:`altharm.noise`.

:func:`count_planewaves` answers "how many plane waves does QE give this
k-point" **without calling QE**, by a direct sphere count. It exists to (a)
predict where the discontinuities fall before spending a DFT calculation, and
(b) cross-check :func:`altharm.qe.reader` indirectly: an eigenvalue run and
this count are two independent readings of the same underlying G-vector list,
computed by different code, and they must place the jumps at the same angles.
"""

from __future__ import annotations

import re

import numpy as np

from .units import BOHR_ANG


def _reciprocal(cell_ang):
    """Rows are b_i, in 1/Angstrom, with b_i . a_j = 2 pi delta_ij."""
    return 2.0 * np.pi * np.linalg.inv(np.asarray(cell_ang, dtype=float)).T


def count_planewaves(kcart_invang, cell_ang, ecutwfc_ry, margin: int = 2):
    """Number of G with :math:`|\\mathbf k+\\mathbf G|^2 \\le E_\\mathrm{cut}`, per k-point.

    A direct sphere count in Rydberg atomic units (where an inverse-Bohr
    wavevector squared is an energy in Ry), enumerating lattice points by
    brute force over a bounding box sized from the exact relation
    :math:`n_j = \\mathbf G \\cdot \\mathbf a_j / 2\\pi`, so the box is never
    too small no matter how skewed the cell is.

    Parameters
    ----------
    kcart_invang : (n, 3) array
        k-points, Cartesian, 1/Angstrom (this package's convention throughout).
    cell_ang : (3, 3) array
        Lattice vectors as rows, Angstrom.
    ecutwfc_ry : float
        The wavefunction cutoff, Rydberg -- QE's ``ecutwfc``.
    margin : int
        Extra shells added to the bounding box, guarding the derivation's
        rounding against points sitting exactly on the cutoff sphere.

    Returns
    -------
    (n,) int array
    """
    kcart_invang = np.atleast_2d(np.asarray(kcart_invang, dtype=float))
    cell_ang = np.asarray(cell_ang, dtype=float)
    B = _reciprocal(cell_ang)
    a_norms = np.linalg.norm(cell_ang, axis=1)

    # R_cut in 1/Angstrom: |k+G|_Bohr^2 <= ecut  <=>  |k+G|_Ang <= sqrt(ecut)/BOHR_ANG
    r_cut = np.sqrt(ecutwfc_ry) / BOHR_ANG
    k_max = np.linalg.norm(kcart_invang, axis=1).max() if kcart_invang.size else 0.0

    # n_j = G.a_j / (2 pi); |G.a_j| <= |G||a_j| <= (r_cut + k_max) |a_j|
    n_max = np.ceil((r_cut + k_max) * a_norms / (2.0 * np.pi)).astype(int) + margin
    grids = [np.arange(-m, m + 1) for m in n_max]
    N1, N2, N3 = np.meshgrid(*grids, indexing="ij")
    n = np.stack([N1.ravel(), N2.ravel(), N3.ravel()], axis=1)
    G = n @ B                                   # (ngrid, 3)

    counts = np.empty(len(kcart_invang), dtype=int)
    for i, k in enumerate(kcart_invang):
        counts[i] = int(np.sum(np.sum((k + G) ** 2, axis=1) <= r_cut ** 2))
    return counts


def parse_pw_count(text) -> np.ndarray:
    """The PW count QE itself printed, in k-point order, from ``pw.x`` stdout.

    Matches the ``( N PWs)`` QE prints next to each k-point in a ``bands`` or
    ``scf`` run. Independent of :mod:`altharm.qe.reader`, which never touches
    the text output -- the two are meant to be compared, not to share code.
    """
    return np.array([int(m) for m in re.findall(r"\(\s*(\d+)\s*PWs\)", text)])


def discontinuity_spectrum(basis, dirs, npw):
    """The angular harmonics carried by the plane-wave-count step function.

    Fits ``npw(theta)`` on the same basis used for the spin splitting, via
    :func:`altharm.project.project`. A jump of :math:`\\delta N` plane waves
    contaminates every harmonic; this says how much of it lands on each
    :math:`l`, so that a peak in ``noise_report``'s forbidden floor at
    the same :math:`l` can be attributed to this mechanism rather than left
    unexplained.
    """
    from .project import project

    return project(basis, dirs, np.asarray(npw, dtype=float))
