"""Writing Quantum ESPRESSO ``pw.x`` inputs for fixed-k-list band runs.

Two things are easy to get silently wrong when driving QE from a script
instead of by hand, and both would corrupt a ~1 meV spin splitting before it
ever reaches the fitting code in :mod:`altharm.basis`:

* the k-point list has to arrive at the diagonalizer in the same order it
  was requested, so that "the n-th k-point" still means "angle theta_n"
  downstream (see :func:`bands_input` for how QE's own symmetry machinery
  fights this by default);
* the eigenvalue convergence threshold has to be tight enough that its noise
  does not swamp the splitting (see the ``diago_thr_init``/``conv_thr``
  discussion there too).

Everywhere a length or a k-vector crosses the boundary into QE's world it is
converted explicitly via :mod:`altharm.units`, never left for QE to
interpret in whatever unit it feels like defaulting to.
"""

from __future__ import annotations

import numpy as np

from ..units import BOHR_ANG


def cart_to_tpiba(k_invang, alat_bohr: float):
    """Cartesian k, 1/Angstrom -> QE ``tpiba`` units (2*pi/alat).

    ``tpiba`` is Cartesian, exactly like the input -- the only thing that
    changes is the length scale, which is why this conversion is safe to
    apply to an angularly-ordered k-list without disturbing the ordering.
    ``alat_bohr`` must be the exact ``celldm(1)`` used in the QE input (see
    :func:`bands_input`), since that is what fixes the scale of "tpiba" on
    the QE side.

    Parameters
    ----------
    k_invang : array_like, shape (..., 3)
        Cartesian k-vector(s), 1/Angstrom.
    alat_bohr : float
        QE's alat, Bohr (QE's native length unit -- this is why the
        Bohr->Angstrom factor from :mod:`altharm.units` is needed here and
        nowhere else in this module).

    Returns
    -------
    ndarray, same shape as ``k_invang``
        k in units of 2*pi/alat.
    """
    k_invang = np.asarray(k_invang, dtype=float)
    alat_ang = alat_bohr * BOHR_ANG
    two_pi_over_alat = 2.0 * np.pi / alat_ang  # 1/Angstrom
    return k_invang / two_pi_over_alat


def tpiba_to_cart(k_tpiba, alat_bohr: float):
    """Inverse of :func:`cart_to_tpiba`: ``tpiba`` units -> 1/Angstrom.

    Needed on the reader side, since QE's XML output reports k-points back
    in tpiba units regardless of what unit they were specified in.
    """
    k_tpiba = np.asarray(k_tpiba, dtype=float)
    alat_ang = alat_bohr * BOHR_ANG
    two_pi_over_alat = 2.0 * np.pi / alat_ang  # 1/Angstrom
    return k_tpiba * two_pi_over_alat


def kpoints_card(kcart_tpiba, weights=None) -> str:
    """Build a QE ``K_POINTS tpiba`` card from an explicit k-point list.

    ``tpiba`` k-points are Cartesian, in units of 2*pi/alat. That is the
    entire reason to use them here instead of ``crystal`` (fractional)
    coordinates: fractional coordinates carry the lattice metric, which for
    anything but a cubic cell stretches the circle of k-points this package
    fits harmonics on into an ellipse (see the module docstring of
    :mod:`altharm.basis`). ``tpiba`` bypasses the metric entirely.

    Parameters
    ----------
    kcart_tpiba : array_like, shape (nk, 3)
        k-points in tpiba units (2*pi/alat), Cartesian, in the order they
        should be computed and later read back in.
    weights : array_like, shape (nk,), optional
        QE weights. For a bands run these do not need to integrate to
        anything physical -- QE just wants a number per k-point -- so the
        default of 1.0 for every point is fine.

    Returns
    -------
    str
        A complete ``K_POINTS`` card, newline-terminated.
    """
    kcart_tpiba = np.atleast_2d(np.asarray(kcart_tpiba, dtype=float))
    if kcart_tpiba.ndim != 2 or kcart_tpiba.shape[1] != 3:
        raise ValueError(f"expected kcart_tpiba of shape (nk, 3), got {kcart_tpiba.shape}")
    nk = kcart_tpiba.shape[0]

    if weights is None:
        weights = np.ones(nk)
    else:
        weights = np.asarray(weights, dtype=float)
        if weights.shape != (nk,):
            raise ValueError(f"weights shape {weights.shape} != ({nk},)")

    lines = ["K_POINTS tpiba", str(nk)]
    for (kx, ky, kz), w in zip(kcart_tpiba, weights):
        lines.append(f"{kx:.10f} {ky:.10f} {kz:.10f} {w:.10f}")
    return "\n".join(lines) + "\n"


def bands_input(
    *,
    prefix: str,
    outdir: str,
    pseudo_dir: str,
    alat_bohr: float,
    cell_ang,
    species,
    positions_ang,
    ecutwfc: float,
    kcart_tpiba,
    weights=None,
    ecutrho: float | None = None,
    nbnd: int | None = None,
    nspin: int = 1,
    starting_magnetization=None,
    diago_thr_init: float = 1.0e-4,
    conv_thr: float = 1.0e-8,
) -> str:
    """Write a ``calculation='bands'`` pw.x input for a fixed, explicit k-list.

    This assumes ``outdir``/``prefix`` already holds a converged charge
    density from a prior ``scf`` run -- that is what makes a ``bands``
    calculation non-self-consistent: QE reads the density once and
    diagonalizes at the k-points given here without ever changing it.

    Parameters
    ----------
    alat_bohr : float
        Written as ``celldm(1)``. This fixes what "tpiba" means for the
        ``K_POINTS`` card below, so it must be the same value later passed
        to :func:`cart_to_tpiba`/:func:`tpiba_to_cart`.
    cell_ang : array_like, shape (3, 3)
        Lattice vectors as rows, Angstrom. Written out in ``alat`` units,
        divided by ``alat_bohr``, because QE stops with "lattice parameter
        specified twice" if ``celldm(1)`` and an explicitly-dimensioned
        ``CELL_PARAMETERS`` are both given. ``celldm(1)`` is the one kept, since
        it is what defines ``tpiba`` for the k-point card.
    species : sequence of (symbol, mass_amu, pseudo_filename)
    positions_ang : sequence of (symbol, x, y, z), Cartesian, Angstrom
    kcart_tpiba, weights : see :func:`kpoints_card`.
    nspin, starting_magnetization : LSDA controls; ``starting_magnetization``
        is one value per species, in the same order as ``species``.
    diago_thr_init, conv_thr : eigenvalue/SCF convergence thresholds -- see
        the comment below on why both matter even though this is not an SCF
        run.

    Returns
    -------
    str
        A complete pw.x input file, newline-terminated.
    """
    cell_ang = np.atleast_2d(np.asarray(cell_ang, dtype=float))
    if cell_ang.shape != (3, 3):
        raise ValueError(f"cell_ang must have shape (3, 3), got {cell_ang.shape}")
    species = list(species)
    positions_ang = list(positions_ang)
    ntyp = len(species)
    nat = len(positions_ang)
    if ntyp == 0:
        raise ValueError("species must contain at least one entry")
    if nat == 0:
        raise ValueError("positions_ang must contain at least one entry")
    if nspin == 2 and starting_magnetization is not None and len(starting_magnetization) != ntyp:
        raise ValueError("starting_magnetization must have one entry per species")

    lines = []

    lines.append("&CONTROL")
    lines.append("  calculation = 'bands'")
    lines.append("  verbosity = 'high'")
    lines.append(f"  prefix = '{prefix}'")
    lines.append(f"  outdir = '{outdir}'")
    lines.append(f"  pseudo_dir = '{pseudo_dir}'")
    lines.append("/")

    lines.append("&SYSTEM")
    lines.append("  ibrav = 0")
    lines.append(f"  celldm(1) = {alat_bohr:.12f}")
    lines.append(f"  nat = {nat}")
    lines.append(f"  ntyp = {ntyp}")
    lines.append(f"  ecutwfc = {ecutwfc:.6f}")
    if ecutrho is not None:
        lines.append(f"  ecutrho = {ecutrho:.6f}")
    if nbnd is not None:
        lines.append(f"  nbnd = {int(nbnd)}")
    if nspin == 2:
        lines.append("  nspin = 2")
        if starting_magnetization is not None:
            for i, m in enumerate(starting_magnetization, start=1):
                lines.append(f"  starting_magnetization({i}) = {m:.6f}")
    # nosym/noinv are not optional. QE's default behaviour for an explicit
    # k-point list -- including for calculation='bands' -- is to use the
    # crystal's symmetry operations (plus time reversal) to map k-points onto
    # each other and only actually diagonalize an irreducible subset,
    # re-expanding afterwards. In practice that means the eigenvalues QE
    # hands back are not in the order they were requested in: exactly the
    # ordering this package depends on to attach "the n-th k-point" to
    # "angle theta_n" around the circle. It is also not merely a bookkeeping
    # nuisance here: time reversal maps (k, spin up) to (-k, spin down) in a
    # collinear calculation, not (k, up) to (-k, up), so folding a spin
    # channel's k-points onto their star using ordinary spatial symmetry is
    # not even the right symmetry operation for what is being measured.
    # Setting both flags disables all of this and returns exactly the
    # k-points given, in the order given.
    lines.append("  nosym = .true.")
    lines.append("  noinv = .true.")
    lines.append("/")

    lines.append("&ELECTRONS")
    # diago_thr_init is the whole story for eigenvalue noise here: a 'bands'
    # run has no SCF loop to progressively tighten the diagonalization
    # threshold, so this single number is the noise floor under the
    # splitting for the entire run. conv_thr governs SCF energy convergence,
    # which is moot for a fixed-density bands run, but QE's namelist parser
    # still wants a sensible value, and leaving it loose from a previous scf
    # run in the same input template is a common way to silently accept a
    # noisier density than intended -- so it is set explicitly here too.
    lines.append(f"  diago_thr_init = {diago_thr_init:.3e}")
    lines.append(f"  conv_thr = {conv_thr:.3e}")
    lines.append("/")

    lines.append("ATOMIC_SPECIES")
    for symbol, mass, pseudo in species:
        lines.append(f"  {symbol} {mass:.6f} {pseudo}")

    lines.append("ATOMIC_POSITIONS angstrom")
    for symbol, x, y, z in positions_ang:
        lines.append(f"  {symbol} {x:.10f} {y:.10f} {z:.10f}")

    # QE refuses "lattice parameter specified twice": with celldm(1) present,
    # CELL_PARAMETERS must be in alat units, not angstrom.  Real pw.x 7.2 stops
    # at cell_base_init with that exact message otherwise.  We keep celldm(1)
    # rather than dropping it, because it is what fixes the meaning of "tpiba"
    # in the K_POINTS card below -- and an ambiguous tpiba would silently put
    # the circle at the wrong radius.
    lines.append("CELL_PARAMETERS alat")
    alat_ang = alat_bohr * BOHR_ANG
    for row in np.asarray(cell_ang, dtype=float) / alat_ang:
        lines.append(f"  {row[0]:.12f} {row[1]:.12f} {row[2]:.12f}")

    lines.append(kpoints_card(kcart_tpiba, weights=weights).rstrip("\n"))

    return "\n".join(lines) + "\n"
