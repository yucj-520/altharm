"""One call from a finished DFT run to the harmonics of its spin splitting.

Everything this module does could already be done by hand with the pieces in
:mod:`altharm.basis`, :mod:`altharm.symmetry`, :mod:`altharm.project` and
:mod:`altharm.qe` -- and in the phases recorded in ``PLAN.md`` it was, one
driver script per run.  Assembling it here is not a convenience wrapper: every
one of those scripts had to re-decide the same four things, and three of them
are the traps this project exists to avoid.

* **The circle's plane.**  Read from the k-points, never assumed.  A circle in
  the y-z plane handed to a basis left in its default x-y frame returns angles,
  a design matrix and a plausible allowed set -- thirteen scattered ``l`` on
  Cr2SeAs where the correct frame gives the even ones -- and raises nothing.
* **The structure.**  Read from the same XML as the eigenvalues, so the cell,
  the positions and which sublattice is which cannot drift apart from the run
  they describe.
* **``symprec``.**  Scanned, not defaulted.  On Cr2SeAs spglib's default drops
  the mirror that makes the material an altermagnet, and the result is a report
  in which every ``l`` including ``l = 0`` is allowed.
* **Band identity.**  Reported against ``E_F`` with the gap to its neighbours,
  because "band 44" means a different state in a 38-band model, in an 80-band
  run, and at a different ``U``.

The number this prints next to each coefficient is the forbidden-harmonic noise
floor, which is an error bar on the *fit* and not on the coefficient.  See
``README.md`` -- the inputs to the Kohn-Sham Hamiltonian (``U``, the functional,
the pseudopotential, the smearing) move the coefficients by 7-21% without the
floor registering any of it.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .basis import CircularBasis
from .noise import NoiseReport, noise_report
from .project import HarmonicFit, project
from .qe.reader import QEBands, QEStructure, read_bands, read_structure
from .symmetry import (
    SpinOp,
    SymmetryToleranceError,
    allowed_l,
    format_symprec_scan,
    spin_operations,
    symprec_scan,
)
from .units import BOHR_ANG

__all__ = ["CircleGeometry", "HarmonicAnalysis", "circle_frame", "analyze_qe"]


@dataclass
class CircleGeometry:
    """The circle the k-points actually lie on, measured rather than assumed."""

    center: np.ndarray          # 1/Angstrom, Cartesian
    normal: np.ndarray          # unit, the circle's axis
    e1: np.ndarray              # unit, in the plane; theta = 0 at the first k-point
    e2: np.ndarray              # unit, = normal x e1
    radius: float               # 1/Angstrom, mean
    radius_spread: float        # 1/Angstrom, max - min; a sanity check on the mesh
    planarity: float            # 1/Angstrom, out-of-plane spread of the points

    @property
    def kz_over_radius(self) -> float:
        """Offset of the circle's plane from the origin, in units of its radius."""
        return float(self.center @ self.normal) / self.radius


def circle_frame(kcart) -> CircleGeometry:
    """Fit the plane of a set of k-points and return a frame spanning it.

    The frame is what :class:`~altharm.basis.CircularBasis` has to be built in.
    ``e1`` is aimed at the first k-point so that ``theta = 0`` there, which makes
    the cos/sin split reproducible across runs that requested the same circle.
    """
    k = np.atleast_2d(np.asarray(kcart, dtype=float))
    if k.shape[0] < 3:
        raise ValueError(f"need at least 3 k-points to fit a circle, got {k.shape[0]}")

    center = k.mean(axis=0)
    rel = k - center
    # Smallest right-singular vector of the centred points is the plane normal.
    _, svals, vt = np.linalg.svd(rel, full_matrices=False)
    normal = vt[-1]
    # Deterministic orientation: largest-magnitude component positive.
    if normal[int(np.argmax(np.abs(normal)))] < 0:
        normal = -normal

    e1 = rel[0] - (rel[0] @ normal) * normal
    n1 = np.linalg.norm(e1)
    if n1 < 1e-12:
        raise ValueError("the first k-point sits on the circle's axis; cannot set theta = 0")
    e1 = e1 / n1
    e2 = np.cross(normal, e1)

    radii = np.linalg.norm(rel - np.outer(rel @ normal, normal), axis=1)
    return CircleGeometry(
        center=center,
        normal=normal,
        e1=e1,
        e2=e2,
        radius=float(radii.mean()),
        radius_spread=float(np.ptp(radii)),
        planarity=float(np.ptp(rel @ normal)),
    )


@dataclass
class HarmonicAnalysis:
    """The result of :func:`analyze_qe`, with the provenance needed to read it."""

    source: str
    bands: QEBands = field(repr=False)
    structure: QEStructure = field(repr=False)
    geometry: CircleGeometry
    basis: CircularBasis = field(repr=False)
    fit: HarmonicFit = field(repr=False)
    noise: NoiseReport
    allowed: set
    ops: list = field(repr=False)
    symprec: float
    symprec_rows: list = field(repr=False)
    band: int
    band_energy_rel_fermi: float
    band_gap_below: float   # closest approach to the band below, eV
    band_gap_above: float   # closest approach to the band above, eV
    lmax: int
    dE: np.ndarray = field(repr=False)

    def amplitudes_mev(self) -> dict:
        """Amplitude per ``l``, meV, for the symmetry-allowed ``l`` only."""
        amp = self.fit.amplitude_by_l()
        return {l: amp[l] / 1e-3 for l in sorted(self.allowed)}

    def summary(self) -> str:
        g = self.geometry
        n_pres = sum(o.chi == +1 for o in self.ops)
        n_exch = sum(o.chi == -1 for o in self.ops)
        allowed = sorted(self.allowed)
        lines = [
            f"source                 : {self.source}",
            f"k-points on the circle : {len(self.dE)}",
            f"circle radius          : {g.radius:.6f} 1/A  "
            f"(spread {g.radius_spread:.2e}, planarity {g.planarity:.2e})",
            f"circle axis            : [{g.normal[0]:+.4f} {g.normal[1]:+.4f} "
            f"{g.normal[2]:+.4f}], offset {g.kz_over_radius:+.4f} x radius",
            f"band                   : {self.band}  "
            f"(E - E_F = {self.band_energy_rel_fermi:+.4f} eV, "
            f"closest approach below {self.band_gap_below * 1e3:.1f} meV, "
            f"above {self.band_gap_above * 1e3:.1f} meV)",
            f"symprec                : {self.symprec:g}",
            f"operations kept        : {len(self.ops)}  "
            f"({n_pres} sublattice-preserving, {n_exch} exchanging)",
            f"allowed l (<= {self.lmax})     : {allowed if allowed else 'nothing allowed'}",
            f"|dE| max               : {np.abs(self.dE).max() / 1e-3:.4f} meV",
            "",
            self.noise.summary(unit_name="meV", unit_scale=1e-3),
        ]
        return "\n".join(lines)


def _pick_band(qb: QEBands, band) -> int:
    nb = min(qb.nbnd_up, qb.nbnd_dw)
    if band in ("fermi", None):
        mean_up = qb.energies_ev[0, :, :nb].mean(axis=0)
        return int(np.argmin(np.abs(mean_up - qb.fermi_ev)))
    band = int(band)
    if not 0 <= band < nb:
        raise ValueError(
            f"band {band} is out of range: this run has {nb} bands in both spin "
            "channels (0-based indexing)"
        )
    return band


def _min_gap(qb: QEBands, lower: int, upper: int, nb: int) -> float:
    """Closest approach of two adjacent bands over the circle, either spin, eV.

    Returns ``inf`` when the pair falls outside the band range.
    """
    if lower < 0 or upper >= nb:
        return float("inf")
    gaps = qb.energies_ev[:, :, upper] - qb.energies_ev[:, :, lower]
    return float(np.nanmin(gaps))


def _resolve_symprec(structure: QEStructure, symprec, require_exchange: bool):
    """Scan symprec, and refuse a value whose answer is not stable around it.

    Cr2SeAs is the case this exists for: the sublattice-exchanging mirror maps
    atom onto atom to within 1.0e-5 A, exactly spglib's default tolerance, so
    the default silently describes an altermagnet as something that cannot be
    one.  Scanning is cheap and the failure it catches is silent.
    """
    rows = symprec_scan(
        structure.cell_ang, structure.scaled_positions, structure.numbers, structure.magmoms
    )
    if symprec is not None:
        return float(symprec), rows

    # Prefer the smallest tolerance that still finds an exchanging operation:
    # looser tolerances can only merge more atoms, and merging too many invents
    # symmetry the structure does not have.
    usable = [r for r in rows if r.get("n_exchanging", 0) > 0] if rows else []
    if not usable and require_exchange:
        raise SymmetryToleranceError(
            "no scanned symprec finds a sublattice-exchanging operation, so this "
            "structure cannot be classified as an altermagnet.\n\n"
            + format_symprec_scan(rows)
        )
    if not usable:
        return 1e-5, rows
    chosen = min(r["symprec"] for r in usable)
    return float(chosen), rows


def analyze_qe(
    xml_path,
    band="fermi",
    lmax=None,
    symprec=None,
    include_time_reversal: bool = True,
    require_exchange: bool = True,
    structure: QEStructure = None,
) -> HarmonicAnalysis:
    """Harmonics of the spin splitting, from a finished QE bands run.

    Parameters
    ----------
    xml_path : path
        ``<outdir>/<prefix>.save/data-file-schema.xml`` from a spin-polarised
        ``calculation = 'bands'`` run over a circle of k-points.  The run must
        have used ``nosym``/``noinv``, or QE will have folded the circle into
        the irreducible wedge and the k-points will no longer be the ones asked
        for -- :mod:`altharm.qe.writer` sets both.
    band : int or ``"fermi"``
        0-based band index, or ``"fermi"`` for the band whose mean energy is
        closest to ``E_F``.  **A band index is not portable** between models,
        between band counts, or between values of ``U``; the returned object
        reports the band's energy relative to ``E_F`` so the choice can be
        checked rather than trusted.
    lmax : int, optional
        Highest harmonic fitted.  Defaults to ``n_kpoints // 4``, which keeps
        the fit over-determined by about 2:1 and stays a factor of two below the
        Nyquist index ``n_kpoints // 2``.  Going closer to Nyquist aliases the
        signal's own high-``l`` tail into the forbidden harmonics and inflates
        the noise floor -- measured in ``PLAN.md`` P3.
    symprec : float, optional
        Passed to spglib.  Default scans and picks; see
        :func:`altharm.symmetry.symprec_scan`.
    include_time_reversal : bool
        Leave ``True`` for a run without spin-orbit coupling.  **Set it False
        for an SOC run**, where ``E(k) = E(-k)`` per spin channel no longer
        holds and odd harmonics are not forbidden.
    structure : QEStructure, optional
        Override the structure read from the XML -- for the rare run whose
        species labels or ``starting_magnetization`` do not describe the
        magnetic order actually converged to.

    Returns
    -------
    HarmonicAnalysis
    """
    qb = read_bands(xml_path)
    if not qb.lsda:
        raise ValueError(
            f"{xml_path} is not a spin-polarised run, so there is no spin splitting "
            "to expand (nspin = 1)"
        )
    st = structure if structure is not None else read_structure(xml_path)

    # tpiba -> 1/Angstrom.  alat is in Bohr; 2*pi/alat is the tpiba unit.
    tpiba = 2.0 * np.pi / (qb.alat_bohr * BOHR_ANG)
    kcart = qb.kpoints_tpiba * tpiba

    geom = circle_frame(kcart)
    n = kcart.shape[0]
    if lmax is None:
        lmax = max(1, n // 4)
    lmax = int(lmax)
    if 2 * lmax + 1 > n:
        raise ValueError(
            f"lmax = {lmax} needs {2 * lmax + 1} coefficients but there are only "
            f"{n} k-points; the fit would be under-determined"
        )

    sp, rows = _resolve_symprec(st, symprec, require_exchange)
    ops = spin_operations(
        st.cell_ang, st.scaled_positions, st.numbers, st.magmoms,
        symprec=sp,
        include_time_reversal=include_time_reversal,
        require_exchange=require_exchange,
    )

    basis = CircularBasis(lmax, e1=geom.e1, e2=geom.e2)
    dirs = kcart / np.linalg.norm(kcart, axis=1)[:, None]
    allowed = allowed_l(basis, ops, dirs)

    b = _pick_band(qb, band)
    nb = min(qb.nbnd_up, qb.nbnd_dw)
    mean_up = qb.energies_ev[0, :, :nb].mean(axis=0)
    # The gap that decides band identity is the CLOSEST APPROACH over the
    # circle, in either spin channel -- not the separation of the mean
    # energies, which stays comfortably large while the bands touch at one
    # angle and swap character there. Cr2SeAs's bands 43/44 are 171 meV apart
    # on average and 12-19 meV apart at their closest.
    gap_below = _min_gap(qb, b - 1, b, nb)
    gap_above = _min_gap(qb, b, b + 1, nb)

    dE = qb.energies_ev[0, :, b] - qb.energies_ev[1, :, b]
    fit = project(basis, dirs, dE)

    return HarmonicAnalysis(
        source=str(xml_path),
        bands=qb,
        structure=st,
        geometry=geom,
        basis=basis,
        fit=fit,
        noise=noise_report(fit, allowed),
        allowed=allowed,
        ops=ops,
        symprec=sp,
        symprec_rows=rows,
        band=b,
        band_energy_rel_fermi=float(mean_up[b] - qb.fermi_ev),
        band_gap_below=gap_below,
        band_gap_above=gap_above,
        lmax=lmax,
        dE=dE,
    )
