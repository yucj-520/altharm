"""Parsing eigenvalues out of Quantum ESPRESSO's XML output.

``pw.x`` writes ``<outdir>/<prefix>.save/data-file-schema.xml``. Every energy
in that file -- eigenvalues, the Fermi level -- is in Hartree, QE's internal
atomic unit, and every k-point under ``ks_energies`` is Cartesian, in units
of 2*pi/alat (``tpiba``), matching what :mod:`altharm.qe.writer` sent in.
Both conversions go through :mod:`altharm.units` explicitly; nothing here
guesses a unit, because a silent factor-of-two or Hartree/eV slip would look
exactly like the ~1 meV altermagnetic splitting this package is trying to
measure.

The XML schema itself is a moving target across QE versions in ways that
don't matter for this module (extra elements, attributes, a namespace on the
root that its children may or may not inherit). Rather than depend on that
namespace behaviour, every lookup below matches on the local (unprefixed)
element name only.
"""

from __future__ import annotations

from dataclasses import dataclass
import xml.etree.ElementTree as ET

import numpy as np

from ..units import HARTREE_EV, BOHR_ANG


def _local(tag: str) -> str:
    """Strip any ``{namespace}`` prefix ElementTree attaches to a tag."""
    return tag.rsplit("}", 1)[-1]


def _child(elem, name):
    """First direct child of ``elem`` whose local name is ``name``, or None."""
    for c in elem:
        if _local(c.tag) == name:
            return c
    return None


def _require(elem, name, xml_path, where):
    c = _child(elem, name)
    if c is None:
        raise ValueError(
            f"<{name}> not found under <{where}> in {xml_path}; refusing to "
            "guess -- check that this is a pw.x data-file-schema.xml and "
            "that the run actually completed."
        )
    return c


def _floats(elem):
    return np.array([float(x) for x in elem.text.split()])


@dataclass
class QEBands:
    """Band structure read back from a QE ``bands``/``nscf`` run.

    Attributes
    ----------
    kpoints_tpiba : ndarray, shape (nk, 3)
        Cartesian k-points, units of 2*pi/alat, in the order QE wrote them
        -- the order requested, provided the run used ``nosym``/``noinv``
        (see :mod:`altharm.qe.writer`).
    energies_ev : ndarray, shape (nspin, nk, nbnd)
        eV. ``nspin`` is 1 for a non-spin-polarised run, 2 (up, then down)
        for LSDA. If ``nbnd_up != nbnd_dw``, ``nbnd`` here is
        ``max(nbnd_up, nbnd_dw)`` and the shorter spin channel is padded
        with NaN, rather than silently truncating the longer one.
    fermi_ev : float
        eV.
    alat_bohr : float
        QE's alat, Bohr -- needed to undo the tpiba scaling of the k-points.
    cell_ang : ndarray, shape (3, 3)
        Real-space lattice vectors as rows, Angstrom (converted from the
        Bohr the XML stores them in).
    lsda : bool
        Whether the run was spin-polarised.
    nbnd_up, nbnd_dw : int
        True band counts per spin channel (equal, and equal to the single
        ``nbnd``, for a non-spin-polarised run). Use these, not
        ``energies_ev.shape[-1]``, to know how far a channel is real rather
        than NaN padding.
    """

    kpoints_tpiba: np.ndarray
    energies_ev: np.ndarray
    fermi_ev: float
    alat_bohr: float
    cell_ang: np.ndarray
    lsda: bool
    nbnd_up: int
    nbnd_dw: int


def read_bands(xml_path) -> QEBands:
    """Parse ``data-file-schema.xml`` into a :class:`QEBands`.

    Raises
    ------
    ValueError
        If any element this function relies on is missing, or the file's
        internal bookkeeping (``nks`` vs. the number of ``ks_energies``
        blocks, eigenvalue array lengths vs. declared band counts) is
        inconsistent. A silent fallback here -- returning None, or padding
        without complaint -- would corrupt a physics result rather than
        just fail loudly, so every assumption is checked.
    """
    xml_path = str(xml_path)
    tree = ET.parse(xml_path)
    root = tree.getroot()

    output = _require(root, "output", xml_path, _local(root.tag))
    atomic_structure = _require(output, "atomic_structure", xml_path, "output")
    band_structure = _require(output, "band_structure", xml_path, "output")

    alat_attr = atomic_structure.get("alat")
    if alat_attr is None:
        raise ValueError(f"<atomic_structure> has no 'alat' attribute in {xml_path}")
    alat_bohr = float(alat_attr)

    cell = _require(atomic_structure, "cell", xml_path, "atomic_structure")
    cell_bohr = np.array(
        [_floats(_require(cell, name, xml_path, "cell")) for name in ("a1", "a2", "a3")]
    )
    cell_ang = cell_bohr * BOHR_ANG

    lsda_elem = _require(band_structure, "lsda", xml_path, "band_structure")
    lsda = lsda_elem.text.strip().lower() == "true"

    if lsda:
        nbnd_up = int(_require(band_structure, "nbnd_up", xml_path, "band_structure").text)
        nbnd_dw = int(_require(band_structure, "nbnd_dw", xml_path, "band_structure").text)
    else:
        nbnd = int(_require(band_structure, "nbnd", xml_path, "band_structure").text)
        nbnd_up = nbnd_dw = nbnd

    nks = int(_require(band_structure, "nks", xml_path, "band_structure").text)

    fermi_elem = _child(band_structure, "fermi_energy")
    if fermi_elem is not None:
        fermi_ev = float(fermi_elem.text) * HARTREE_EV
    else:
        # A fixed-spin-moment metallic run reports two Fermi levels instead
        # of one; average them rather than refusing outright, since nothing
        # downstream of this module actually depends on fermi_ev.
        two_fermi = _child(band_structure, "two_fermi_energies")
        if two_fermi is None:
            raise ValueError(
                "neither <fermi_energy> nor <two_fermi_energies> found under "
                f"<band_structure> in {xml_path}"
            )
        fermi_ev = float(np.mean(_floats(two_fermi))) * HARTREE_EV

    ks_list = [c for c in band_structure if _local(c.tag) == "ks_energies"]
    if len(ks_list) != nks:
        raise ValueError(
            f"<band_structure><nks> says {nks} but {len(ks_list)} "
            f"<ks_energies> elements are present in {xml_path}"
        )
    if nks == 0:
        raise ValueError(f"<band_structure> has zero <ks_energies> elements in {xml_path}")

    nspin = 2 if lsda else 1
    nbnd_max = max(nbnd_up, nbnd_dw)
    kpoints = np.empty((nks, 3))
    energies = np.full((nspin, nks, nbnd_max), np.nan)

    for ik, ks in enumerate(ks_list):
        k_point = _require(ks, "k_point", xml_path, "ks_energies")
        kvec = _floats(k_point)
        if kvec.size != 3:
            raise ValueError(
                f"<k_point> at ks_energies[{ik}] has {kvec.size} components, expected 3 "
                f"({xml_path})"
            )
        kpoints[ik] = kvec

        eig_elem = _require(ks, "eigenvalues", xml_path, "ks_energies")
        eig = _floats(eig_elem)
        expected = nbnd_up + nbnd_dw if lsda else nbnd_up
        if eig.size != expected:
            raise ValueError(
                f"<eigenvalues> at ks_energies[{ik}] has {eig.size} values, "
                f"expected {expected} ({xml_path})"
            )
        if lsda:
            # Up-spin values first, then down-spin -- the two blocks are
            # concatenated, not interleaved. Getting this backwards silently
            # negates the splitting rather than erroring, which is exactly
            # the failure mode this module is written to avoid.
            energies[0, ik, :nbnd_up] = eig[:nbnd_up] * HARTREE_EV
            energies[1, ik, :nbnd_dw] = eig[nbnd_up:] * HARTREE_EV
        else:
            energies[0, ik, :nbnd_up] = eig * HARTREE_EV

    return QEBands(
        kpoints_tpiba=kpoints,
        energies_ev=energies,
        fermi_ev=fermi_ev,
        alat_bohr=alat_bohr,
        cell_ang=cell_ang,
        lsda=lsda,
        nbnd_up=nbnd_up,
        nbnd_dw=nbnd_dw,
    )


def spin_splitting(qebands: QEBands, band_index: int) -> np.ndarray:
    """dE(k) = E_up(k) - E_dn(k) for one band, eV, shape (nk,).

    This is the raw quantity the rest of the package expands in cos(l*theta)
    harmonics. It is only defined for a spin-polarised run: a
    non-spin-polarised calculation has no up/down channels to subtract, so
    that case is an error, not a silent zero.
    """
    if not qebands.lsda:
        raise ValueError(
            "spin_splitting requires a spin-polarised (LSDA) run; this "
            "QEBands has lsda=False, so there is no up/down difference to take"
        )
    nbnd_common = min(qebands.nbnd_up, qebands.nbnd_dw)
    if not (0 <= band_index < nbnd_common):
        raise ValueError(
            f"band_index={band_index} out of range: nbnd_up={qebands.nbnd_up}, "
            f"nbnd_dw={qebands.nbnd_dw} (only bands 0..{nbnd_common - 1} exist in both)"
        )
    return qebands.energies_ev[0, :, band_index] - qebands.energies_ev[1, :, band_index]
