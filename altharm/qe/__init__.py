"""Quantum ESPRESSO I/O: writing bands-run inputs, reading eigenvalues back.

This subpackage is the boundary between the rest of ``altharm`` -- Cartesian
directions, eV, Angstrom, per the conventions in :mod:`altharm.basis` and
:mod:`altharm.units` -- and QE's own conventions (Hartree energies,
alat/tpiba-scaled lengths and k-points, Bohr). All conversions happen here,
explicitly, via :mod:`altharm.units`.

See :mod:`altharm.qe.writer` for building the pw.x input (and why
``nosym``/``noinv`` are mandatory) and :mod:`altharm.qe.reader` for parsing
``data-file-schema.xml`` back into eigenvalues.
"""

from __future__ import annotations

from .writer import kpoints_card, cart_to_tpiba, tpiba_to_cart, bands_input
from .reader import hubbard_is_active, require_hubbard, QEBands, read_bands, spin_splitting

__all__ = [
    "kpoints_card",
    "cart_to_tpiba",
    "tpiba_to_cart",
    "bands_input",
    "QEBands",
    "read_bands", "hubbard_is_active", "require_hubbard",
    "spin_splitting",
]
