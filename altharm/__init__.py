"""altharm -- angular harmonics of altermagnetic spin splitting, direct from DFT."""

from __future__ import annotations

from .basis import CircularBasis, SphericalBasis, real_sph_harm
from .kpath import circle, sphere_grid, angles_of
from .project import HarmonicFit, project, evaluate
from .symmetry import SpinOp, spin_operations, symmetrizer, allowed_harmonics, allowed_l
from .noise import NoiseReport, noise_report
from .symmetry import SymmetryToleranceError, symprec_scan, format_symprec_scan
from . import wannier, qe

__version__ = "0.0.1"

__all__ = [
    "CircularBasis", "SphericalBasis", "real_sph_harm",
    "circle", "sphere_grid", "angles_of",
    "HarmonicFit", "project", "evaluate",
    "SpinOp", "spin_operations", "symmetrizer", "allowed_harmonics", "allowed_l",
    "NoiseReport", "noise_report",
    "SymmetryToleranceError", "symprec_scan", "format_symprec_scan",
    "wannier", "qe",
]
