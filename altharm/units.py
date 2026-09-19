"""Physical constants and unit conversions.

Single source of truth: every conversion factor used anywhere in the package is
defined here, with the CODATA value it comes from.
"""

# CODATA 2018
HARTREE_EV = 27.211386245988      # 1 Ha in eV
RYDBERG_EV = HARTREE_EV / 2.0     # 1 Ry in eV
BOHR_ANG = 0.529177210903         # 1 Bohr in Angstrom
ANG_BOHR = 1.0 / BOHR_ANG

MEV = 1.0e-3                      # 1 meV in eV
