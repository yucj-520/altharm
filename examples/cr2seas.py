"""Cr2SeAs: the harmonics of the spin splitting, from the group's Wannier model.

This is the first real material the package has been pointed at, and it doubles
as a cross-check: the symmetry module predicts which harmonics may be non-zero
from the structure alone, and the Wannier Hamiltonian -- fitted to DFT, knowing
nothing about that prediction -- supplies the actual numbers.  They have to
agree.

The model is the v03 one, PBE+U (U = 3.0 eV), 38 p-d orbitals per spin channel.
"""

from __future__ import annotations

import sys
import numpy as np

from altharm.basis import CircularBasis
from altharm.kpath import circle
from altharm.project import project
from altharm.noise import noise_report
from altharm.symmetry import (
    spin_operations, allowed_l, symprec_scan, format_symprec_scan,
)
from altharm import wannier

ROOT = "/scratch/work/yuc11/HIGH_THROUGHT_WANNIERcr2X3/1.cr2se2as/Cr2SeAs_wannier_v03"
E_FERMI = 8.1169  # eV, from the model's own evaluate_report.txt

CELL = np.array([
    [3.7258396285, 0.0000002055, 0.0000000080],
    [-1.8629196361, 3.2266718606, -0.0000000204],
    [0.0000000476, -0.0000000907, 16.9179784926],
])

# order as in c.up.win; the README's table fixes which Cr is which sublattice
POSITIONS = np.array([
    [1 / 3, 2 / 3, 0.5896431867],   # Cr1  up
    [1 / 3, 2 / 3, 0.4103569463],   # Cr2  down
    [2 / 3, 1 / 3, 0.0896259946],   # Cr3  down
    [2 / 3, 1 / 3, 0.9103736386],   # Cr4  up
    [0.0, 0.0, 0.3332384847],       # Se1
    [0.0, 0.0, 0.6667617671],       # Se2
    [1 / 3, 2 / 3, 0.8334101415],   # Se3
    [1 / 3, 2 / 3, 0.1665899648],   # Se4
    [2 / 3, 1 / 3, 0.5000002984],   # As1
    [0.0, 0.0, 0.0],                # As2
])
NUMBERS = np.array([24, 24, 24, 24, 34, 34, 34, 34, 33, 33])
MAGMOMS = np.array([+3.0, -3.0, -3.0, +3.0, 0, 0, 0, 0, 0, 0])


def reciprocal(cell):
    return 2.0 * np.pi * np.linalg.inv(np.asarray(cell, float)).T


def main():
    lmax = 12
    n = 4 * lmax + 8
    B = reciprocal(CELL)
    b_inplane = np.linalg.norm(B[0])
    b_z = np.linalg.norm(B[2])

    # symprec matters here: see the scan.  The horizontal mirror that exchanges
    # the two Cr sublattices maps atom onto atom to within 1.0e-5 A, which is
    # spglib's default tolerance, so the default drops the one operation that
    # makes this material an altermagnet.
    print(format_symprec_scan(symprec_scan(CELL, POSITIONS, NUMBERS, MAGMOMS)))
    print()
    ops = spin_operations(CELL, POSITIONS, NUMBERS, MAGMOMS, symprec=1e-4)
    print(f"non-magnetic parent operations kept : {len(ops)}")
    print(f"  sublattice-preserving (chi = +1)  : {sum(o.chi == +1 for o in ops)}")
    print(f"  sublattice-exchanging (chi = -1)  : {sum(o.chi == -1 for o in ops)}")

    basis = CircularBasis(lmax)

    print("\nallowed harmonics, by k_z:")
    for kz_frac in (0.0, 0.1, 0.25, 0.5):
        pts = circle(0.25 * b_inplane, n, center=(0.0, 0.0, kz_frac * b_z))
        got = allowed_l(basis, ops, pts)
        print(f"  k_z = {kz_frac:4.2f} b3 : {sorted(got) if got else 'nothing allowed'}")

    # --- the numbers, from the Wannier model -------------------------------
    print("\nreading the Wannier Hamiltonians ...")
    up = wannier.read_hr(f"{ROOT}/hamiltonian/c.up_hr.dat", cell=CELL)
    dn = wannier.read_hr(f"{ROOT}/hamiltonian/c.dn_hr.dat", cell=CELL)
    print(f"  num_wann = {up.num_wann}, R-vectors = {up.num_rpts}, "
          f"reach = {wannier.hopping_range(up):.2f} A")

    for kz_frac in (0.0, 0.25):
        pts = circle(0.25 * b_inplane, n, center=(0.0, 0.0, kz_frac * b_z))
        allowed = allowed_l(basis, ops, pts)
        frac = wannier.cart_to_frac_k(pts, CELL)
        e_up = wannier.eigenvalues(up, frac)
        e_dn = wannier.eigenvalues(dn, frac)
        band = int(np.argmin(np.abs(e_up.mean(axis=0) - E_FERMI)))
        dE = e_up[:, band] - e_dn[:, band]

        print(f"\n--- k_z = {kz_frac} b3, band {band} "
              f"(mean E - E_F = {e_up[:, band].mean() - E_FERMI:+.3f} eV) ---")
        print(f"  |dE| max = {np.abs(dE).max() / 1e-3:.4f} meV, "
              f"mean = {dE.mean() / 1e-3:+.4f} meV")
        fit = project(basis, pts, dE)
        rep = noise_report(fit, allowed)
        print(f"  allowed l from symmetry: {sorted(allowed) if allowed else '{}'}")
        print(rep.summary(unit_name="meV", unit_scale=1e-3))


if __name__ == "__main__":
    sys.exit(main())
