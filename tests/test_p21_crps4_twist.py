"""P21: a 90-degree twist turns CrPS4 from nodal into a d-wave altermagnet.

**Symmetry only.**  Nothing here was run through DFT and nothing here says how
LARGE any splitting is.  It is in the suite because the symmetry question
decides whether a DFT run is worth its cost, and for CrPS4 it is decisive: the
cheapest commensurate twisted bilayer is 288 atoms, roughly 100-400 hours for
one bands run calibrated on this project's own 10-atom timings.

**The monolayer is nodal in-plane, in every AFM configuration tried.**  All four
Cr sit at z = 0, so the lattice has M_z; under sublattice exchange it carries
chi = -1, which forces dE(k_x, k_y, 0) = -dE(k_x, k_y, 0) = 0.

    configuration                         ops  chi=-1   allowed l at k_z=0
    FM (all up)                             8     0     0,2,4,6,8,10,12
    AFM Cr0,2 up / Cr1,3 down              12     6     empty
    AFM Cr0,1 up / Cr2,3 down               8     4     empty
    AFM Cr0,3 up / Cr1,2 down               8     4     empty

**The twist removes exactly that mirror.**  a/b = 1.4847 so 2a ~ 3b (1.03%
mismatch); the smallest square supercell is 21.890 A with each layer strained
+-0.5% to fit -- an idealisation, a real twisted bilayer being a moire.

    bilayer                                ops  chi=-1   allowed l at k_z=0
    untwisted, same configuration           52    26     empty
    90 deg twist, same configuration        --    --     NO exchanging op
    90 deg twist, layer 2 spins reversed     4     2     2,4,6,8,10,12

So the twist opens the in-plane splitting with leading l = 2 -- d-wave.  The
spin reversal is not optional: with both layers alike the twisted structure has
no operation of the non-magnetic parent that exchanges the sublattices, so it
is not an altermagnet at all.  That is a design rule.

Robust against symprec from 1e-4 to 1e-1 -- 4 operations, 2 with chi = -1, same
allowed set throughout.  No cliff, unlike Cr2Se2As.

**Unverified assumption**: which intralayer order is the ground state.  CrPS4 is
usually described with ferromagnetic layers coupled antiferromagnetically, in
which case the monolayer is a ferromagnet and these configurations are not its
ground state.  Everything here is conditional on the configuration existing.
"""

import numpy as np
import pytest

from altharm import CircularBasis, allowed_l, spin_operations

A0, B0, C0 = 10.88896, 7.334202, 6.66998626
BETA = np.deg2rad(91.95294845)
Z = {"Cr": 24, "P": 15, "S": 16}
LABELS = ["Cr"] * 4 + ["P"] * 4 + ["S"] * 16

FRAC = np.array([
    [0.5, 0.75467, 0], [0.5, 0.24533, 0], [0.0, 0.25467, 0], [0.0, 0.74533, 0],
    [0.297353, 0, 0.155247], [0.702647, 0, 0.844753],
    [0.797353, 0.5, 0.155247], [0.202647, 0.5, 0.844753],
    [0.6333465, 0.7672895, 0.721824], [0.3666535, 0.7672895, 0.278176],
    [0.3666535, 0.2327105, 0.278176], [0.6333465, 0.2327105, 0.721824],
    [0.104926, 0, 0.176922], [0.895074, 0, 0.823078],
    [0.629053, 0, 0.131548], [0.370947, 0, 0.868452],
    [0.1333465, 0.2672895, 0.721824], [0.8666535, 0.2672895, 0.278176],
    [0.8666535, 0.7327105, 0.278176], [0.1333465, 0.7327105, 0.721824],
    [0.604926, 0.5, 0.176922], [0.395074, 0.5, 0.823078],
    [0.129053, 0.5, 0.131548], [0.870947, 0.5, 0.868452]])

AFM_A = [1, -1, 1, -1]          # Cr0,2 up / Cr1,3 down -- the y -> -y pairing
INPLANE = np.stack([np.cos(t := 2 * np.pi * np.arange(360) / 360),
                    np.sin(t), np.zeros(360)], axis=1)


def _monolayer(spins):
    cell = np.array([[A0, 0, 0], [0, B0, 0],
                     [C0 * np.cos(BETA), 0, C0 * np.sin(BETA)]])
    nums = [Z[s] for s in LABELS]
    mag = list(spins) + [0.0] * 20
    return cell, FRAC, nums, mag


def _bilayer(twist_deg, spins1, spins2, gap=6.6, c_tot=24.0):
    """Two layers in the smallest commensurate square cell.

    Each layer is strained anisotropically so 2a and 3b both equal the square
    side exactly; without that the rotated layer does not land on the supercell
    lattice and spglib finds no symmetry at all.
    """
    side = 0.5 * (2 * A0 + 3 * B0)
    a1, b1 = side / 2.0, side / 3.0
    zc = (FRAC[:, 2] + 0.5) % 1.0 - 0.5
    pos, mag = [], []
    for i in range(2):                               # layer 1: a along x, 2 x 3
        for j in range(3):
            for n in range(24):
                pos.append([(FRAC[n, 0] + i) * a1, (FRAC[n, 1] + j) * b1,
                            zc[n] * C0 - gap / 2])
                mag.append(spins1[n] if n < 4 else 0.0)
    th = np.deg2rad(twist_deg)
    rot = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    ninety = abs(twist_deg - 90.0) < 1e-9
    nx, ny = (3, 2) if ninety else (2, 3)
    for i in range(nx):
        for j in range(ny):
            for n in range(24):
                xy = rot @ np.array([FRAC[n, 0] * a1, FRAC[n, 1] * b1])
                pos.append([xy[0] + i * (b1 if ninety else a1),
                            xy[1] + j * (a1 if ninety else b1),
                            zc[n] * C0 + gap / 2])
                mag.append(spins2[n] if n < 4 else 0.0)
    pos = np.array(pos)
    cell = np.array([[side, 0, 0], [0, side, 0], [0, 0, c_tot]])
    nums = [Z[s] for s in LABELS] * 12
    return cell, (pos @ np.linalg.inv(cell)) % 1.0, nums, mag


def _allowed_inplane(args, symprec, lmax=12):
    cell, frac, nums, mag = args
    ops = spin_operations(cell, frac, nums, mag, symprec=symprec)
    return sorted(allowed_l(CircularBasis(lmax), ops, INPLANE)), ops


def test_the_cr_atoms_all_sit_in_one_plane():
    """Which is why M_z exists and kills the in-plane splitting."""
    assert np.allclose(FRAC[:4, 2], 0.0)


def test_the_monolayer_is_nodal_in_plane_for_every_afm_arrangement():
    for spins in ([1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]):
        al, ops = _allowed_inplane(_monolayer(spins), 1e-3)
        assert al == [], (spins, al)
        assert any(o.chi == -1 for o in ops), spins


def test_the_ferromagnetic_monolayer_is_not_nodal():
    """The contrast: with no sublattice exchange there is nothing to forbid."""
    al, ops = _allowed_inplane(_monolayer([1, 1, 1, 1]), 1e-3)
    assert al != []
    assert all(o.chi == +1 for o in ops)


def test_the_untwisted_bilayer_stays_nodal():
    al, _ = _allowed_inplane(_bilayer(0.0, AFM_A, AFM_A), 5e-2)
    assert al == [], al


def test_the_twist_without_a_spin_flip_is_not_an_altermagnet():
    """No operation of the parent exchanges the sublattices -- it raises."""
    with pytest.raises(Exception):
        _allowed_inplane(_bilayer(90.0, AFM_A, AFM_A), 5e-2)


def test_the_twist_with_a_spin_flip_opens_a_d_wave_splitting():
    al, ops = _allowed_inplane(_bilayer(90.0, AFM_A, [-s for s in AFM_A]), 5e-2)
    assert al == [2, 4, 6, 8, 10, 12], al
    assert al[0] == 2                              # leading harmonic: d-wave
    assert sum(1 for o in ops if o.chi == -1) == 2


@pytest.mark.parametrize("symprec", [1e-4, 1e-3, 1e-2, 1e-1])
def test_the_result_has_no_symprec_cliff(symprec):
    """Cr2Se2As taught this project that the tolerance can decide the answer."""
    al, ops = _allowed_inplane(
        _bilayer(90.0, AFM_A, [-s for s in AFM_A]), symprec)
    assert al == [2, 4, 6, 8, 10, 12], (symprec, al)
    assert len(ops) == 4
