"""Tests for altharm.qe.

Rule for this package: a test must check against something independent of
the code it exercises, not just re-run that code and compare to itself.
Where possible below, expected numbers are worked out by hand with explicit
literals (the CODATA Bohr radius, a hand-multiplied Hartree->eV factor) or
parsed with plain string operations, never by calling the function under
test to produce its own reference value.
"""

from __future__ import annotations

import math
import textwrap

import numpy as np
import pytest

from altharm.qe.writer import cart_to_tpiba, tpiba_to_cart, kpoints_card, bands_input
from altharm.qe.reader import read_bands, spin_splitting


# CODATA 2018, same source as altharm/units.py -- duplicated here as a
# literal, deliberately, so this test does not depend on units.py agreeing
# with itself.
_BOHR_ANG = 0.529177210903
_HARTREE_EV = 27.211386245988


# ---------------------------------------------------------------------------
# cart_to_tpiba / tpiba_to_cart
# ---------------------------------------------------------------------------


def test_tpiba_roundtrip():
    rng = np.random.default_rng(0)
    k = rng.uniform(-2.0, 2.0, size=(7, 3))
    alat_bohr = 12.345
    back = tpiba_to_cart(cart_to_tpiba(k, alat_bohr), alat_bohr)
    np.testing.assert_allclose(back, k, atol=1e-12, rtol=1e-12)


def test_cart_to_tpiba_hand_value():
    # Hand-worked reference, independent of cart_to_tpiba's own arithmetic:
    # alat = 10 Bohr -> alat in Angstrom -> 2*pi/alat in 1/Angstrom -> a
    # k of 1/Angstrom along x divided by that factor.
    alat_bohr = 10.0
    alat_ang = alat_bohr * _BOHR_ANG  # 5.29177210903 Angstrom
    two_pi_over_alat_invang = 2.0 * math.pi / alat_ang  # ~1.187349942 1/Ang
    expected_kx_tpiba = 1.0 / two_pi_over_alat_invang  # ~0.842211689

    k_invang = np.array([1.0, 0.0, 0.0])
    got = cart_to_tpiba(k_invang, alat_bohr)

    assert got[1] == 0.0
    assert got[2] == 0.0
    assert got[0] == pytest.approx(expected_kx_tpiba, rel=1e-12)
    # Pin the literal too, so a future accidental change to _BOHR_ANG above
    # can't quietly drag the "hand" value along with the code under test.
    assert got[0] == pytest.approx(0.8422116888679485, rel=1e-9)


def test_tpiba_to_cart_hand_value():
    # Same hand-worked factor, inverse direction.
    alat_bohr = 10.0
    alat_ang = alat_bohr * _BOHR_ANG
    two_pi_over_alat_invang = 2.0 * math.pi / alat_ang
    k_tpiba = np.array([0.5, -0.25, 0.0])
    expected = k_tpiba * two_pi_over_alat_invang
    got = tpiba_to_cart(k_tpiba, alat_bohr)
    np.testing.assert_allclose(got, expected, rtol=1e-12)


# ---------------------------------------------------------------------------
# kpoints_card
# ---------------------------------------------------------------------------


def test_kpoints_card_parses_back_with_plain_string_splitting():
    kcart = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.1, 0.2, 0.3],
            [-0.4, 0.5, -0.6],
        ]
    )
    weights = np.array([1.0, 2.5, 0.25])
    text = kpoints_card(kcart, weights=weights)

    lines = text.strip("\n").split("\n")
    assert lines[0].split() == ["K_POINTS", "tpiba"]
    nk = int(lines[1].strip())
    assert nk == 3
    assert len(lines) == 2 + nk

    for i, line in enumerate(lines[2:]):
        fields = [float(x) for x in line.split()]
        assert len(fields) == 4
        assert fields[:3] == pytest.approx(list(kcart[i]))
        assert fields[3] == pytest.approx(weights[i])


def test_kpoints_card_default_weights_are_one():
    kcart = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    text = kpoints_card(kcart)
    lines = text.strip("\n").split("\n")
    for line in lines[2:]:
        fields = [float(x) for x in line.split()]
        assert fields[3] == pytest.approx(1.0)


def test_kpoints_card_rejects_bad_weight_shape():
    kcart = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    with pytest.raises(ValueError):
        kpoints_card(kcart, weights=[1.0, 2.0, 3.0])


# ---------------------------------------------------------------------------
# bands_input
# ---------------------------------------------------------------------------


def _minimal_bands_input(**overrides):
    kwargs = dict(
        prefix="test",
        outdir="./out",
        pseudo_dir="./pseudo",
        alat_bohr=10.0,
        cell_ang=np.eye(3) * 5.0,
        species=[("Fe", 55.845, "Fe.upf")],
        positions_ang=[("Fe", 0.0, 0.0, 0.0)],
        ecutwfc=60.0,
        kcart_tpiba=np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]]),
    )
    kwargs.update(overrides)
    return bands_input(**kwargs)


def test_bands_input_sets_nosym_and_noinv():
    text = _minimal_bands_input()
    system_block = text.split("&SYSTEM")[1].split("/")[0]
    assert "nosym = .true." in system_block
    assert "noinv = .true." in system_block


def test_bands_input_exposes_convergence_controls():
    text = _minimal_bands_input(diago_thr_init=1e-6, conv_thr=1e-12)
    electrons_block = text.split("&ELECTRONS")[1].split("/")[0]
    assert "1.000e-06" in electrons_block or "1e-06" in electrons_block.lower()
    assert "diago_thr_init" in electrons_block
    assert "conv_thr" in electrons_block


def test_bands_input_calculation_is_bands():
    text = _minimal_bands_input()
    assert "calculation = 'bands'" in text


# ---------------------------------------------------------------------------
# reader: synthetic data-file-schema.xml
# ---------------------------------------------------------------------------

_NON_LSDA_XML = """<?xml version="1.0" encoding="UTF-8"?>
<qes:espresso xmlns:qes="http://www.quantum-espresso.org/ns/qes/qes-1.0">
  <output>
    <atomic_structure alat="8.0" bravais_index="0">
      <cell>
        <a1>8.0 0.0 0.0</a1>
        <a2>0.0 8.0 0.0</a2>
        <a3>0.0 0.0 8.0</a3>
      </cell>
    </atomic_structure>
    <band_structure>
      <lsda>false</lsda>
      <nbnd>4</nbnd>
      <nks>2</nks>
      <fermi_energy>0.20</fermi_energy>
      <ks_energies>
        <k_point weight="1.0">0.0 0.0 0.0</k_point>
        <eigenvalues size="4">0.10 0.20 0.30 0.40</eigenvalues>
      </ks_energies>
      <ks_energies>
        <k_point weight="1.0">0.5 0.0 0.0</k_point>
        <eigenvalues size="4">0.11 0.21 0.31 0.41</eigenvalues>
      </ks_energies>
    </band_structure>
  </output>
</qes:espresso>
"""

# LSDA case with nbnd_up != nbnd_dw, the case the task explicitly calls out
# as easy to get wrong (must split, not interleave; must not assume equal
# lengths).
_LSDA_XML = """<?xml version="1.0" encoding="UTF-8"?>
<qes:espresso xmlns:qes="http://www.quantum-espresso.org/ns/qes/qes-1.0">
  <output>
    <atomic_structure alat="10.0" bravais_index="0">
      <cell>
        <a1>10.0 0.0 0.0</a1>
        <a2>0.0 10.0 0.0</a2>
        <a3>0.0 0.0 10.0</a3>
      </cell>
    </atomic_structure>
    <band_structure>
      <lsda>true</lsda>
      <nbnd_up>3</nbnd_up>
      <nbnd_dw>2</nbnd_dw>
      <nks>2</nks>
      <fermi_energy>0.5</fermi_energy>
      <ks_energies>
        <k_point weight="1.0">0.1 0.2 0.3</k_point>
        <eigenvalues size="5">0.01 0.02 0.03 0.04 0.05</eigenvalues>
      </ks_energies>
      <ks_energies>
        <k_point weight="1.0">0.4 0.5 0.6</k_point>
        <eigenvalues size="5">0.06 0.07 0.08 0.09 0.10</eigenvalues>
      </ks_energies>
    </band_structure>
  </output>
</qes:espresso>
"""


def test_read_bands_non_lsda(tmp_path):
    p = tmp_path / "data-file-schema.xml"
    p.write_text(_NON_LSDA_XML)

    qb = read_bands(p)

    assert qb.lsda is False
    assert qb.nbnd_up == 4
    assert qb.nbnd_dw == 4
    assert qb.energies_ev.shape == (1, 2, 4)
    np.testing.assert_allclose(qb.alat_bohr, 8.0)
    np.testing.assert_allclose(qb.cell_ang, np.eye(3) * 8.0 * _BOHR_ANG)
    np.testing.assert_allclose(qb.kpoints_tpiba[0], [0.0, 0.0, 0.0])
    np.testing.assert_allclose(qb.kpoints_tpiba[1], [0.5, 0.0, 0.0])

    expected_k0 = np.array([0.10, 0.20, 0.30, 0.40]) * _HARTREE_EV
    expected_k1 = np.array([0.11, 0.21, 0.31, 0.41]) * _HARTREE_EV
    np.testing.assert_allclose(qb.energies_ev[0, 0], expected_k0, rtol=1e-12)
    np.testing.assert_allclose(qb.energies_ev[0, 1], expected_k1, rtol=1e-12)
    assert qb.fermi_ev == pytest.approx(0.20 * _HARTREE_EV, rel=1e-12)


def test_read_bands_lsda_different_nbnd(tmp_path):
    p = tmp_path / "data-file-schema.xml"
    p.write_text(_LSDA_XML)

    qb = read_bands(p)

    assert qb.lsda is True
    assert qb.nbnd_up == 3
    assert qb.nbnd_dw == 2
    # nbnd dimension is padded up to max(nbnd_up, nbnd_dw) = 3.
    assert qb.energies_ev.shape == (2, 2, 3)

    # Up channel: first nbnd_up=3 values of each eigenvalues block.
    expected_up_k0 = np.array([0.01, 0.02, 0.03]) * _HARTREE_EV
    expected_up_k1 = np.array([0.06, 0.07, 0.08]) * _HARTREE_EV
    np.testing.assert_allclose(qb.energies_ev[0, 0], expected_up_k0, rtol=1e-12)
    np.testing.assert_allclose(qb.energies_ev[0, 1], expected_up_k1, rtol=1e-12)

    # Down channel: remaining nbnd_dw=2 values, in the last slot NaN-padded.
    expected_dw_k0 = np.array([0.04, 0.05]) * _HARTREE_EV
    expected_dw_k1 = np.array([0.09, 0.10]) * _HARTREE_EV
    np.testing.assert_allclose(qb.energies_ev[1, 0, :2], expected_dw_k0, rtol=1e-12)
    np.testing.assert_allclose(qb.energies_ev[1, 1, :2], expected_dw_k1, rtol=1e-12)
    assert np.isnan(qb.energies_ev[1, 0, 2])
    assert np.isnan(qb.energies_ev[1, 1, 2])

    np.testing.assert_allclose(qb.kpoints_tpiba[0], [0.1, 0.2, 0.3])
    np.testing.assert_allclose(qb.kpoints_tpiba[1], [0.4, 0.5, 0.6])
    assert qb.fermi_ev == pytest.approx(0.5 * _HARTREE_EV, rel=1e-12)
    np.testing.assert_allclose(qb.cell_ang, np.eye(3) * 10.0 * _BOHR_ANG)


def test_spin_splitting_values(tmp_path):
    p = tmp_path / "data-file-schema.xml"
    p.write_text(_LSDA_XML)
    qb = read_bands(p)

    # Hand-worked from the same literals used to build the XML above:
    # band 0: up (0.01, 0.06), down (0.04, 0.09) Hartree.
    dE0 = spin_splitting(qb, 0)
    expected0 = (np.array([0.01, 0.06]) - np.array([0.04, 0.09])) * _HARTREE_EV
    np.testing.assert_allclose(dE0, expected0, rtol=1e-12)

    # band 1: up (0.02, 0.07), down (0.05, 0.10) Hartree.
    dE1 = spin_splitting(qb, 1)
    expected1 = (np.array([0.02, 0.07]) - np.array([0.05, 0.10])) * _HARTREE_EV
    np.testing.assert_allclose(dE1, expected1, rtol=1e-12)


def test_spin_splitting_requires_lsda(tmp_path):
    p = tmp_path / "data-file-schema.xml"
    p.write_text(_NON_LSDA_XML)
    qb = read_bands(p)
    with pytest.raises(ValueError):
        spin_splitting(qb, 0)


def test_spin_splitting_band_index_out_of_range(tmp_path):
    p = tmp_path / "data-file-schema.xml"
    p.write_text(_LSDA_XML)
    qb = read_bands(p)
    # nbnd_dw = 2, so band index 2 exists for "up" but not for "down".
    with pytest.raises(ValueError):
        spin_splitting(qb, 2)


def test_read_bands_missing_band_structure_raises(tmp_path):
    xml = """<?xml version="1.0" encoding="UTF-8"?>
    <qes:espresso xmlns:qes="http://www.quantum-espresso.org/ns/qes/qes-1.0">
      <output>
        <atomic_structure alat="8.0" bravais_index="0">
          <cell><a1>8 0 0</a1><a2>0 8 0</a2><a3>0 0 8</a3></cell>
        </atomic_structure>
      </output>
    </qes:espresso>
    """
    p = tmp_path / "data-file-schema.xml"
    p.write_text(textwrap.dedent(xml))
    with pytest.raises(ValueError, match="band_structure"):
        read_bands(p)


def test_read_bands_nks_mismatch_raises(tmp_path):
    # nks says 2 but only one ks_energies block is present.
    xml = _NON_LSDA_XML.replace("<nks>2</nks>", "<nks>2</nks>").replace(
        """      <ks_energies>
        <k_point weight="1.0">0.5 0.0 0.0</k_point>
        <eigenvalues size="4">0.11 0.21 0.31 0.41</eigenvalues>
      </ks_energies>
""",
        "",
    )
    p = tmp_path / "data-file-schema.xml"
    p.write_text(xml)
    with pytest.raises(ValueError, match="nks"):
        read_bands(p)
