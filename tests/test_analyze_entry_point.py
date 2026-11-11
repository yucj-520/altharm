"""The single-call entry point (``altharm.analyze_qe`` and the ``altharm`` command) on stored
Quantum ESPRESSO runs.

Independence rule: the reference values were established earlier in this project by the
step-by-step route (reader -> basis -> symmetry -> projection, see PLAN.md and the README of
that time), before ``analyze_qe`` existed. ``analyze_qe`` has to reproduce them from the XML
file alone, without being told the band, the plane of the circle or the symmetry tolerance.
"""
from pathlib import Path

import pytest

from altharm import analyze_qe
from altharm.__main__ import main

DATA = Path(__file__).parent / "data"
CR2SE2AS = DATA / "qe_cr2seas_sm24" / "data-file-schema.xml"
# Cr2Se2As, band 44, 152-point circle at k_z = 0.25 b3, 24x24x8, degauss 0.0150, U = 3.0 eV (meV)
REFERENCE_MEV = {3: 64.293, 9: 43.447, 15: 3.603, 21: 8.662, 27: 1.876, 33: 2.583}


def test_cr2se2as_reference_coefficients():
    a = analyze_qe(CR2SE2AS)
    assert a.band == 44
    assert len(a.dE) == 152
    assert sorted(a.allowed) == [3, 9, 15, 21, 27, 33]
    amplitudes = a.amplitudes_mev()
    for l, ref in REFERENCE_MEV.items():
        assert amplitudes[l] == pytest.approx(ref, abs=6e-4)
    # every allowed harmonic stands far above the forbidden ones
    assert a.noise.floor / 1e-3 < 0.05
    assert min(amplitudes.values()) > 50 * a.noise.floor / 1e-3


def test_cof2_leading_harmonic_is_l2():
    a = analyze_qe(DATA / "cof2_u4_kz025" / "data-file-schema.xml")
    amplitudes = a.amplitudes_mev()
    assert sorted(a.allowed)[:3] == [2, 6, 10]            # d-wave: l = 2 (mod 4)
    assert max(amplitudes, key=amplitudes.get) == 2


def test_mnte_leading_harmonic_is_l3():
    a = analyze_qe(DATA / "mnte_kz025_u4" / "data-file-schema.xml")
    amplitudes = a.amplitudes_mev()
    assert sorted(a.allowed)[:2] == [3, 9]                # g-wave: l = 3 on the circle
    assert max(amplitudes, key=amplitudes.get) == 3


def test_explicit_band_matches_fermi_choice():
    assert analyze_qe(CR2SE2AS, band=44).amplitudes_mev() == analyze_qe(CR2SE2AS).amplitudes_mev()


def test_command_line(capsys):
    assert main([str(CR2SE2AS)]) == 0
    out = capsys.readouterr().out
    assert "allowed l" in out and "noise floor" in out
    assert "64.29" in out
