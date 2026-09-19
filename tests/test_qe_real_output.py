"""The QE reader against a file pw.x actually wrote.

Independence rule: the reference is the eigenvalue table QE *printed*, parsed
here by this test's own regex.  The text and the XML are written by different
code paths inside QE, so agreeing with one says nothing about having read the
other correctly by accident.

`tests/test_qe.py` covers the same reader against a synthetic XML built in the
test.  That checks the parsing logic; this checks that the logic was aimed at
the right schema.  Both are needed: the synthetic test passed for a reader that
real pw.x would never have fed, because the input it generated was rejected
outright -- `celldm(1)` together with `CELL_PARAMETERS angstrom` makes QE stop
with "lattice parameter specified twice".  Only running it found that.
"""

import re
from pathlib import Path

import numpy as np
import pytest

from altharm.qe.reader import read_bands

DATA = Path(__file__).parent / "data" / "qe_fe_bands"
XML = DATA / "data-file-schema.xml"
OUT = DATA / "bands.out"

ALAT_BOHR = 4.696901          # celldm(1) in bands.in
RADIUS_TPIBA = 0.40           # the circle the k-list was built on
NK, NBND = 24, 16
TEXT_PRECISION = 1e-4         # pw.x prints eigenvalues to four decimals


def text_eigenvalues():
    """Parse the printed bands.  Nothing from altharm is used here."""
    txt = OUT.read_text()
    parts = re.split(r"------ SPIN (UP|DOWN) ------*", txt)
    out = {}
    for tag, body in zip(parts[1::2], parts[2::2]):
        rows = [
            [float(v) for v in m.group(2).split()]
            for m in re.finditer(
                r"k =([-\s\d.]+?)\(\s*\d+ PWs\)\s*bands \(ev\):\s*\n\n"
                r"((?:\s*[-\d.]+\n?)+)",
                body,
            )
        ]
        out[tag] = np.array(rows)
    return out


def text_kpoints():
    """QE prints adjacent negative components with no space between them, so
    split on the numbers rather than on whitespace."""
    heads = re.findall(r"k =([-\s\d.]+?)\(\s*\d+ PWs\)", OUT.read_text())
    return np.array([[float(v) for v in re.findall(r"-?\d+\.\d+", h)] for h in heads])[:NK]


@pytest.fixture(scope="module")
def bands():
    return read_bands(XML)


def test_shapes_and_spin_channels(bands):
    assert bands.lsda is True
    assert bands.nbnd_up == NBND and bands.nbnd_dw == NBND
    assert bands.energies_ev.shape == (2, NK, NBND)


def test_alat_and_fermi_match_the_input_and_the_scf(bands):
    assert np.isclose(bands.alat_bohr, ALAT_BOHR, atol=1e-6)
    assert np.isclose(bands.fermi_ev, 18.1457, atol=5e-5)


def test_eigenvalues_match_what_qe_printed(bands):
    """Hartree in the XML, eV in the text: this is where a wrong conversion shows."""
    ref = text_eigenvalues()
    for i, tag in enumerate(("UP", "DOWN")):
        assert ref[tag].shape == (NK, NBND)
        dev = np.abs(bands.energies_ev[i] - ref[tag]).max()
        assert dev < TEXT_PRECISION, f"spin {tag}: max deviation {dev:.2e} eV"


def test_kpoints_match_what_qe_printed(bands):
    dev = np.abs(bands.kpoints_tpiba - text_kpoints()).max()
    assert dev < TEXT_PRECISION, f"max deviation {dev:.2e} (2 pi / alat)"


def test_the_k_list_came_back_whole_and_in_order(bands):
    """nosym and noinv did their job.

    Without them QE folds an explicit k-list into the irreducible wedge and
    returns fewer points, reordered.  The angular fit would then be reading a
    scrambled circle.  Here all 24 came back, each at exactly the requested
    radius, in the order they were given.
    """
    assert len(bands.kpoints_tpiba) == NK
    radii = np.linalg.norm(bands.kpoints_tpiba, axis=1)
    assert np.allclose(radii, RADIUS_TPIBA, atol=1e-6), (
        f"radii span {radii.min():.6f}..{radii.max():.6f}"
    )
    angles = np.arctan2(bands.kpoints_tpiba[:, 1], bands.kpoints_tpiba[:, 0])
    steps = np.diff(np.unwrap(angles))
    assert np.allclose(steps, 2 * np.pi / NK, atol=1e-5), "the k-points were reordered"


def test_generated_input_does_not_specify_the_lattice_parameter_twice():
    """The bug real pw.x caught: celldm(1) with CELL_PARAMETERS angstrom.

    Regression on the generated file itself, so that the writer cannot drift
    back into an input QE refuses.
    """
    text = (DATA / "bands.in").read_text()
    assert "celldm(1)" in text
    assert "CELL_PARAMETERS alat" in text
    assert "CELL_PARAMETERS angstrom" not in text
    assert "nosym = .true." in text and "noinv = .true." in text
