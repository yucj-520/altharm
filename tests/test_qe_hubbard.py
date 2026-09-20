"""The HUBBARD card, and checking it survived into the run that used it.

Independence rule: the card's text is checked against the syntax QE documents
and against the group's own working `scf.in`, written by hand before this code
existed; the "no U applied" case is checked against a real pw.x XML from a run
that genuinely had no U.
"""

from pathlib import Path

import numpy as np
import pytest

from altharm.kpath import circle
from altharm.qe.writer import bands_input
from altharm.qe.reader import hubbard_is_active, require_hubbard

FE_XML = Path(__file__).parent / "data" / "qe_fe_bands" / "data-file-schema.xml"

BASE = dict(
    prefix="AB", outdir="./tmp", pseudo_dir="/pseudo",
    alat_bohr=7.040816482142,
    cell_ang=np.array([[3.7258, 0.0, 0.0], [-1.8629, 3.2267, 0.0], [0.0, 0.0, 16.918]]),
    species=[("CrU", 51.9961, "cr.UPF"), ("CrD", 51.9961, "cr.UPF")],
    ecutwfc=60.0, nspin=2, starting_magnetization=[1.0, -1.0],
    kcart_tpiba=circle(0.3, 12),
)


def make(**kw):
    args = dict(BASE)
    args.setdefault("positions_crystal",
                    [("CrU", 1 / 3, 2 / 3, 0.59), ("CrD", 1 / 3, 2 / 3, 0.41)])
    args.setdefault("positions_ang", None)
    args.update(kw)
    return bands_input(**args)


def test_hubbard_card_matches_the_documented_syntax():
    """`HUBBARD (projector)` then one `U <label> <value>` per manifold.

    The reference is the group's own working scf.in, which reads
        HUBBARD (ortho-atomic)
          U CrU-3d 3.0
          U CrD-3d 3.0
    """
    text = make(hubbard=[("CrU-3d", 3.0), ("CrD-3d", 3.0)],
                hubbard_projectors="ortho-atomic")
    assert "HUBBARD (ortho-atomic)" in text
    lines = [l.strip() for l in text.splitlines()]
    assert "U CrU-3d 3.000000" in lines
    assert "U CrD-3d 3.000000" in lines
    # the card must come after the namelists, not inside one
    assert lines.index("HUBBARD (ortho-atomic)") > lines.index("/")


def test_no_hubbard_card_when_none_asked_for():
    assert "HUBBARD" not in make()


def test_extra_system_entries_are_written():
    text = make(extra_system={"input_dft": "PBE", "vdw_corr": "grimme-d3",
                              "lspinorb": False})
    assert "input_dft = 'PBE'" in text
    assert "vdw_corr = 'grimme-d3'" in text
    assert "lspinorb = .false." in text


def test_extra_system_refuses_to_duplicate_an_entry():
    """QE takes the last occurrence, so a duplicate would silently win."""
    with pytest.raises(ValueError, match="duplicates an entry"):
        make(extra_system={"nspin": 1})
    with pytest.raises(ValueError, match="duplicates an entry"):
        make(extra_system={"ecutwfc": 30.0})


def test_positions_must_be_given_exactly_once():
    with pytest.raises(ValueError, match="exactly one of"):
        make(positions_ang=None, positions_crystal=None)
    with pytest.raises(ValueError, match="exactly one of"):
        make(positions_ang=[("CrU", 0.0, 0.0, 0.0)],
             positions_crystal=[("CrU", 0.0, 0.0, 0.0)])


def test_crystal_positions_keep_their_exact_thirds():
    """1/3 written to twelve places, not rounded through an Angstrom detour.

    Rounding a third is how a near-degenerate symmetry operation gets lost, and
    losing the sublattice-exchanging one is what turns an altermagnet into
    something this package would describe as a ferromagnet.
    """
    text = make(positions_crystal=[("CrU", 1 / 3, 2 / 3, 0.59),
                                   ("CrD", 1 / 3, 2 / 3, 0.41)])
    assert "ATOMIC_POSITIONS crystal" in text
    assert "0.333333333333" in text and "0.666666666667" in text


# ------------------------------------------------- the inheritance trap

def test_a_real_run_without_u_reports_no_u():
    """The bcc-iron run had no Hubbard U; the check must say so."""
    active, found = hubbard_is_active(FE_XML)
    assert active is False
    assert found == {}


def test_require_hubbard_refuses_and_says_why():
    """The message has to name the failure mode, not just fail.

    A bands run inherits the density but not the HUBBARD card, so the silent
    outcome is plain PBE on a PBE+U density -- no warning from QE at all.
    """
    with pytest.raises(ValueError, match="no Hubbard U applied"):
        require_hubbard(FE_XML)
    with pytest.raises(ValueError, match="silently fell back"):
        require_hubbard(FE_XML)


# ---------------------------------------------- the positive path, on real data

CR2SEAS_XML = Path(__file__).parent / "data" / "qe_cr2seas_p3" / "data-file-schema.xml"


@pytest.mark.skipif(not CR2SEAS_XML.exists(), reason="Cr2SeAs fixture not present")
def test_a_real_pbe_u_run_reports_its_u_in_eV_not_hartree():
    """The regression the bug fix never had.

    QE's XML stores Hubbard_U in **Hartree**: a card asking for 3.0 eV comes
    back as 1.102479665269648e-01. Reading that number without converting
    reports a U about 27 times too small and fails every comparison against
    the eV value the input asked for. Until now this conversion was only ever
    checked by a throwaway script -- the committed tests exercised the
    no-Hubbard case alone, so reverting the fix would have gone unnoticed.
    """
    active, found = hubbard_is_active(CR2SEAS_XML)
    assert active is True
    for key in ("CrU-3d", "CrD-3d"):
        assert key in found, f"expected {key} in {sorted(found)}"
        assert found[key] == pytest.approx(3.0, abs=1e-6), (
            f"{key} = {found[key]} eV -- if this is ~0.11, the Hartree->eV "
            "conversion was lost"
        )


@pytest.mark.skipif(not CR2SEAS_XML.exists(), reason="Cr2SeAs fixture not present")
def test_per_species_entries_do_not_collapse_onto_one_key():
    """The second half of the same bug.

    QE writes one Hubbard_U element per species, and both carry label="3d"
    (specie="CrU" and specie="CrD"). Keying on `label` alone silently drops one
    of them; it went unnoticed on this run only because the two U values happen
    to be equal, so the collision was invisible in the numbers.
    """
    _, found = hubbard_is_active(CR2SEAS_XML)
    assert len(found) == 2, f"expected one entry per species, got {sorted(found)}"
    assert set(found) == {"CrU-3d", "CrD-3d"}


@pytest.mark.skipif(not CR2SEAS_XML.exists(), reason="Cr2SeAs fixture not present")
def test_require_hubbard_accepts_the_matching_expectation():
    found = require_hubbard(CR2SEAS_XML, expected=[("CrU-3d", 3.0), ("CrD-3d", 3.0)])
    assert set(found) == {"CrU-3d", "CrD-3d"}


@pytest.mark.skipif(not CR2SEAS_XML.exists(), reason="Cr2SeAs fixture not present")
def test_require_hubbard_refuses_a_wrong_value_and_a_missing_manifold():
    """Both `expected=` failure branches, which nothing exercised before."""
    with pytest.raises(ValueError, match="expected 4.0 eV"):
        require_hubbard(CR2SEAS_XML, expected=[("CrU-3d", 4.0)])
    with pytest.raises(ValueError, match="expected a Hubbard U on 'Fe-3d'"):
        require_hubbard(CR2SEAS_XML, expected=[("Fe-3d", 3.0)])
