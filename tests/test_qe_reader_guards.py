"""Malformed-XML paths in the QE reader that nothing had ever triggered.

Coverage showed every one of these `raise ValueError` branches unexecuted --
the reader had only ever been handed well-formed files. A guard nobody has
run is a guard nobody knows works; on a silently truncated or
wrong-schema-version file, a reader that crashes on a NoneType three
functions later is far worse than one that says which element is missing.

Independence rule: each fixture is a hand-written XML string broken in one
named way, and the test asserts the message names that way.
"""

import numpy as np
import pytest

from altharm.qe.reader import read_bands, spin_splitting

HEAD = '<?xml version="1.0"?><qes:espresso xmlns:qes="http://www.quantum-espresso.org/ns/qes/qes-1.0">'


def _xml(body, tmp_path, name="data-file-schema.xml"):
    p = tmp_path / name
    p.write_text(HEAD + body + "</qes:espresso>")
    return p


def _structure(alat='alat="10.0"'):
    return (f'<atomic_structure {alat}><cell>'
            '<a1>10.0 0 0</a1><a2>0 10.0 0</a2><a3>0 0 10.0</a3>'
            '</cell></atomic_structure>')


def _bands(nks=1, lsda="false", nbnd=2, ks=None, fermi="<fermi_energy>0.2</fermi_energy>"):
    ks = ks if ks is not None else (
        '<ks_energies><k_point>0.0 0.0 0.0</k_point>'
        '<eigenvalues>0.1 0.2</eigenvalues></ks_energies>'
    )
    return (f'<band_structure><lsda>{lsda}</lsda><nbnd>{nbnd}</nbnd>'
            f'<nks>{nks}</nks>{fermi}{ks}</band_structure>')


def test_missing_output_element_is_named(tmp_path):
    p = _xml("<input/>", tmp_path)
    with pytest.raises(ValueError, match="output"):
        read_bands(p)


def test_missing_alat_attribute_is_named(tmp_path):
    p = _xml("<output>" + _structure(alat="") + _bands() + "</output>", tmp_path)
    with pytest.raises(ValueError, match="no 'alat' attribute"):
        read_bands(p)


def test_nks_disagreeing_with_the_ks_energies_count_is_refused(tmp_path):
    p = _xml("<output>" + _structure() + _bands(nks=3) + "</output>", tmp_path)
    with pytest.raises(ValueError, match="says 3 but 1"):
        read_bands(p)


def test_zero_kpoints_is_refused(tmp_path):
    body = ('<band_structure><lsda>false</lsda><nbnd>2</nbnd><nks>0</nks>'
            '<fermi_energy>0.2</fermi_energy></band_structure>')
    p = _xml("<output>" + _structure() + body + "</output>", tmp_path)
    with pytest.raises(ValueError, match="zero <ks_energies>"):
        read_bands(p)


def test_a_k_point_with_the_wrong_number_of_components_is_refused(tmp_path):
    ks = ('<ks_energies><k_point>0.0 0.0</k_point>'
          '<eigenvalues>0.1 0.2</eigenvalues></ks_energies>')
    p = _xml("<output>" + _structure() + _bands(ks=ks) + "</output>", tmp_path)
    with pytest.raises(ValueError, match="has 2 components, expected 3"):
        read_bands(p)


def test_an_eigenvalue_count_that_contradicts_nbnd_is_refused(tmp_path):
    ks = ('<ks_energies><k_point>0.0 0.0 0.0</k_point>'
          '<eigenvalues>0.1 0.2 0.3</eigenvalues></ks_energies>')
    p = _xml("<output>" + _structure() + _bands(ks=ks) + "</output>", tmp_path)
    with pytest.raises(ValueError, match="has 3 values, expected 2"):
        read_bands(p)


def test_two_fermi_energies_are_averaged_rather_than_refused(tmp_path):
    """A fixed-spin-moment metallic run reports two Fermi levels. The reader
    averages them -- checked here against the hand-computed mean in eV."""
    from altharm.units import HARTREE_EV
    fermi = "<two_fermi_energies>0.10 0.30</two_fermi_energies>"
    p = _xml("<output>" + _structure() + _bands(fermi=fermi) + "</output>", tmp_path)
    qb = read_bands(p)
    assert qb.fermi_ev == pytest.approx(0.20 * HARTREE_EV, rel=1e-12)


def test_neither_fermi_element_present_is_refused(tmp_path):
    p = _xml("<output>" + _structure() + _bands(fermi="") + "</output>", tmp_path)
    with pytest.raises(ValueError, match="neither <fermi_energy> nor"):
        read_bands(p)


def test_spin_splitting_refuses_a_non_spin_polarised_run(tmp_path):
    p = _xml("<output>" + _structure() + _bands() + "</output>", tmp_path)
    qb = read_bands(p)
    assert qb.lsda is False
    with pytest.raises(ValueError):
        spin_splitting(qb, 0)


def test_hubbard_key_falls_back_when_specie_and_label_are_not_both_present(tmp_path):
    """Older QE schemas label Hubbard_U differently -- by specie alone, by an
    index, or by neither. The fallback key chain had never been taken."""
    from altharm.qe.reader import hubbard_is_active
    from altharm.units import HARTREE_EV
    body = ('<output><dftU>'
            '<Hubbard_U specie="Ni">0.1102479665269648</Hubbard_U>'
            '<Hubbard_U index="2">0.2204959330539296</Hubbard_U>'
            '</dftU></output>')
    active, found = hubbard_is_active(_xml(body, tmp_path))
    assert active is True
    assert found["Ni"] == pytest.approx(0.1102479665269648 * HARTREE_EV)
    assert found["2"] == pytest.approx(0.2204959330539296 * HARTREE_EV)


def test_an_unparseable_hubbard_value_is_skipped_not_fatal(tmp_path):
    """A malformed or empty Hubbard_U element must be stepped over, so one bad
    entry cannot take down a reader that is only being asked whether U was on."""
    from altharm.qe.reader import hubbard_is_active
    body = ('<output><dftU>'
            '<Hubbard_U specie="Cr" label="3d">not-a-number</Hubbard_U>'
            '<Hubbard_U specie="Fe" label="3d"></Hubbard_U>'
            '<Hubbard_U specie="Ni" label="3d">0.1102479665269648</Hubbard_U>'
            '</dftU></output>')
    active, found = hubbard_is_active(_xml(body, tmp_path))
    assert active is True
    assert set(found) == {"Ni-3d"}, f"only the parseable entry should survive: {found}"


def test_a_zero_hubbard_u_does_not_count_as_active(tmp_path):
    """U = 0 is the same as no U; it must not make the check pass."""
    from altharm.qe.reader import hubbard_is_active
    body = '<output><dftU><Hubbard_U specie="Cr" label="3d">0.0</Hubbard_U></dftU></output>'
    active, found = hubbard_is_active(_xml(body, tmp_path))
    assert active is False and found == {}
