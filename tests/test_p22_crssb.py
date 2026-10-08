"""P22: Cr2S2Sb -- the multi-harmonic character is structural, and my
explanation for it was wrong.

Job 20820707. Cr2S2Sb is the sister compound of Cr2Se2As: same structure type
(4 Cr + 4 chalcogen + 2 pnictogen, same Wyckoff pattern), different chemistry,
structure taken from the group's own tree. Settings identical to P12's
Cr2Se2As runs except two deliberate improvements from this project's own
findings -- degauss 0.0150 rather than 0.0037 (P3) and N = 151 rather than 152,
coprime to the allowed-set period (P14).

**It is a semiconductor and still multi-harmonic.** No band crosses E_F: bands
53 (-0.56 eV) and 54 (+0.27 eV) straddle it with a 712 meV gap, where Cr2Se2As
is metallic with bands 43/44 inside +-130 meV. The allowed set is identical --
[3, 9, 15, 21, 27, 33, 39, 45] -- as the shared structure requires.

    band  E - E_F    l=3        l=9       next/leading
     53   -0.561   20.415 meV   9.112      44.63%
     54   +0.274    2.537       2.698     106.35%   (l=9 exceeds l=3)

So the multi-harmonic character belongs to the structure type, not to Cr2Se2As.

**And it falsifies the mechanism I had been offering.** Asked why MnTe is nearly
single-harmonic, I proposed that richness tracks how close the band sits to its
neighbours along the circle -- avoided crossings make sharp features in theta,
which is P10's mechanism. Across four materials:

    material (band)      nearest-neighbour gap   next/leading
    Cr2Se2As (44)              20.1 meV             67.58%
    CoF2 (30)                  92.0 meV              0.02%
    Cr2S2Sb (54)               96.1 meV            106.35%
    Cr2S2Sb (53)              112.6 meV             44.63%
    MnTe (20)                 388.3 meV              5.63%

Not monotone; log-log correlation -0.22. CoF2's gap is SMALLER than either
Cr2S2Sb band's and its harmonic content is two thousand times lower.

P10's mechanism stands -- level repulsion generates the harmonics, measured on
Cr2Se2As. What is withdrawn is the claim that neighbour spacing predicts
richness across materials. Four points cannot support it and they contradict it.
What the data group by is structure type.
"""

from pathlib import Path

import numpy as np
import pytest

from altharm import CircularBasis, noise_report, project
from altharm.qe.reader import read_bands, read_structure
from altharm.symmetry import allowed_l, spin_operations

XML = Path(__file__).parent / "data" / "crssb" / "data-file-schema.xml"
EF = 9.7926
LMAX = 45

pytestmark = pytest.mark.skipif(not XML.exists(),
                                reason="the Cr2S2Sb fixture is not present")


def _fit(band):
    q = read_bands(XML)
    st = read_structure(XML)
    k = q.kpoints_tpiba
    th = np.arctan2(k[:, 1], k[:, 0])
    rad = np.linalg.norm(k[0][:2])
    kz = k[0, 2]
    dirs = np.stack([np.cos(th), np.sin(th), (kz / rad) * np.ones_like(th)], axis=1)
    ops = spin_operations(st.cell_ang, st.scaled_positions, st.numbers,
                          st.magmoms, symprec=1e-3)
    allowed = sorted(allowed_l(CircularBasis(LMAX), ops, dirs))
    f = project(CircularBasis(LMAX), dirs,
                q.energies_ev[0, :, band] - q.energies_ev[1, :, band])
    return f.amplitude_by_l(), allowed, noise_report(f, set(allowed)).floor, q


def test_the_sample_count_is_coprime_to_the_period():
    """P14 applied to a new run rather than only recorded."""
    from math import gcd
    q = read_bands(XML)
    assert len(q.kpoints_tpiba) == 151
    assert gcd(151, 6) == 1


def test_the_allowed_set_matches_cr2se2as():
    """Same structure type, so symmetry must give the same answer."""
    _, allowed, _, _ = _fit(53)
    assert allowed == [3, 9, 15, 21, 27, 33, 39, 45]


def test_it_is_a_semiconductor_unlike_cr2se2as():
    q = read_bands(XML)
    e = q.energies_ev[0]
    assert not (e[:, 53].min() < EF < e[:, 53].max())
    assert not (e[:, 54].min() < EF < e[:, 54].max())
    assert float(np.min(e[:, 54] - e[:, 53])) > 0.5        # 712 meV


@pytest.mark.parametrize("band,floor_ratio", [(53, 0.3), (54, 0.8)])
def test_both_bands_are_multi_harmonic(band, floor_ratio):
    amp, allowed, _, _ = _fit(band)
    ratio = amp[allowed[1]] / amp[allowed[0]]
    assert ratio > floor_ratio, (band, ratio)


def test_the_conduction_band_has_l9_above_l3():
    """106% -- the second allowed harmonic exceeds the first."""
    amp, allowed, _, _ = _fit(54)
    assert amp[9] > amp[3]


def test_neighbour_gap_does_not_predict_harmonic_richness():
    """The explanation I had been offering, falsified on its own terms.

    Keep this as a test so the withdrawn claim cannot quietly return.
    """
    q = read_bands(XML)
    e = q.energies_ev[0]
    gap53 = min(float(np.min(e[:, 53] - e[:, 52])),
                float(np.min(e[:, 54] - e[:, 53])))
    amp53, allowed, _, _ = _fit(53)
    ratio53 = amp53[allowed[1]] / amp53[allowed[0]]

    # CoF2, measured in P8: a SMALLER gap with a FAR smaller ratio
    cof2_gap_mev, cof2_ratio = 92.0, 0.00021
    assert gap53 * 1e3 > cof2_gap_mev          # Cr2S2Sb's band is further apart
    assert ratio53 > 100 * cof2_ratio          # yet vastly more multi-harmonic
