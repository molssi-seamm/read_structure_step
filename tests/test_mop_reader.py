# -*- coding: utf-8 -*-

"""The MOPAC (.mop) reader (2026-11-19): the geometry is converted as MOPAC reads it
(Cartesians, Z-matrices with dummy atoms, or both mixed) without running MOPAC, and
the bonds are perceived from the distances, with metal complexes ionic + dative:
ligand donor -> metal (bond order 9), closed-shell ligands, the metals carrying the
rest of the charge (Paul, 2026-11-19). Before, Open Babel's bond perception gave
e.g. CO as [C+2][O-], and files with dummy atoms needed a MOPAC run."""

from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from molsystem import dative_donor_acceptor
from molsystem.system_db import SystemDB
import read_structure_step
from read_structure_step.formats.mop.geometry import (
    mopac_geometry,
    UnsupportedGeometry,
)

from . import build_filenames


@pytest.fixture()
def configuration():
    db = SystemDB(filename="file:seamm_db?mode=memory&cache=shared")
    system = db.create_system(name="default")
    configuration = system.create_configuration(name="default")
    yield configuration
    db.close()


def _read(configuration, name):
    read_structure_step.read(build_filenames.build_data_filename(name), configuration)
    return configuration


def _datives(configuration):
    """Counter of (donor element, acceptor element) of the dative bonds."""
    symbols = configuration.atoms.symbols
    result = Counter()
    for row in configuration.bonds.bonds():
        pair = dative_donor_acceptor(row["i"], row["j"], row["bondorder"])
        if pair is not None:
            result[(symbols[pair[0] - 1], symbols[pair[1] - 1])] += 1
    return result


def _bonds(configuration):
    symbols = configuration.atoms.symbols
    return Counter(
        tuple(sorted((symbols[r["i"] - 1], symbols[r["j"] - 1])))
        for r in configuration.bonds.bonds()
    )


# -- the geometry, checked against MOPAC 23 (0SCF, AUX) --------------------------


def test_z_matrix_frame_and_dihedral_sign():
    """Atom 2 along +x, atom 3 in the xy plane on the +y side, and the sign of the
    dihedral, as MOPAC places them."""
    symbols, xyz = mopac_geometry(
        [
            "C 0 +0 0 +0 0 +0",
            "C 1.5 +1 0 +0 0 +0 1 0 0",
            "H 1.1 +1 110 +1 0 +0 2 1 0",
            "H 1.1 +1 110 +1 60 +1 1 2 3",
        ]
    )
    assert symbols == ["C", "C", "H", "H"]
    mopac = [[0, 0, 0], [1.5, 0, 0], [1.8762, 1.0337, 0], [-0.3762, 0.5168, -0.8952]]
    assert np.allclose(xyz, mopac, atol=1e-4)
    # atom 3 referring to atoms 1 and 2 is on the +y side too
    _, xyz = mopac_geometry(
        [
            "O 0 +0 0 +0 0 +0",
            "H 0.96 +1 0 +0 0 +0 1 0 0",
            "H 0.96 +1 104.5 +1 0 +0 1 2 0",
        ]
    )
    assert np.allclose(xyz[2], [-0.2404, 0.9294, 0.0], atol=1e-4)


def test_atoms_without_connectivity_are_cartesian():
    """In MOPAC an atom with no connectivity, or 0 0 0, has Cartesian
    coordinates, atoms 2 and 3 included (a mixed geometry)."""
    _, xyz = mopac_geometry(
        [
            "O 0 +0 0 +0 0 +0",
            "H 0.96 +1 0 +0 0 +0 1 0 0",
            "H 0.96 +1 104.5 +1 0 +0 0 0 0",
            "H 1.5 +1 2.0 +1 0.5 +1",
        ]
    )
    assert np.allclose(xyz[2:], [[0.96, 104.5, 0.0], [1.5, 2.0, 0.5]])


def test_dummy_atoms_are_dropped_and_labels_read():
    symbols, xyz = mopac_geometry(
        [
            " XX     1.00000000 +0    0.0000000 +0    0.0000000 +0     0     0     0",
            "  C     1.00000000 +0    0.0000000 +0    0.0000000 +0     1     0     0",
            "  N(HETATM    2  N1  HET B   2)  1.15 +1  90.0 +1  0.0 +0  2  1  0",
        ]
    )
    assert symbols == ["C", "N"]
    assert np.linalg.norm(xyz[1] - xyz[0]) == pytest.approx(1.15)


def test_translation_vectors_are_not_converted():
    with pytest.raises(UnsupportedGeometry):
        mopac_geometry(["C 0 +0 0 +0 0 +0", "Tv 5.0 +0 0.0 +0 0.0 +0"])


# -- reading files ------------------------------------------------------------


def test_acetonitrile(configuration):
    """A Z-matrix with a SYMMETRY section and charges in the last column."""
    _read(configuration, "acetonitrile.mop")
    assert configuration.atoms.symbols == ["C", "C", "H", "H", "H", "N"]
    assert configuration.canonical_smiles == "CC#N"
    assert configuration.bonds.n_bonds == 5


def test_methylidyne_with_a_dummy_atom(configuration):
    """Read without MOPAC: the dummy atom is dropped."""
    _read(configuration, "methylidyne.mop")
    assert configuration.atoms.symbols == ["C", "H"]
    assert configuration.spin_multiplicity == 2
    assert configuration.canonical_smiles == "[CH]"


@pytest.mark.parametrize(
    "name, smiles",
    [
        ("lithium fluoride, trimer.mop", "[F-].[F-].[F-].[Li+].[Li+].[Li+]"),
        ("potassium fluoride, dimer.mop", "[F-].[F-].[K+].[K+]"),
    ],
)
def test_alkali_halides_are_ions(configuration, name, smiles):
    """Alkali and alkaline-earth metals are never bonded by the perception."""
    _read(configuration, name)
    assert configuration.bonds.n_bonds == 0
    assert configuration.canonical_smiles == smiles


def test_chromium_acetate(configuration):
    """Cr2(OAc)4(HOAc)2, mixed Cartesian/internal geometry: eight O -> Cr dative
    bonds, a Cr-Cr single bond, acetates -1 and Cr2+."""
    _read(configuration, "Cr_ACETCR.mop")
    assert configuration.n_atoms == 46
    assert _datives(configuration) == Counter({("O", "Cr"): 8})
    assert _bonds(configuration)[("Cr", "Cr")] == 1
    symbols = configuration.atoms.symbols
    charges = configuration.atoms.get_column_data("formal_charge")
    assert [q for s, q in zip(symbols, charges) if s == "Cr"] == [2, 2]
    assert sum(charges) == 0 == configuration.charge
    assert configuration.spin_multiplicity == 3


def test_europium_complex(configuration):
    """Eu3+ with six N and four O (two kappa-2 acetates) donors; the acetate C,
    reached by the wider metal tolerance, is no donor."""
    _read(configuration, "RM1 63 Europium, CCDC entry LAPJAN.mop")
    assert configuration.n_atoms == 69
    assert _datives(configuration) == Counter({("N", "Eu"): 6, ("O", "Eu"): 4})
    charges = configuration.atoms.get_column_data("formal_charge")
    assert charges[0] == 3 and sum(charges) == 1 == configuration.charge


def test_copper_bromide_cluster(configuration):
    """[Cu4Br10]2-: every Cu-Br bond is dative from Br, and the charges add up.
    The Br...Br contacts of 2.71 Å in this geometry are within 1.2 times the
    covalent radii (2.74 Å), so they are perceived as bonds: a known limit of
    distance-based perception for clusters."""
    _read(configuration, "Cu(II)4Br10(2-) (CIVNAW10).mop")
    assert configuration.atoms.symbols == ["Cu"] * 4 + ["Br"] * 10
    bonds = _bonds(configuration)
    assert _datives(configuration)[("Br", "Cu")] == bonds[("Br", "Cu")] == 15
    assert sum(configuration.atoms.get_column_data("formal_charge")) == -2


# -- the reference data ---------------------------------------------------------

BIB = """\
@misc{NIST,
  title = {NIST Chemistry WebBook},
  mopac_text = {NIST WebBook},
}

@article{LLNBS82,
  author = {Lias, S. G.},
  title = {Ionization energies},
  year = {1982},
  mopac_locator = {Table 2},
  mopac_text = {Lias et al.},
}

@misc{NLM1967,
  title = {Dipole moments},
  mopac_derivation = {from the Stark effect},
  mopac_text = {Nelson, Lide and Maryott},
}

@misc{PW91D,
  title = {PW91 density functional},
  mopac_kind = {calculated},
  mopac_text = {PW91 calculations},
}
"""


def _copy(tmp_path, name, bib=True):
    """The test file in tmp_path/set/normal, with references.bib in tmp_path/set."""
    directory = tmp_path / "set" / "normal"
    directory.mkdir(parents=True)
    target = directory / name
    target.write_bytes(Path(build_filenames.build_data_filename(name)).read_bytes())
    if bib:
        (tmp_path / "set" / "references.bib").write_text(BIB)
    return target


def test_reference_data_with_the_bibliography(configuration, tmp_path):
    """Acetonitrile: D=3.92 I=12.21 IR=LLNBS82 HR=NIST DR=NLM1967 H=17.7 S=58.17
    CP=12.48, in SEAMM's units, with the bibliography in a parent directory."""
    path = _copy(tmp_path, "acetonitrile.mop")
    read_structure_step.read(str(path), configuration)
    got = {k: v["value"] for k, v in configuration.properties.get().items()}
    assert got["enthalpy of formation#experiment"] == pytest.approx(17.7 * 4.184)
    assert got["enthalpy of formation,reference#experiment"] == "NIST"
    assert got["entropy#experiment"] == pytest.approx(58.17 * 4.184)
    assert got["constant pressure heat capacity#experiment"] == pytest.approx(
        12.48 * 4.184
    )
    assert got["dipole moment#experiment"] == pytest.approx(3.92)
    assert got["dipole moment,note#experiment"] == "from the Stark effect"
    assert got["ionization energy#experiment"] == pytest.approx(12.21)
    assert got["ionization energy,locator#experiment"] == "Table 2"
    assert got["MOPAC reference line"].split() == (
        "D=3.92 I=12.21 IR=LLNBS82 HR=NIST DR=NLM1967 H=17.7 S=58.17 CP=12.48".split()
    )
    assert got["MOPAC keywords"].split() == ["SYMMETRY", "PM7"]
    assert got["MOPAC geometry"].startswith("  C ")  # a Z-matrix, kept verbatim
    cited = got["bibliography"]
    assert "@article{LLNBS82," in cited and "@misc{NIST," in cited
    assert "@misc{NLM1967," in cited and "PW91D" not in cited


def test_a_calculated_source_tags_the_values(configuration, tmp_path):
    """[Cu4Br10]2-: H=-248.8 HR=PW91D, which the bibliography marks calculated;
    the MOPAC heat of formation of its comment lines is PM7's."""
    path = _copy(tmp_path, "Cu(II)4Br10(2-) (CIVNAW10).mop")
    read_structure_step.read(str(path), configuration)
    got = {k: v["value"] for k, v in configuration.properties.get().items()}
    assert got["enthalpy of formation#PW91D"] == pytest.approx(-248.8 * 4.184)
    assert got["enthalpy of formation,reference#PW91D"] == "PW91D"
    assert "enthalpy of formation#experiment" not in got
    assert got["enthalpy of formation#MOPAC#PM7"] == pytest.approx(-179.345 * 4.184)


def test_reference_data_without_a_bibliography(configuration, tmp_path):
    """Without references.bib the values and keys are stored, but no locators,
    notes or bibliography, and every source counts as experimental."""
    path = _copy(tmp_path, "Cu(II)4Br10(2-) (CIVNAW10).mop", bib=False)
    read_structure_step.read(str(path), configuration)
    got = {k: v["value"] for k, v in configuration.properties.get().items()}
    assert got["enthalpy of formation,reference#experiment"] == "PW91D"
    assert "bibliography" not in got
    assert not any(",locator#" in k or ",note#" in k for k in got)
