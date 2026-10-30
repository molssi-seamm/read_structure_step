"""SDF files written by SEAMM carry the system and configuration names in their data
block. With several such records and a configuration name that is not taken from the
file (e.g. the default canonical SMILES), every record after the first went into the
first record's configuration, leaving the later systems empty (thermo, ChemAI 5511:
2-methylpyridine's atoms in 2-butanone's system, 2026-10-30)."""

from types import SimpleNamespace

import pytest
from molsystem.system_db import SystemDB
import seamm

import read_structure_step

pytest.importorskip("openbabel")

STEP = SimpleNamespace(
    parameters=SimpleNamespace(
        current_values_to_dict=lambda context=None: {"save properties": True}
    )
)
SPECIES = {"CCC(C)=O": 13, "Cc1ccccn1": 14}


@pytest.fixture()
def db(monkeypatch):
    if seamm.flowchart_variables is None:
        monkeypatch.setattr(seamm, "flowchart_variables", seamm.Variables())
    db = SystemDB(filename="file:sdf_seamm_names?mode=memory&cache=shared")
    yield db
    db.close()


def _seamm_sdf(tmp_path):
    """Two SEAMM-written SDF records, titled 'SEAMM=<SMILES>/initial'."""
    source = SystemDB(filename="file:sdf_seamm_names_source?mode=memory&cache=shared")
    text = ""
    for smiles in SPECIES:
        configuration = source.create_system(name=smiles).create_configuration(
            name="initial"
        )
        configuration.from_smiles(smiles)
        text += configuration.to_sdf_text()
    source.close()
    path = tmp_path / "two.sdf"
    path.write_text(text)
    return path


@pytest.mark.parametrize(
    "system_name, configuration_name",
    [
        ("keep current name", "use Canonical SMILES string"),
        ("title", "sequential"),
        ("title", "title"),
    ],
)
def test_each_record_gets_its_own_system(db, tmp_path, system_name, configuration_name):
    path = _seamm_sdf(tmp_path)
    configuration = db.create_system().create_configuration()
    configurations = read_structure_step.read(
        str(path),
        configuration,
        system_db=db,
        system=configuration.system,
        system_name=system_name,
        configuration_name=configuration_name,
        step=STEP,
    )
    assert len(configurations) == 2
    for conf, (smiles, n_atoms) in zip(configurations, SPECIES.items()):
        assert conf.system.name == smiles
        assert conf.n_atoms == n_atoms
    # Every system holds the atoms of its own species, none is left empty
    for smiles, n_atoms in SPECIES.items():
        system = db.get_system(smiles)
        assert [c.n_atoms for c in system.configurations] == [n_atoms]
