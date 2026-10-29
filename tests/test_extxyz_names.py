"""Names with '=', spaces, quotes or backslashes in extended XYZ files.

mlff's SEAMM job 5433 (2026-10-29): a configuration named 'LiDMCct_chel_ctmin_toO='
was written by ASE without quotes, and read back as 'LiDMCct_chel_ctmin_toO=pbc=F F
F'. The names this writer puts on the comment line are quoted and escaped, so they
read back as themselves, here and in ASE.
"""

from types import SimpleNamespace

import pytest
from molsystem.system_db import SystemDB
import seamm

from read_structure_step.formats.extxyz.extxyz import load_extxyz, quoted, write_extxyz

NAMES = [
    "LiDMCct_chel_ctmin_toO=",
    "two words",
    'a "quoted" name',
    "back\\slash",
    "plain",
]


def test_quoted():
    assert quoted("plain") == '"plain"'
    assert quoted("a_toO=") == '"a_toO="'
    assert quoted('say "hi"') == '"say \\"hi\\""'
    assert quoted("a\\b") == '"a\\\\b"'


@pytest.fixture()
def db():
    db = SystemDB(filename="file:extxyz_names?mode=memory&cache=shared")
    yield db
    db.close()


def _configurations(db):
    result = []
    for name in NAMES:
        system = db.create_system(name=f"sys {name}")
        configuration = system.create_configuration(name=name)
        configuration.atoms.append(
            x=[0.0, 0.0, 0.76], y=[0.0, 0.76, 0.0], z=[0.0, 0.0, 0.0],
            symbol=["O", "H", "H"],
        )  # fmt: skip
        result.append(configuration)
    return result


def test_names_read_back_by_ase(db, tmp_path):
    ase_io = pytest.importorskip("ase.io")
    path = tmp_path / "names.extxyz"
    write_extxyz(str(path), _configurations(db))
    frames = ase_io.read(str(path), index=":", format="extxyz")
    assert [a.info["SEAMM/configuration_name"] for a in frames] == NAMES
    assert [a.info["SEAMM/system_name"] for a in frames] == [f"sys {n}" for n in NAMES]


STEP = SimpleNamespace(
    parameters=SimpleNamespace(
        current_values_to_dict=lambda context=None: {"save properties": True}
    )
)


def test_names_read_back_here(db, tmp_path, monkeypatch):
    if seamm.flowchart_variables is None:
        monkeypatch.setattr(seamm, "flowchart_variables", seamm.Variables())
    path = tmp_path / "names.extxyz"
    write_extxyz(str(path), _configurations(db))
    other = SystemDB(filename="file:extxyz_names_back?mode=memory&cache=shared")
    try:
        system = other.create_system()
        configuration = system.create_configuration()
        load_extxyz(
            str(path),
            configuration,
            system_db=other,
            system=system,
            subsequent_as_configurations=False,
            system_name="title",
            configuration_name="title",
            perceive_bonds=False,
            step=STEP,
        )
        names = {c.name for s in other.systems for c in s.configurations}
        assert names - {"", None} == set(NAMES)
    finally:
        other.close()
