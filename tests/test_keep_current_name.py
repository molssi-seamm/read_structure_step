"""'keep current name' keeps a system's or configuration's name, and gives one
without a name the name in the file (its title, else the file's name) -- the
readers stored the words themselves as the name (science's NMS files, all named
'keep current name', 2026-10-09)."""

from types import SimpleNamespace

import pytest
from molsystem.system_db import SystemDB
import seamm

import read_structure_step

STEP = SimpleNamespace(
    parameters=SimpleNamespace(
        current_values_to_dict=lambda context=None: {"save properties": True}
    )
)
WATER = "O 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n"


@pytest.fixture()
def db(monkeypatch):
    if seamm.flowchart_variables is None:
        monkeypatch.setattr(seamm, "flowchart_variables", seamm.Variables())
    db = SystemDB(filename="file:keep_current_name?mode=memory&cache=shared")
    yield db
    db.close()


def _read(path, configuration, **kwargs):
    return read_structure_step.read(
        str(path),
        configuration,
        system_db=configuration.system_db,
        system=configuration.system,
        system_name="keep current name",
        configuration_name="keep current name",
        step=STEP,
        **kwargs,
    )


def test_a_named_system_keeps_its_name(db, tmp_path):
    configuration = db.create_system(name="mine").create_configuration(name="c")
    path = tmp_path / "w.xyz"
    path.write_text(f"3\nwater from a file\n{WATER}")
    (conf,) = _read(path, configuration)
    assert conf.system.name == "mine" and conf.name == "c"


@pytest.mark.parametrize(
    "title, expected", [("DMC cis-cis", "DMC cis-cis"), ("", "dmc_ct")]
)
def test_an_unnamed_system_takes_the_name_in_the_file(db, tmp_path, title, expected):
    configuration = db.create_system().create_configuration()
    path = tmp_path / "dmc_ct.xyz"
    path.write_text(f"3\n{title}\n{WATER}")
    (conf,) = _read(path, configuration)
    assert conf.system.name == expected
    assert conf.system.name != "keep current name"
    assert conf.name != "keep current name"


def test_two_files_into_new_systems_get_distinct_names(db, tmp_path):
    names = []
    for stem in ("dmc_cc", "dmc_ct"):
        configuration = db.create_system().create_configuration()
        path = tmp_path / f"{stem}.xyz"
        path.write_text(f"3\n\n{WATER}")
        (conf,) = _read(path, configuration)
        names.append(conf.system.name)
    assert names == ["dmc_cc", "dmc_ct"]


def test_smiles_file(db, tmp_path):
    configuration = db.create_system().create_configuration()
    path = tmp_path / "methanol.smi"
    path.write_text("CO methanol\n")
    (conf,) = _read(path, configuration)
    assert conf.system.name == "methanol"


SDF_WATER = """{title}
  test

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9600    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2400    0.9300    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END
$$$$
"""


def test_sdf_records_are_named_by_their_titles(db, tmp_path):
    """The configurations of a multi-record SDF without SEAMM's data block are named
    by their title lines (else the file and record), never "None" (science, job
    5363, 2026-10-21: every configuration was named "None")."""
    configuration = db.create_system().create_configuration()
    path = tmp_path / "minima.sdf"
    path.write_text(
        "".join(SDF_WATER.format(title=t) for t in ("min_naw_00_xm", "", "min_naw_02"))
    )
    confs = _read(path, configuration, subsequent_as_configurations=True)
    assert [c.name for c in confs] == ["min_naw_00_xm", "minima_2", "min_naw_02"]


def test_sdf_keeps_a_current_name(db, tmp_path):
    configuration = db.create_system(name="mine").create_configuration(name="c")
    path = tmp_path / "one.sdf"
    path.write_text(SDF_WATER.format(title="from the file"))
    (conf,) = _read(path, configuration)
    assert conf.system.name == "mine" and conf.name == "c"
