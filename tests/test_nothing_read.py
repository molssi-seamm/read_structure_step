"""A file that gives no structure stops the step, instead of the flowchart going on
with an empty structure (science, 2026-10-20: an SDF entry that could not be read,
then a Liquid Structure step without atoms)."""

from types import SimpleNamespace

import pytest
from molsystem.system_db import SystemDB
import seamm

import read_structure_step
from read_structure_step.read_structure import check_something_read

STEP = SimpleNamespace(
    parameters=SimpleNamespace(
        current_values_to_dict=lambda context=None: {"save properties": True}
    )
)


@pytest.fixture()
def db(monkeypatch):
    if seamm.flowchart_variables is None:
        monkeypatch.setattr(seamm, "flowchart_variables", seamm.Variables())
    db = SystemDB(filename="file:nothing_read?mode=memory&cache=shared")
    yield db
    db.close()


def test_an_unreadable_sdf_entry_reads_nothing(db, tmp_path, monkeypatch):
    configuration = db.create_system(name="s").create_configuration(name="c")
    path = tmp_path / "bad.sdf"
    path.write_text(
        "\n test\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "M  END\n$$$$\n"
    )

    def fail(self, *args, **kwargs):
        raise ValueError("too many values to unpack (expected 4)")

    monkeypatch.setattr(type(configuration), "from_OBMol", fail)
    configurations = read_structure_step.read(
        str(path),
        configuration,
        system_db=db,
        system=configuration.system,
        printer=lambda *a, **k: None,
        step=STEP,
    )
    assert configurations == []
    with pytest.raises(RuntimeError, match="No structure could be read"):
        check_something_read(configurations, configuration, str(path))


def test_something_read_passes(db):
    configuration = db.create_system(name="s").create_configuration(name="c")
    configuration.atoms.append(x=[0.0], y=[0.0], z=[0.0], atno=[8])
    check_something_read([configuration], configuration, "x.sdf")
    check_something_read(None, configuration, "x.mop")
    empty = db.create_system(name="e").create_configuration(name="e")
    with pytest.raises(RuntimeError):
        check_something_read(None, empty, "x.mop")
