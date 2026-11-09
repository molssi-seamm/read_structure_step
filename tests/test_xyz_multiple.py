"""Several structures in one XYZ file, in each of its layouts: with atom counts
(blank lines between structures or not), without counts (a blank line in place of
the count, or none at all), and the Minnesota variant."""

import io
from types import SimpleNamespace

import pytest
from molsystem.system_db import SystemDB
import seamm

import read_structure_step
from read_structure_step.formats.xyz.xyz import _records

STEP = SimpleNamespace(
    parameters=SimpleNamespace(
        current_values_to_dict=lambda context=None: {"save properties": True}
    )
)

LIF = ["Li  0.0 0.0 0.0\n", "F   1.56 0.0 0.0\n"]
CO2 = ["C 0.0 0.0 0.0\n", "O -1.16 0.0 0.0\n", "O 1.16 0.0 0.0\n"]


def counted(blank_between):
    sep = "\n" if blank_between else ""
    return (
        "2\nLithium fluoride\n"
        + "".join(LIF)
        + sep
        + "3\nCarbon dioxide\n"
        + "".join(CO2)
    )


# The count line left blank: the blank line before each structure stands for it
BLANK_COUNT = (
    "\nLithium fluoride\n" + "".join(LIF) + "\nCarbon dioxide\n" + "".join(CO2)
)
# No count line at all, structures separated by a blank line
NO_COUNT = "Lithium fluoride\n" + "".join(LIF) + "\nCarbon dioxide\n" + "".join(CO2)
MINNESOTA = (
    "Lithium fluoride\n\n0 1\n"
    + "".join(LIF)
    + "\nCarbon dioxide\n\n0 1\n"
    + "".join(CO2)
)
# A comment line that is just a number does not fool the count
NUMBER_TITLE = "2\n42\n" + "".join(LIF) + "3\n7\n" + "".join(CO2)


@pytest.mark.parametrize(
    "text",
    [counted(False), counted(True), BLANK_COUNT, NO_COUNT, MINNESOTA, NUMBER_TITLE],
    ids=["counts", "counts+blank", "blank-count", "no-count", "minnesota", "number"],
)
def test_records(text):
    records = list(_records(io.StringIO(text)))
    assert len(records) == 2
    atoms = [[x for x in r if len(x.split()) == 4] for r in records]
    assert [len(a) for a in atoms] == [2, 3]


@pytest.fixture()
def db(monkeypatch):
    if seamm.flowchart_variables is None:
        monkeypatch.setattr(seamm, "flowchart_variables", seamm.Variables())
    db = SystemDB(filename="file:xyz_multiple?mode=memory&cache=shared")
    yield db
    db.close()


@pytest.mark.parametrize(
    "text",
    [counted(False), counted(True), BLANK_COUNT, NO_COUNT, MINNESOTA],
    ids=["counts", "counts+blank", "blank-count", "no-count", "minnesota"],
)
def test_read_two_structures(db, tmp_path, text):
    path = tmp_path / "two.xyz"
    path.write_text(text)
    system = db.create_system(name="default")
    configuration = system.create_configuration(name="default")
    configurations = read_structure_step.read(
        str(path),
        configuration,
        system_db=db,
        system=system,
        add_hydrogens=False,
        step=STEP,
        system_name="title",
        configuration_name="title",
    )
    assert [c.n_atoms for c in configurations] == [2, 3]
    assert [c.system.name for c in configurations] == [
        "Lithium fluoride",
        "Carbon dioxide",
    ]


def test_selected_structures(db, tmp_path):
    """Indices count structures, not blank lines."""
    path = tmp_path / "many.xyz"
    path.write_text(counted(False) + "".join(["2\nLiF again\n"] + LIF))
    system = db.create_system(name="default")
    configuration = system.create_configuration(name="default")
    configurations = read_structure_step.read(
        str(path),
        configuration,
        system_db=db,
        system=system,
        add_hydrogens=False,
        step=STEP,
        indices="2:3",
        system_name="title",
        configuration_name="title",
    )
    assert [c.system.name for c in configurations] == ["Carbon dioxide", "LiF again"]
