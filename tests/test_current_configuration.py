"""After reading several structures, the step says which one is now current: the last
one read, which the following steps use (2026-11-05: five VPT2 jobs ran on a stray
second SDF record, ethyl vinyl ether, instead of the alkane that came first)."""

from pathlib import Path
from types import SimpleNamespace

import pytest
from molsystem.system_db import SystemDB
import seamm

import read_structure_step
from read_structure_step.read_structure import current_configuration_text

DATA = Path(__file__).parent / "data"
STEP = SimpleNamespace(
    parameters=SimpleNamespace(
        current_values_to_dict=lambda context=None: {"save properties": False}
    )
)


@pytest.fixture()
def db(monkeypatch):
    if seamm.flowchart_variables is None:
        monkeypatch.setattr(seamm, "flowchart_variables", seamm.Variables())
    db = SystemDB(filename="file:current_configuration?mode=memory&cache=shared")
    yield db
    db.close()


def read_two(db, indices="1:end"):
    configuration = db.create_system(name="s").create_configuration(name="c")
    return read_structure_step.read(
        str(DATA / "two_structures.sdf"),
        configuration,
        system_db=db,
        system=configuration.system,
        indices=indices,
        printer=lambda *a, **k: None,
        step=STEP,
    )


def test_the_last_structure_read_is_current(db):
    configurations = read_two(db)
    assert len(configurations) == 2
    system = db.system
    configuration = system.configuration
    assert configuration.formula[0] == "C4 H8 O"
    text = current_configuration_text(configurations, system, configuration)
    assert "is now structure 2 of 2, C4 H8 O" in text
    assert f"Configuration name = {configuration.name}" in text


def test_the_first_structure_only(db):
    configurations = read_two(db, indices="1")
    assert len(configurations) == 1
    assert configurations[0].formula[0] == "C H4"
