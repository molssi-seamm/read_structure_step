"""Extended XYZ: reading over an existing structure (seamm_exec#41) and files
without Properties= (#81)."""

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

WATER = [("O", 0, 0, 0), ("H", 0.96, 0, 0), ("H", -0.24, 0.93, 0)]
METHANE = [
    ("C", 0, 0, 0),
    ("H", 0.63, 0.63, 0.63),
    ("H", -0.63, -0.63, 0.63),
    ("H", -0.63, 0.63, -0.63),
    ("H", 0.63, -0.63, -0.63),
]


def _write(path, atoms, comment):
    with open(path, "w") as fd:
        fd.write(f"{len(atoms)}\n{comment}\n")
        for sym, x, y, z in atoms:
            fd.write(f"{sym} {x} {y} {z}\n")
    return str(path)


@pytest.fixture()
def configuration(monkeypatch):
    if seamm.flowchart_variables is None:
        monkeypatch.setattr(seamm, "flowchart_variables", seamm.Variables())
    db = SystemDB(filename="file:extxyz_overwrite?mode=memory&cache=shared")
    system = db.create_system(name="default")
    configuration = system.create_configuration(name="default")
    yield configuration
    db.close()


def _read(path, configuration):
    return read_structure_step.read(
        path,
        configuration,
        system_db=configuration.system_db,
        system=configuration.system,
        step=STEP,
    )


def test_overwrite_with_a_different_size(configuration, tmp_path):
    """Reading a second file into the same configuration replaces the first
    structure; it failed with an IndexError when the sizes differed."""
    properties = 'Lattice="10 0 0 0 10 0 0 0 10" pbc="T T T" '
    properties += "Properties=species:S:1:pos:R:3"
    _read(_write(tmp_path / "w.extxyz", WATER, properties), configuration)
    assert configuration.n_atoms == 3
    assert configuration.periodicity == 3
    _read(
        _write(tmp_path / "m.extxyz", METHANE, "Properties=species:S:1:pos:R:3"),
        configuration,
    )
    assert configuration.n_atoms == 5
    assert configuration.atoms.symbols == ["C", "H", "H", "H", "H"]
    # The second frame is not periodic, so neither is the configuration now.
    assert configuration.periodicity == 0


def test_no_properties_key(configuration, tmp_path):
    """Without Properties= the columns are species and positions, as ASE does."""
    path = _write(
        tmp_path / "w.extxyz", WATER, 'Lattice="10 0 0 0 10 0 0 0 10" pbc="T T T"'
    )
    confs = _read(path, configuration)
    assert len(confs) == 1
    assert configuration.n_atoms == 3
    assert configuration.periodicity == 3


def test_plain_comment_line(configuration, tmp_path):
    path = _write(tmp_path / "w.extxyz", WATER, "water from somewhere")
    assert len(_read(path, configuration)) == 1
    assert configuration.atoms.symbols == ["O", "H", "H"]


def test_no_structures(configuration, tmp_path):
    path = tmp_path / "empty.extxyz"
    path.write_text("\n")
    with pytest.raises(ValueError, match="contains no structures"):
        _read(str(path), configuration)
