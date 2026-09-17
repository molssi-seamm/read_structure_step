"""Reading an extended XYZ file perceives the bonds from the geometry."""

import math
from types import SimpleNamespace

import numpy as np
import pytest
from molsystem.system_db import SystemDB
import seamm

import read_structure_step

R_OH = 0.9572
HALF = math.radians(104.52 / 2)


def _water(o, rot=np.eye(3)):
    o = np.asarray(o, dtype=float)
    h1 = rot @ np.array([R_OH * math.sin(HALF), 0.0, R_OH * math.cos(HALF)])
    h2 = rot @ np.array([-R_OH * math.sin(HALF), 0.0, R_OH * math.cos(HALF)])
    return [("O", o), ("H", o + h1), ("H", o + h2)]


def _write_frames(path, frames, L=10.0):
    with open(path, "w") as fd:
        for atoms in frames:
            fd.write(f"{len(atoms)}\n")
            fd.write(
                f'Lattice="{L} 0 0 0 {L} 0 0 0 {L}" pbc="T T T" '
                "Properties=species:S:1:pos:R:3\n"
            )
            for sym, xyz in atoms:
                fd.write(f"{sym} {xyz[0]:.6f} {xyz[1]:.6f} {xyz[2]:.6f}\n")


STEP = SimpleNamespace(
    parameters=SimpleNamespace(
        current_values_to_dict=lambda context=None: {"save properties": True}
    )
)


@pytest.fixture()
def configuration(monkeypatch):
    if seamm.flowchart_variables is None:
        monkeypatch.setattr(seamm, "flowchart_variables", seamm.Variables())
    db = SystemDB(filename="file:extxyz_bonds?mode=memory&cache=shared")
    system = db.create_system(name="default")
    configuration = system.create_configuration(name="default")
    yield configuration
    db.close()


@pytest.fixture()
def water_frames(tmp_path):
    L = 10.0
    # one water wrapped across the far face, one in the middle
    across = _water(
        [L - 0.1, 5.0, 5.0], rot=np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    )
    across = [(s, xyz % L) for s, xyz in across]
    frame1 = across + _water([5.0, 5.0, 5.0])
    frame2 = _water([2.0, 2.0, 2.0]) + _water([7.0, 7.0, 7.0]) + _water([2.0, 7.0, 2.0])
    path = tmp_path / "waters.extxyz"
    _write_frames(path, [frame1, frame2], L=L)
    return path


def test_bonds_perceived_by_default(configuration, water_frames):
    confs = read_structure_step.read(
        str(water_frames),
        configuration,
        system_db=configuration.system_db,
        system=configuration.system,
        subsequent_as_configurations=True,
        step=STEP,
    )
    assert len(confs) == 2
    c1, c2 = confs
    assert c1.periodicity == 3
    assert c1.bonds.n_bonds == 4
    assert c1.find_molecules(as_indices=True) == [[0, 1, 2], [3, 4, 5]]
    # the wrapped water's bonds are through the boundary and still 0.957 A
    assert np.allclose(c1.bonds.get_lengths(), R_OH, atol=1e-4)
    assert c2.bonds.n_bonds == 6
    assert len(c2.find_molecules()) == 3


def test_perception_can_be_turned_off(configuration, water_frames):
    confs = read_structure_step.read(
        str(water_frames),
        configuration,
        system_db=configuration.system_db,
        system=configuration.system,
        subsequent_as_configurations=True,
        step=STEP,
        perceive_bonds=False,
    )
    assert [c.bonds.n_bonds for c in confs] == [0, 0]


def test_summary_line(configuration, water_frames):
    lines = []
    read_structure_step.read(
        str(water_frames),
        configuration,
        system_db=configuration.system_db,
        system=configuration.system,
        subsequent_as_configurations=True,
        step=STEP,
        printer=lines.append,
    )
    assert any(
        "Perceived 10 bonds" in line and "2 structures" in line for line in lines
    )
