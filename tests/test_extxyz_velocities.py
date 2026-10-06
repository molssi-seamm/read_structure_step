"""Extended XYZ velocities: read in (eV/amu)^0.5, as the writer writes them.

Reading any file with a velocities column failed: the conversion was from
Å*amu^0.5/eV^0.5 (ASE's unit of time) to Å/fs.
"""

from types import SimpleNamespace

import numpy as np
import pytest
from molsystem.system_db import SystemDB
import seamm
from seamm_util import Q_

import read_structure_step

STEP = SimpleNamespace(
    parameters=SimpleNamespace(
        current_values_to_dict=lambda context=None: {"save properties": True}
    )
)

#: 1 (eV/amu)^0.5 in Å/fs
SQRT_EV_AMU = 0.0982269475


@pytest.fixture()
def configuration(monkeypatch):
    if seamm.flowchart_variables is None:
        monkeypatch.setattr(seamm, "flowchart_variables", seamm.Variables())
    db = SystemDB(filename="file:extxyz_velocities?mode=memory&cache=shared")
    system = db.create_system(name="default")
    configuration = system.create_configuration(name="default")
    yield configuration
    db.close()


def _read(path, configuration):
    return read_structure_step.read(
        str(path),
        configuration,
        system_db=configuration.system_db,
        system=configuration.system,
        step=STEP,
    )


def test_the_unit(configuration):
    assert Q_("eV^0.5/amu^0.5").m_as("Å/fs") == pytest.approx(SQRT_EV_AMU)


def test_velocities_are_read(configuration, tmp_path):
    path = tmp_path / "v.extxyz"
    path.write_text(
        "2\nProperties=species:S:1:pos:R:3:velocities:R:3\n"
        "O 0.0 0.0 0.0 1.0 0.0 0.0\n"
        "H 0.96 0.0 0.0 0.0 -2.0 0.5\n"
    )
    _read(path, configuration)
    velocities = configuration.atoms.get_velocities(fractionals=False, as_array=True)
    expected = np.array([[1.0, 0.0, 0.0], [0.0, -2.0, 0.5]]) * SQRT_EV_AMU
    assert velocities == pytest.approx(expected)


def test_round_trip(configuration, tmp_path):
    """What the writer writes, the reader reads back."""
    configuration.atoms.append(
        x=[0.0, 0.96], y=[0.0, 0.0], z=[0.0, 0.0], symbol=["O", "H"]
    )
    velocities = np.array([[0.01, 0.0, -0.02], [0.0, 0.03, 0.0]])  # Å/fs
    configuration.atoms.set_velocities(velocities.tolist(), fractionals=False)
    path = tmp_path / "out.extxyz"
    read_structure_step.write(str(path), [configuration], extension=".extxyz")

    db = configuration.system_db
    other = db.create_system(name="other").create_configuration(name="other")
    _read(path, other)
    read_back = other.atoms.get_velocities(fractionals=False, as_array=True)
    assert read_back == pytest.approx(velocities)
