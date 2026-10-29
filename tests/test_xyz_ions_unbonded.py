"""A plain XYZ file read with Open Babel leaves alkali and alkaline-earth ions
unbonded, as molsystem's own bond perception does."""

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

# Li+ 1.90 Å from the carbonyl O of formaldehyde (mbe's stage C cluster, 2026-10-29)
LI_FORMALDEHYDE = """5
Li+ and formaldehyde
Li   0.000000   0.000000  -1.900000
O    0.000000   0.000000   0.000000
C    0.000000   0.000000   1.208000
H    0.000000   0.943000   1.795000
H    0.000000  -0.943000   1.795000
"""


@pytest.fixture()
def configuration(monkeypatch):
    if seamm.flowchart_variables is None:
        monkeypatch.setattr(seamm, "flowchart_variables", seamm.Variables())
    db = SystemDB(filename="file:xyz_ions?mode=memory&cache=shared")
    system = db.create_system(name="default")
    configuration = system.create_configuration(name="default")
    yield configuration
    db.close()


def test_li_is_not_bonded(configuration, tmp_path):
    path = tmp_path / "li_formaldehyde.xyz"
    path.write_text(LI_FORMALDEHYDE)
    (conf,) = read_structure_step.read(
        str(path),
        configuration,
        system_db=configuration.system_db,
        system=configuration.system,
        add_hydrogens=False,
        step=STEP,
    )
    symbols = conf.atoms.symbols
    ids = conf.atoms.ids
    li = ids[symbols.index("Li")]
    bonded = {(b["i"], b["j"]) for b in conf.bonds.bonds()}
    assert all(li not in pair for pair in bonded)
    # formaldehyde keeps its own three bonds
    assert conf.bonds.n_bonds == 3
