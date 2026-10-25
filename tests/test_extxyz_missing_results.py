"""Writing extxyz when one structure's calculation failed among good ones.

SEAMM job 5413 (2026-10-25): the other structures had DDEC6 charges, so the
'charge' column existed for all of them, and the failed one's values were None.
The writer stopped with a TypeError after 19 of 36 structures.
"""

import logging

import pytest
from molsystem.system_db import SystemDB

from read_structure_step.formats.extxyz.extxyz import write_extxyz

ENERGY = "DfE0#ORCA#DFT@test"


@pytest.fixture()
def configurations():
    db = SystemDB(filename="file:extxyz_missing?mode=memory&cache=shared")
    db.properties.add(ENERGY, units="kJ/mol", noerror=True)
    result = []
    for name in ("good1", "failed", "good2"):
        system = db.create_system(name=name)
        configuration = system.create_configuration(name=name)
        configuration.atoms.append(
            x=[0.0, 0.0, 0.76], y=[0.0, 0.76, 0.0], z=[0.0, 0.0, 0.0],
            symbol=["O", "H", "H"],
        )  # fmt: skip
        result.append(configuration)
    result[0].atoms.add_attribute("charge", coltype="float", default=None)
    for configuration in (result[0], result[2]):
        configuration.atoms["charge"] = [-0.8, 0.4, 0.4]
        configuration.properties.put(ENERGY, -100.0)
    yield result
    db.close()


def test_every_structure_is_written(configurations, tmp_path):
    path = tmp_path / "out.extxyz"
    write_extxyz(str(path), configurations)
    text = path.read_text()
    assert text.count('SEAMM/configuration_name="') == 3
    frames = text.split('SEAMM/configuration_name="')[1:]
    failed = next(f for f in frames if f.startswith("failed"))
    assert "REF_energy" not in failed.splitlines()[0]
    assert "charge" in text.splitlines()[1]  # the good ones keep their charges


def test_the_failed_structure_is_reported(configurations, tmp_path, caplog):
    printed = []
    with caplog.at_level(logging.WARNING):
        write_extxyz(
            str(tmp_path / "out.extxyz"), configurations, printer=printed.append
        )
    report = "\n".join(printed)
    assert "1 of the 3 structures" in report
    assert "failed: no REF_energy, charges" in report
    assert "failed: no REF_energy, charges" in caplog.text


def test_no_report_when_all_have_the_same_data(configurations, tmp_path):
    printed = []
    write_extxyz(
        str(tmp_path / "out.extxyz"),
        [configurations[0], configurations[2]],
        printer=printed.append,
    )
    assert not any("lack data" in line for line in printed)
