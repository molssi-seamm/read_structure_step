#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `read.py` module."""

import pytest  # noqa: F401
import read_structure_step  # noqa: F401
from . import build_filenames

from molsystem.system_db import SystemDB


@pytest.fixture(scope="module")
def configuration():
    """Create a system db, system and configuration."""
    db = SystemDB(filename="file:seamm_db?mode=memory&cache=shared")
    system = db.create_system(name="default")
    configuration = system.create_configuration(name="default")

    yield configuration

    db.close()
    try:
        del db
    except:  # noqa: E722
        print("Caught error deleting the database")


@pytest.mark.parametrize("file_name", [1, [], {}, 1.0])
def test_read_filename_type(configuration, file_name):
    with pytest.raises(TypeError):
        read_structure_step.read(file_name, configuration)


def test_empty_filename(configuration):
    with pytest.raises(NameError):
        read_structure_step.read("", configuration)


def test_unregistered_reader(configuration):
    with pytest.raises(KeyError):
        xyz_file = build_filenames.build_data_filename("spc.xyz")
        read_structure_step.read(xyz_file, configuration, extension=".mp3")


def test_write_extxyz_charges(tmp_path):
    """Per-atom charges are written as a 'charge:R:1' column when present."""
    from molsystem.system_db import SystemDB
    from read_structure_step.formats.extxyz.extxyz import write_extxyz

    db = SystemDB(filename="file:cq_test?mode=memory&cache=shared")
    system = db.create_system(name="s")
    conf = system.create_configuration(name="c")
    conf.atoms.append(symbol=["H", "H"], x=[0.0, 0.0], y=[0.0, 0.0], z=[0.0, 0.74])
    conf.atoms.add_attribute("charge", coltype="float", configuration_dependent=True)
    conf.atoms["charge"][0:] = [-0.13, 0.13]

    path = tmp_path / "out.extxyz"
    write_extxyz(str(path), [conf])
    lines = path.read_text().strip().splitlines()

    assert "charge:R:1" in lines[1]  # header Properties column
    assert lines[2].split()[-1] == f"{-0.13:14.8f}".strip()  # first atom's charge
    assert lines[3].split()[-1] == f"{0.13:14.8f}".strip()  # second atom's charge


def test_write_extxyz_no_charges(tmp_path):
    """No charge column is written when the structure has no charges."""
    from molsystem.system_db import SystemDB
    from read_structure_step.formats.extxyz.extxyz import write_extxyz

    db = SystemDB(filename="file:cq_test2?mode=memory&cache=shared")
    system = db.create_system(name="s")
    conf = system.create_configuration(name="c")
    conf.atoms.append(symbol=["H", "H"], x=[0.0, 0.0], y=[0.0, 0.0], z=[0.0, 0.74])

    path = tmp_path / "out.extxyz"
    write_extxyz(str(path), [conf])
    assert "charge:R:1" not in path.read_text()


def test_extxyz_charges_round_trip(tmp_path, monkeypatch):
    """Per-atom charges written to .extxyz are read back onto the structure."""
    from types import SimpleNamespace
    import seamm
    from molsystem.system_db import SystemDB
    from read_structure_step.formats.extxyz.extxyz import write_extxyz, load_extxyz

    monkeypatch.setattr(seamm, "flowchart_variables", SimpleNamespace(_data={}))

    db = SystemDB(filename="file:rt_test?mode=memory&cache=shared")
    conf = db.create_system(name="s").create_configuration(name="c")
    conf.atoms.append(symbol=["H", "H"], x=[0.0, 0.0], y=[0.0, 0.0], z=[0.0, 0.74])
    conf.atoms.add_attribute("charge", coltype="float", configuration_dependent=True)
    conf.atoms["charge"][0:] = [-0.13, 0.13]

    path = tmp_path / "rt.extxyz"
    write_extxyz(str(path), [conf])

    step = SimpleNamespace(
        parameters=SimpleNamespace(
            current_values_to_dict=lambda context=None: {"save properties": True}
        )
    )
    conf2 = db.create_system(name="s2").create_configuration(name="c2")
    load_extxyz(
        str(path),
        conf2,
        system_db=db,
        system_name=None,
        configuration_name=None,
        step=step,
    )
    assert "charge" in conf2.atoms
    assert list(conf2.atoms["charge"]) == pytest.approx([-0.13, 0.13])
