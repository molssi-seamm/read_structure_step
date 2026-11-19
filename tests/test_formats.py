#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `utils` module."""

import pytest  # noqa: F401
import read_structure_step  # noqa: F401
from . import build_filenames

from molsystem.system_db import SystemDB

bond_string = """\
    i   j  bondorder symop1 symop2
1   1   2          1      .      .
2   1   5          1      .      .
3   1   7          1      .      .
4   2   3          2      .      .
5   3   4          1      .      .
6   3   6          1      .      .
7   4   5          2      .      .
8   5   8          1      .      .
9   6   9          1      .      .
10  6  10          1      .      ."""

xyz_bond_string = """\
    i   j  bondorder symop1 symop2
1   1   7          1      .      .
2   5   8          1      .      .
3   1   5          1      .      .
4   1   2          1      .      .
5   4   5          2      .      .
6   2   3          2      .      .
7   3   4          1      .      .
8   3   6          1      .      .
9   6  10          1      .      .
10  6   9          1      .      ."""


@pytest.fixture()
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


@pytest.mark.parametrize(
    "structure",
    [
        "3TR_model.mol2",
        "3TR_model.xyz",
        "3TR_model_MN.xyz",
        "3TR_model.pdb",
        "3TR_model.sdf",
        "3TR_model.sdf.gz",
        "3TR_model.sdf.bz2",
    ],
)
def test_format(configuration, structure):
    file_name = build_filenames.build_data_filename(structure)
    read_structure_step.read(file_name, configuration)

    assert configuration.n_atoms == 10
    assert all(
        atom in ["N", "N", "N", "N", "C", "C", "H", "H", "H", "H"]
        for atom in configuration.atoms.symbols
    )
    coordinates = configuration.atoms.coordinates
    assert len(coordinates) == 10
    assert all(len(point) == 3 for point in coordinates)
    assert configuration.bonds.n_bonds == 10
    if "xyz" in structure:
        if str(configuration.bonds) != xyz_bond_string:
            print(configuration.bonds)
        assert str(configuration.bonds) == xyz_bond_string
    else:
        if str(configuration.bonds) != bond_string:
            print(configuration.bonds)
        assert str(configuration.bonds) == bond_string


@pytest.mark.parametrize(
    "structure",
    [
        "tmQM_test.xyz.gz",
        "tmQM_test.xyz.bz2",
    ],
)
def test_compressed(configuration, structure):
    file_name = build_filenames.build_data_filename(structure)
    system = configuration.system
    system_db = system.system_db
    read_structure_step.read(
        file_name,
        configuration,
        system_db=system_db,
        system=system,
        subsequent_as_configurations=False,
    )

    assert system_db.n_systems == 7


@pytest.mark.parametrize(
    "structure",
    [
        "tmQM_test.cif",
        "tmQM_test.cif.gz",
        "tmQM_test.cif.bz2",
    ],
)
def test_cif(configuration, structure):
    file_name = build_filenames.build_data_filename(structure)
    system = configuration.system
    system_db = system.system_db
    read_structure_step.read(
        file_name,
        configuration,
        system_db=system_db,
        system=system,
        subsequent_as_configurations=False,
    )

    assert system_db.n_systems == 3
