#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for CIF files."""

import pytest  # noqa: F401
import read_structure_step  # noqa: F401
from . import build_filenames

from molsystem.system_db import SystemDB


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
    except Exception:
        print("Caught error deleting the database")


#        "ABENEB01.cif",


@pytest.mark.parametrize(
    "structure",
    [
        "benzene.cif",
    ],
)
def test_cif(configuration, structure):
    """Test the bonds in a benzene crystal."""
    correct = (
        "1.3788 1.3788 1.3788 1.3788 1.3788 1.3788 1.3788 1.3788 "
        "1.3820 1.3820 1.3820 1.3820 1.3820 1.3820 1.3820 1.3820 "
        "0.9300 0.9300 0.9300 0.9300 0.9300 0.9300 0.9300 0.9300 "
        "1.3761 1.3761 1.3761 1.3761 1.3761 1.3761 1.3761 1.3761 "
        "0.9306 0.9306 0.9306 0.9306 0.9306 0.9306 0.9306 0.9306 "
        "0.9296 0.9296 0.9296 0.9296 0.9296 0.9296 0.9296 0.9296"
    )

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

    assert system_db.n_systems == 1
    configuration = system.configuration

    atoms = configuration.atoms
    assert atoms.n_atoms == 48
    bonds = configuration.bonds
    assert bonds.n_bonds == 48

    lengths = " ".join([f"{r:.4f}" for r in bonds.get_lengths()])
    if lengths != correct:
        print(lengths)

    assert lengths == correct
