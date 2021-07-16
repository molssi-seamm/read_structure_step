"""
Implementation of the reader for Tripos MOL2 files using OpenBabel
"""

from pathlib import Path

from openbabel import openbabel

from ..registries import register_format_checker
from ..registries import register_reader
from ..registries import set_format_metadata

set_format_metadata(
    [".mol2"],
    single_structure=False,
    dimensionality=0,
    coordinate_dimensionality=3,
    property_data=True,
    bonds=True,
    is_complete=False,
    add_hydrogens=True,
)


@register_format_checker(".mol2")
def check_format(path):
    """Check if a file is an Tripos MOL2.

    Check for "@<TRIPOS>"

    Parameters
    ----------
    path : str or Path
    """
    result = True
    with open(path, "r") as fd:
        for line in fd:
            line = line.strip()
            if "@<TRIPOS>" in line:
                result = True
                break

    return result


@register_reader(".mol2 -- Tripos MOL2 file")
def load_mol2(
    path,
    configuration,
    extension=".sdf",
    add_hydrogens=True,
    system_db=None,
    system=None,
    indices="1:end",
    subsequent_as_configurations=False,
    system_name="Canonical SMILES",
    configuration_name="sequential",
    **kwargs,
):
    """Read an MDL stucture-data (SDF) file.

    See https://en.wikipedia.org/wiki/Chemical_table_file for a description of the
    format. This function is using Open Babel to handle the file, so trusts that Open
    Babel knows what it is doing.

    Parameters
    ----------
    file_name : str or Path
        The path to the file, as either a string or Path.

    configuration : molsystem.Configuration
        The configuration to put the imported structure into.

    extension : str, optional, default: None
        The extension, including initial dot, defining the format.

    add_hydrogens : bool = True
        Whether to add any missing hydrogen atoms.

    system_db : System_DB = None
        The system database, used if multiple structures in the file.

    system : System = None
        The system to use if adding subsequent structures as configurations.

    indices : str = "1:end"
        The generalized indices (slices, SMARTS, etc.) to select structures
        from a file containing multiple structures.

    subsequent_as_configurations : bool = False
        Normally and subsequent structures are loaded into new systems; however,
        if this option is True, they will be added as configurations.

    system_name : str = "from file"
        The name for systems. Can be directives like "SMILES" or
        "Canonical SMILES". If None, no name is given.

    configuration_name : str = "sequential"
        The name for configurations. Can be directives like "SMILES" or
        "Canonical SMILES". If None, no name is given.

    Returns
    -------
    [Configuration]
        The list of configurations created.
    """
    if isinstance(path, str):
        path = Path(path)

    path.expanduser().resolve()

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol2", "smi")

    configurations = []
    structure_no = 1
    while True:
        if structure_no == 1:
            obMol = openbabel.OBMol()
            not_done = obConversion.ReadFile(obMol, str(path))
        else:
            obMol = openbabel.OBMol()
            not_done = obConversion.Read(obMol)

        if not not_done:
            break

        if add_hydrogens:
            obMol.AddHydrogens()

        if structure_no > 1:
            if subsequent_as_configurations:
                configuration = system.create_configuration()
            else:
                system = system_db.create_system()
                configuration = system.create_configuration()

        configuration.from_OBMol(obMol)
        configurations.append(configuration)

        # Set the system name
        if system_name is not None and system_name != "":
            lower_name = system_name.lower()
            if "from file" in lower_name:
                system.name = obMol.GetTitle()
            elif "canonical smiles" in lower_name:
                system.name = configuration.canonical_smiles
            elif "smiles" in lower_name:
                system.name = configuration.smiles
            else:
                system.name = system_name

        # And the configuration name
        if configuration_name is not None and configuration_name != "":
            lower_name = configuration_name.lower()
            if "from file" in lower_name:
                configuration.name = obMol.GetTitle()
            elif "canonical smiles" in lower_name:
                configuration.name = configuration.canonical_smiles
            elif "smiles" in lower_name:
                configuration.name = configuration.smiles
            elif lower_name == "sequential":
                configuration.name = str(structure_no)
            else:
                configuration.name = configuration_name

        structure_no += 1

    return configurations