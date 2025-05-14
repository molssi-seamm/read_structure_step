"""Implementation of the chemical file reader/write using Open Babel"""

from pathlib import Path

from openbabel import openbabel

import molsystem

# Get the list of file formats from Open Babel
obConversion = openbabel.OBConversion()
known_input_formats = obConversion.GetSupportedInputFormat()
known_output_formats = obConversion.GetSupportedOutputFormat()
del obConversion


def load_file(
    path,
    configuration,
    extension=".sdf",
    add_hydrogens=False,
    system_db=None,
    system=None,
    indices="1:end",
    subsequent_as_configurations=False,
    system_name="Canonical SMILES",
    configuration_name="sequential",
    printer=None,
    references=None,
    bibliography=None,
    **kwargs,
):
    """Use Open Babel for reading any of the formats it supports.

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

    printer : Logger or Printer
        A function that prints to the appropriate place, used for progress.

    references : ReferenceHandler = None
        The reference handler object or None

    bibliography : dict
        The bibliography as a dictionary.

    Returns
    -------
    [Configuration]
        The list of configurations created.
    """
    if isinstance(path, str):
        path = Path(path)

    path.expanduser().resolve()

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(extension.lstrip("."), "smi")

    obMol = openbabel.OBMol()
    obConversion.ReadFile(obMol, str(path))

    if add_hydrogens:
        obMol.AddHydrogens()

    configuration.from_OBMol(obMol)

    # Set the system name
    if system_name is not None and system_name != "":
        lower_name = system_name.lower()
        if lower_name == "title":
            tmp = obMol.GetTitle()
            if tmp != "":
                system.name = tmp
            else:
                system.name = path.stem
        elif "canonical smiles" in lower_name:
            system.name = configuration.canonical_smiles
        elif "smiles" in lower_name:
            system.name = configuration.smiles
        elif "iupac" in lower_name:
            system.name = configuration.PC_iupac_name
        elif "inchikey" in lower_name:
            system.name = configuration.inchikey
        elif "inchi" in lower_name:
            system.name = configuration.inchi
        elif "formula" in lower_name:
            system.name = configuration.formula[0]
        else:
            system.name = system_name

    # And the configuration name
    if configuration_name is not None and configuration_name != "":
        lower_name = configuration_name.lower()
        if lower_name == "title":
            tmp = obMol.GetTitle()
            if tmp != "":
                configuration.name = tmp
            else:
                configuration.name = path.stem
        elif "canonical smiles" in lower_name:
            configuration.name = configuration.canonical_smiles
        elif "smiles" in lower_name:
            configuration.name = configuration.smiles
        elif "iupac" in lower_name:
            configuration.name = configuration.PC_iupac_name
        elif "inchikey" in lower_name:
            configuration.name = configuration.inchikey
        elif "inchi" in lower_name:
            configuration.name = configuration.inchi
        elif lower_name == "sequential":
            configuration.name = "1"
        elif "formula" in lower_name:
            configuration.name = configuration.formula[0]
        else:
            configuration.name = configuration_name

    if references:
        # Add the citations for Open Babel
        citations = molsystem.openbabel_citations()
        for i, citation in enumerate(citations, start=1):
            references.cite(
                raw=citation,
                alias=f"openbabel_{i}",
                module="read_structure_step",
                level=1,
                note=f"The principle citation #{i} for OpenBabel.",
            )

    return [configuration]


def write_file(
    path,
    configurations,
    extension=".sdf",
    remove_hydrogens="no",
    printer=None,
    references=None,
    bibliography=None,
    **kwargs,
):
    """Use Open Babel for reading any of the formats it supports.

    See https://en.wikipedia.org/wiki/Chemical_table_file for a description of the
    format. This function is using Open Babel to handle the file, so trusts that Open
    Babel knows what it is doing.

    Parameters
    ----------
    file_name : str or Path
        The path to the file, as either a string or Path.

    configurations : [molsystem.Configuration]
        The configurations to write -- should be one for this module

    extension : str, optional, default: None
        The extension, including initial dot, defining the format.

    remove_hydrogens : str = "no"
        Whether to remove any hydrogen atoms before writing the file.

    printer : Logger or Printer
        A function that prints to the appropriate place, used for progress.

    references : ReferenceHandler = None
        The reference handler object or None

    bibliography : dict
        The bibliography as a dictionary.

    Returns
    -------
    [Configuration]
        The list of configurations created.
    """
    if isinstance(path, str):
        path = Path(path)

    path.expanduser().resolve()

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi", extension.lstrip("."))

    configuration = configurations[0]
    system = configuration.system
    obMol = configuration.to_OBMol()
    title = f"{system.name}/{configuration.name}"
    obMol.SetTitle(title)

    if remove_hydrogens == "nonpolar":
        obMol.DeleteNonPolarHydrogens()
    elif remove_hydrogens == "all":
        obMol.DeleteHydrogens()

    obMol.SetTitle(f"{system.name}/{configuration.name}")

    obConversion.WriteFile(obMol, str(path))

    if references:
        # Add the citations for Open Babel
        citations = molsystem.openbabel_citations()
        for i, citation in enumerate(citations, start=1):
            references.cite(
                raw=citation,
                alias=f"openbabel_{i}",
                module="read_structure_step",
                level=1,
                note=f"The principle citation #{i} for OpenBabel.",
            )

    return [configuration]
