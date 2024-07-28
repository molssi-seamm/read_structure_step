"""
Implementation of the reader for Tripos MOL2 files using OpenBabel
"""

from pathlib import Path
import shutil
import string
import subprocess
import time

from openbabel import openbabel

from ..registries import register_format_checker
from ..registries import register_reader
from ..registries import set_format_metadata

if "OpenBabel_version" not in globals():
    OpenBabel_version = None

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
    printer=None,
    references=None,
    bibliography=None,
    **kwargs,
):
    """Read a Tripos MOL2 file.

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

    system_name : str = "title"
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
    global OpenBabel_version

    if isinstance(path, str):
        path = Path(path)

    path.expanduser().resolve()

    # Get the information for progress output, if requested.
    if printer is not None:
        n_structures = 0
        with path.open() as fd:
            for line in fd:
                if line[0:17] == "@<TRIPOS>MOLECULE":
                    n_structures += 1
        printer(f"The Tripos MOL2 file contains {n_structures} structures.")
        last_percent = 0
        t0 = time.time()
        last_t = t0

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
            if lower_name == "title":
                tmp = obMol.GetTitle()
                if tmp != "":
                    system.name = tmp
                else:
                    system.name = f"{path.stem} {structure_no}"
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
                    configuration.name = f"{path.stem} {structure_no}"
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
                configuration.name = str(structure_no)
            else:
                configuration.name = configuration_name

        structure_no += 1
        if printer:
            percent = int(100 * structure_no / n_structures)
            if percent > last_percent:
                t1 = time.time()
                if t1 - last_t >= 60:
                    t = int(t1 - t0)
                    rate = structure_no / (t1 - t0)
                    t_left = int((n_structures - structure_no) / rate)
                    printer(
                        f"\t{structure_no:6} ({percent}%) structures read in {t} "
                        f"seconds. About {t_left} seconds remaining."
                    )
                    last_t = t1
                    last_percent = percent

    if printer:
        t1 = time.time()
        rate = structure_no / (t1 - t0)
        printer(
            f"Read {structure_no} structures in {t1 - t0:.1f} seconds = {rate:.2f} "
            "per second"
        )

    if references:
        # Add the citations for Open Babel
        references.cite(
            raw=bibliography["openbabel"],
            alias="openbabel_jcinf",
            module="read_structure_step",
            level=1,
            note="The principle Open Babel citation.",
        )

        # See if we can get the version of obabel
        if OpenBabel_version is None:
            path = shutil.which("obabel")
            if path is not None:
                path = Path(path).expanduser().resolve()
                try:
                    result = subprocess.run(
                        [str(path), "--version"],
                        stdin=subprocess.DEVNULL,
                        capture_output=True,
                        text=True,
                    )
                except Exception:
                    OpenBabel_version = "unknown"
                else:
                    OpenBabel_version = "unknown"
                    lines = result.stdout.splitlines()
                    for line in lines:
                        line = line.strip()
                        tmp = line.split()
                        if len(tmp) == 9 and tmp[0] == "Open":
                            OpenBabel_version = {
                                "version": tmp[2],
                                "month": tmp[4],
                                "year": tmp[6],
                            }
                        break

        if isinstance(OpenBabel_version, dict):
            try:
                template = string.Template(bibliography["obabel"])

                citation = template.substitute(
                    month=OpenBabel_version["month"],
                    version=OpenBabel_version["version"],
                    year=OpenBabel_version["year"],
                )

                references.cite(
                    raw=citation,
                    alias="obabel-exe",
                    module="read_structure_step",
                    level=1,
                    note="The principle citation for the Open Babel executables.",
                )
            except Exception:
                pass

    return configurations
