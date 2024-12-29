"""
Implementation of the reader for SDF files using OpenBabel
"""

import bz2
import gzip
from pathlib import Path
import shutil
import string
import subprocess
import time

from openbabel import openbabel

from ..registries import register_format_checker
from ..registries import register_reader
from ..registries import register_writer
from ..registries import set_format_metadata
from ...utils import parse_indices

if "OpenBabel_version" not in globals():
    OpenBabel_version = None

set_format_metadata(
    [".sd", ".sdf"],
    single_structure=False,
    dimensionality=0,
    coordinate_dimensionality=3,
    property_data=True,
    bonds=True,
    is_complete=False,
    add_hydrogens=True,
)


@register_format_checker(".sdf")
def check_format(path):
    """Check if a file is an MDL SDFile.

    Check if the last line is "$$$$", which is the terminator for a molecule in SDFiles.

    Parameters
    ----------
    path : str or Path
    """
    last = ""
    with open(path, "r") as fd:
        for line in fd:
            line = line.strip()
            if line != "":
                last = line

    return last == "$$$$"


@register_reader(".sd -- MDL structure-data file")
@register_reader(".sdf -- MDL structure-data file")
def load_sdf(
    path,
    configuration,
    extension=".sdf",
    add_hydrogens=True,
    system_db=None,
    system=None,
    indices="1-end",
    subsequent_as_configurations=False,
    system_name="title",
    configuration_name="sequential",
    printer=None,
    references=None,
    bibliography=None,
    **kwargs,
):
    """Read an MDL structure-data (SDF) file.

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

    indices : str = "1-end"
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

    path = path.expanduser().resolve()

    # Get the information for progress output, if requested.
    n_records = 0
    with (
        gzip.open(path, mode="rt")
        if path.suffix == ".gz"
        else bz2.open(path, mode="rt") if path.suffix == ".bz2" else open(path, "r")
    ) as fd:
        for line in fd:
            if line[0:4] == "$$$$":
                n_records += 1
    if printer is not None:
        printer("")
        printer(f"    The SDF file contains {n_records} structures.")
        last_percent = 0
        t0 = time.time()
        last_t = t0

    # Get the indices to pick
    indices = parse_indices(indices, n_records)
    n_structures = len(indices)
    if n_structures == 0:
        return
    stop = indices[-1]

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("sdf", "smi")

    configurations = []
    record_no = 0
    structure_no = 0
    n_errors = 0
    obMol = openbabel.OBMol()
    text = ""
    with (
        gzip.open(path, mode="rt")
        if path.suffix == ".gz"
        else bz2.open(path, mode="rt") if path.suffix == ".bz2" else open(path, "r")
    ) as fd:
        for line in fd:
            text += line

            if line[0:4] != "$$$$":
                continue

            record_no += 1
            if record_no > stop:
                text = ""
                break
            if record_no not in indices:
                text = ""
                continue

            obConversion.ReadString(obMol, text)

            # See if the system and configuration names are encoded in the title
            title = obMol.GetTitle()
            sysname = title
            confname = title
            have_sysname = False
            if "SEAMM=" in title:
                for tmp in title.split("|"):
                    if "SEAMM=" in tmp and "/" in tmp:
                        sysname, confname = tmp.split("=", 1)[1].split("/", 1)
                        sysname = sysname.strip()
                        confname = confname.strip()
                        have_sysname = True

            if add_hydrogens:
                obMol.AddHydrogens()

            structure_no += 1

            # See if either the system or configuration names are "title"
            if (
                system_name is not None
                and system_name.lower() == "title"
                and have_sysname
            ):
                # Reuse the system if it exists
                if system_db.system_exists(sysname):
                    system = system_db.get_system(sysname)
                elif structure_no > 1:
                    system = system_db.create_system()
                if configuration_name.lower() == "title":
                    names = system.configuration_names
                    if confname in names:
                        cid = system.get_configuration_id(confname)
                        configuration = system.get_configuration(cid)
                    elif structure_no > 1:
                        configuration = system.create_configuration()
            elif structure_no > 1:
                if subsequent_as_configurations:
                    configuration = system.create_configuration()
                else:
                    system = system_db.create_system()
                    configuration = system.create_configuration()

            try:
                configuration.from_OBMol(obMol)
            except Exception as e:
                n_errors += 1
                printer("")
                printer(f"    Error handling entry {record_no} in the SDF file:")
                printer("        " + str(e))
                printer("    Text of the entry is")
                printer("    " + 60 * "-")
                for line in text.splitlines():
                    printer("    " + line)
                printer("    " + 60 * "-")
                printer("")
                text = ""
                continue

            configurations.append(configuration)
            text = ""

            # Set the system name
            if system_name is not None and system_name != "":
                lower_name = system_name.lower()
                if lower_name == "title":
                    if sysname != "":
                        system.name = sysname
                    else:
                        system.name = f"{path.stem}_{record_no}"
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
                    if confname != "":
                        configuration.name = confname
                    else:
                        configuration.name = f"{path.stem}_{record_no}"
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
                    configuration.name = str(record_no)
                elif "formula" in lower_name:
                    configuration.name = configuration.formula[0]
                else:
                    configuration.name = configuration_name

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
            f"    Read {structure_no - n_errors} structures in {t1 - t0:.1f} "
            f"seconds = {rate:.2f} per second"
        )
        if n_errors > 0:
            printer(f"    {n_errors} structures could not be read due to errors.")

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


@register_writer(".sd -- MDL structure-data file")
@register_writer(".sdf -- MDL structure-data file")
def write_sdf(
    path,
    configurations,
    extension=None,
    remove_hydrogens="no",
    printer=None,
    references=None,
    bibliography=None,
):
    """Write an MDL structure-data (SDF) file.

    See https://en.wikipedia.org/wiki/Chemical_table_file for a description of the
    format. This function is using Open Babel to handle the file, so trusts that Open
    Babel knows what it is doing.

    Parameters
    ----------
    path : str
        Name of the file

    configurations : [Configuration]
        The SEAMM configurations to write

    extension : str, optional, default: None
        The extension, including initial dot, defining the format.

    remove_hydrogens : str = "no"
        Whether to remove hydrogen atoms before writing the structure to file.

    printer : Logger or Printer
        A function that prints to the appropriate place, used for progress.

    references : ReferenceHandler = None
        The reference handler object or None

    bibliography : dict
        The bibliography as a dictionary.
    """
    global OpenBabel_version

    if isinstance(path, str):
        path = Path(path)

    path.expanduser().resolve()

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi", "sdf")

    n_structures = len(configurations)
    last_percent = 0
    last_t = t0 = time.time()
    structure_no = 1
    compress = path.suffix == ".gz"
    with gzip.open(path, mode="wb") if compress else open(path, "w") as fd:
        for configuration in configurations:
            obMol = configuration.to_OBMol(properties="*")

            system = configuration.system
            title = f"SEAMM={system.name}/{configuration.name}"
            obMol.SetTitle(title)

            if remove_hydrogens == "nonpolar":
                obMol.DeleteNonPolarHydrogens()
            elif remove_hydrogens == "all":
                obMol.DeleteHydrogens()

            if structure_no == 1:
                text = obConversion.WriteString(obMol)
            else:
                text = obConversion.WriteString(obMol)

            # if not ok
            if text is None or text == "":
                raise RuntimeError("Error writing file")

            if compress:
                fd.write(bytes(text, "utf-8"))
            else:
                fd.write(text)

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
                            f"\t{structure_no:6} ({percent}%) structures wrote in {t} "
                            f"seconds. About {t_left} seconds remaining."
                        )
                        last_t = t1
                        last_percent = percent

    if printer:
        t1 = time.time()
        rate = structure_no / (t1 - t0)
        printer(
            f"    Wrote {structure_no - 1} structures in {t1 - t0:.1f} seconds = "
            f"{rate:.2f} per second"
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
