"""
The cif/mmcif reader/writer
"""

import bz2
import gzip
import logging
from pathlib import Path
import time

from ..registries import register_format_checker
from ..registries import register_reader
from ..registries import register_writer
from ..registries import set_format_metadata
from seamm_util.printing import FormattedText as __
from ...utils import parse_indices

logger = logging.getLogger(__name__)

set_format_metadata(
    [".cif"],
    single_structure=False,
    dimensionality=3,
    coordinate_dimensionality=3,
    property_data=True,
    bonds=True,
    is_complete=True,
    add_hydrogens=False,
)


@register_format_checker(".cif")
def check_format(path):
    """Check if a file is a Crystallographic Information File (CIF) file

    Check for "data_..." at the beginning of a line and no dots in item keys

    Parameters
    ----------
    path : str or Path
    """
    result = False
    in_data_block = False
    with open(path, "r") as fd:
        for line in fd:
            line = line.strip()
            if in_data_block:
                if line[0] == "_":
                    result = "." not in line.split()[0]
                    break
            if line[0:5] == "data_":
                in_data_block = True
    return result


@register_reader(".cif -- Crystallographic Information File")
def load_cif(
    path,
    configuration,
    extension=".cif",
    add_hydrogens=False,
    system_db=None,
    system=None,
    indices="1-end",
    subsequent_as_configurations=False,
    system_name="from file",
    configuration_name="sequential",
    printer=None,
    references=None,
    bibliography=None,
    **kwargs,
):
    """Read a Crystallographic Information File

    See https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax for a description
    of the format.

    Parameters
    ----------
    file_name : str or Path
        The path to the file, as either a string or Path.

    configuration : molsystem.Configuration
        The configuration to put the imported structure into.

    extension : str, optional, default: None
        The extension, including initial dot, defining the format.

    add_hydrogens : bool = False
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

    path = path.expanduser().resolve()

    # Get the information for progress output, if requested.
    n_records = 0
    with (
        gzip.open(path, mode="rt")
        if path.suffix == ".gz"
        else bz2.open(path, mode="rt") if path.suffix == ".bz2" else open(path, "r")
    ) as fd:
        for line in fd:
            if line[0:5] == "data_":
                n_records += 1
    if printer is not None:
        printer("")
        printer(f"    The CIF file contains {n_records} data blocks.")
        last_percent = 0
        t0 = time.time()
        last_t = t0

    # Get the indices to pick
    indices = parse_indices(indices, n_records)
    n_structures = len(indices)
    if n_structures == 0:
        return
    stop = indices[-1]

    configurations = []
    record_no = 0
    structure_no = 0
    n_errors = 0
    lines = []
    in_block = False
    block_name = ""
    with (
        gzip.open(path, mode="rt")
        if path.suffix == ".gz"
        else bz2.open(path, mode="rt") if path.suffix == ".bz2" else open(path, "r")
    ) as fd:
        for line in fd:
            if line[0:5] == "data_":
                logger.debug(f"Found block {line}")
                if not in_block:
                    in_block = True
                    block_name = line[5:].strip()
                else:
                    record_no += 1
                    if record_no > stop:
                        lines = []
                        break
                    if record_no not in indices:
                        lines = []
                        lines.append(line)
                        continue

                    structure_no += 1
                    if structure_no > 1:
                        if subsequent_as_configurations:
                            configuration = system.create_configuration()
                        else:
                            system = system_db.create_system()
                            configuration = system.create_configuration()

                    try:
                        text = configuration.from_cif_text("".join(lines))
                        if text != "":
                            printer("\n")
                            printer(__(text, indent=4 * " "))
                    except Exception as e:
                        n_errors += 1
                        printer("")
                        printer(f" Error handling entry {record_no} in the CIF file:")
                        printer("        " + str(e))
                        printer("    Text of the entry is")
                        printer("    " + 60 * "-")
                        for tmp in lines:
                            printer("    " + tmp.rstrip())
                        printer("    " + 60 * "-")
                        printer("")
                        if n_structures <= 1:
                            raise
                    else:
                        configurations.append(configuration)

                        logger.debug(
                            f"   added system {system_db.n_systems}: {block_name}"
                        )

                        # Set the system name
                        if system_name is not None and system_name != "":
                            lower_name = str(system_name).lower()
                            if "from file" in lower_name:
                                system.name = block_name
                            elif "file name" in lower_name:
                                system.name = path.stem
                            elif "formula" in lower_name:
                                system.name = configuration.formula()[0]
                            elif "empirical formula" in lower_name:
                                system.name = configuration.formula()[1]
                            else:
                                system.name = str(system_name)

                        # And the configuration name
                        if configuration_name is not None and configuration_name != "":
                            lower_name = str(configuration_name).lower()
                            if "from file" in lower_name:
                                configuration.name = block_name
                            elif "file name" in lower_name:
                                configuration.name = path.stem
                            elif "formula" in lower_name:
                                configuration.name = configuration.formula()[0]
                            elif "empirical formula" in lower_name:
                                configuration.name = configuration.formula()[1]
                            else:
                                configuration.name = str(configuration_name)
                        logger.debug(
                            f"   added system {system_db.n_systems}: {block_name}"
                        )

                        if printer:
                            percent = int(100 * structure_no / n_structures)
                            if percent > last_percent:
                                t1 = time.time()
                                if t1 - last_t >= 60:
                                    t = int(t1 - t0)
                                    rate = structure_no / (t1 - t0)
                                    t_left = int((n_structures - structure_no) / rate)
                                    printer(
                                        f"\t{structure_no:6} ({percent}%) structures "
                                        f"read in {t} seconds. About {t_left} seconds "
                                        "remaining."
                                    )
                                    last_t = t1
                                    last_percent = percent
                block_name = line[5:].strip()
                lines = []
            lines.append(line)

        if len(lines) > 0:
            # The last block just ends at the end of the file
            record_no += 1
            if record_no in indices:
                structure_no += 1
                if structure_no > 1:
                    if subsequent_as_configurations:
                        configuration = system.create_configuration()
                    else:
                        system = system_db.create_system()
                        configuration = system.create_configuration()

                text = configuration.from_cif_text("".join(lines))
                logger.debug(f"   added system {system_db.n_systems}: {block_name}")
                if text != "":
                    printer("\n")
                    printer(__(text, indent=4 * " "))

                configurations.append(configuration)

                # Set the system name
                if system_name is not None and system_name != "":
                    lower_name = str(system_name).lower()
                    if "from file" in lower_name:
                        system.name = block_name
                    elif "file name" in lower_name:
                        system.name = path.stem
                    elif "formula" in lower_name:
                        system.name = configuration.formula()[0]
                    elif "empirical formula" in lower_name:
                        system.name = configuration.formula()[1]
                    else:
                        system.name = str(system_name)

                # And the configuration name
                if configuration_name is not None and configuration_name != "":
                    lower_name = str(configuration_name).lower()
                    if "from file" in lower_name:
                        configuration.name = block_name
                    elif "file name" in lower_name:
                        configuration.name = path.stem
                    elif "formula" in lower_name:
                        configuration.name = configuration.formula()[0]
                    elif "empirical formula" in lower_name:
                        configuration.name = configuration.formula()[1]
                    else:
                        configuration.name = str(configuration_name)

    if printer:
        t1 = time.time()
        rate = structure_no / (t1 - t0)
        printer(
            f"    Read {structure_no - n_errors} structures in {t1 - t0:.1f} "
            f"seconds = {rate:.2f} per second"
        )
        if n_errors > 0:
            printer(f"    {n_errors} structures could not be read due to errors.")

    return configurations


@register_writer(".cif -- Crystallographic Information File")
def write_cif(
    path,
    configurations,
    extension=None,
    remove_hydrogens="no",
    printer=None,
    references=None,
    bibliography=None,
):
    """Write to CIF files.

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

    if isinstance(path, str):
        path = Path(path)
    path.expanduser().resolve()

    n_structures = len(configurations)
    last_percent = 0
    last_t = t0 = time.time()
    structure_no = 0
    compress = path.suffix in (".gz", ".bz")
    with (
        gzip.open(path, mode="wb")
        if path.suffix == ".gz"
        else bz2.open(path, mode="wb") if path.suffix == ".bz2" else open(path, "w")
    ) as fd:
        for configuration in configurations:
            text = configuration.to_cif_text()

            structure_no += 1

            if compress:
                fd.write(bytes(text, "utf-8"))
            else:
                fd.write(text)

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
            f"    Wrote {structure_no} structures in {t1 - t0:.1f} seconds = "
            f"{rate:.2f} per second"
        )
    return configurations
