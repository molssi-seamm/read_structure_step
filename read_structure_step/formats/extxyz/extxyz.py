"""
The reader/write for the ASE extended XYZ files, .extxyz
"""

import bz2
import gzip
import logging
from pathlib import Path
import shlex
import time

import numpy as np

from ..registries import register_format_checker
from ..registries import register_reader
from ..registries import register_writer
from ..registries import set_format_metadata
from ...utils import parse_indices
import seamm
from seamm_util import Q_

logger = logging.getLogger(__name__)

set_format_metadata(
    [".extxyz"],
    single_structure=False,
    dimensionality=3,
    coordinate_dimensionality=3,
    property_data=True,
    bonds=True,
    is_complete=True,
    add_hydrogens=False,
    append=True,
    extra_attributes=True,
)


@register_format_checker(".extxyz")
def check_format(path):
    """Check if a file is an ASE style extended XYZ file

    The second line must contain something like

        "Properties=species:S:1:pos:R:3:REF_forces:R:3"

    Check for "Properties="

    Parameters
    ----------
    path : str or Path
    """
    result = False
    with open(path, "r") as fd:
        for lineno, line in enumerate(fd, start=1):
            if lineno == 2:
                if "properties=" in line.lower():
                    result = True
                break
    return result


@register_reader(".extxyz -- ASE style extended XYZ file")
def load_extxyz(
    path,
    configuration,
    extension=".extxyz",
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
    step=None,
    **kwargs,
):
    """Read an extended XYZ File

    See https://www.ovito.org/manual/reference/file_formats/input/xyz.html
    and https://ase-lib.org/ase/io/formatoptions.html#extxyz
    for a description of the format.

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

    step : seamm.Node
        The step running, so we can access parameters, etc.

    Returns
    -------
    [Configuration]
        The list of configurations created.
    """
    # Get the values of the parameters, dereferencing any variables
    P = step.parameters.current_values_to_dict(context=seamm.flowchart_variables._data)
    save_properties = P["save properties"]

    configurations = []

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
            if "Properties=" in line:
                n_records += 1
    if printer is not None:
        printer("")
        printer(f"    The .extxyz file contains {n_records} data blocks.")
        last_percent = 0
        t0 = time.time()
        last_t = t0

    # Get the indices to pick
    indices = parse_indices(indices, n_records)
    n_structures = len(indices)
    if n_structures == 0:
        return configurations
    stop = indices[-1]

    record_no = 0
    structure_no = 0
    line_no = 0
    section = None
    with (
        gzip.open(path, mode="rt")
        if path.suffix == ".gz"
        else bz2.open(path, mode="rt") if path.suffix == ".bz2" else open(path, "r")
    ) as fd:
        for line in fd:
            line_no += 1
            line = line.strip()
            if section is None:
                record_no += 1
                if record_no > stop:
                    break
                try:
                    natoms = int(line)
                except Exception as e:
                    print(e)
                    print(f"{line_no}: {line}")
                    raise
                section = "header"
                atom = 0
            elif section == "header":
                # the header line like
                # Lattice="10.81905 0..." Properties=species:S:1:pos:R:3:REF_forces:R:3
                # REF_energy=-12.0786 REF_stress="0.0054339 ..." pbc="T T T"
                s = shlex.shlex(line, posix=True)
                s.whitespace_split = True
                header = {}
                for tmp in s:
                    key, value = tmp.split("=")
                    if key.startswith("REF_"):
                        key = key[4:]
                    header[key] = value
                if "Properties" not in header:
                    raise ValueError(
                        f"No property definition found in the header line for record "
                        f"{record_no}\n\t{line}"
                    )
                tmp = header["Properties"].replace("REF_", "").split(":")
                meta = [
                    (key, kind, int(n))
                    for key, kind, n in zip(tmp[0::3], tmp[1::3], tmp[2::3])
                ]
                data = {key: [] for key, _, _ in meta}
                nvalues = [n for _, _, n in meta]
                ncolumns = sum(nvalues)
                section = "atoms"
            else:
                tmp = line.split()
                if len(tmp) != ncolumns:
                    raise ValueError(
                        f"There is an error in line {line_no}: \n\t{line}\n\t{meta}"
                    )
                i = 0
                for key, kind, n in meta:
                    i0 = i
                    i += n
                    if n == 1:
                        if kind == "S":
                            data[key].append(tmp[i0])
                        elif kind == "R":
                            data[key].append(float(tmp[i0]))
                        elif kind == "I":
                            data[key].append(int(tmp[i0]))
                    else:
                        if kind == "S":
                            data[key].append(tmp[i0:i])
                        elif kind == "R":
                            data[key].append([float(v) for v in tmp[i0:i]])
                        elif kind == "I":
                            data[key].append([int(v) for v in tmp[i0:i]])
                atom += 1
                if atom == natoms:
                    section = None
                    # Last atom ... create the configuration
                    if record_no not in indices:
                        continue

                    structure_no += 1

                    # See if the system and configuration names are given
                    if "SEAMM/system_name" in header:
                        sysname = header["SEAMM/system_name"]
                        have_sysname = True
                    else:
                        have_sysname = False
                        sysname = None
                    if "SEAMM/configuration_name" in header:
                        confname = header["SEAMM/configuration_name"]
                    else:
                        confname = None

                    # See if either the system or configuration names are "title"
                    if (
                        system_name is not None
                        and system_name.lower() in ("keep current name", "title")
                        and have_sysname
                    ):
                        # Reuse the system if it exists
                        if system_db.system_exists(sysname):
                            system = system_db.get_system(sysname)
                        elif structure_no > 1:
                            system = system_db.create_system(name=sysname)
                        if configuration_name.lower() in ("keep current name", "title"):
                            names = system.configuration_names
                            if confname in names:
                                cid = system.get_configuration_id(confname)
                                configuration = system.get_configuration(cid)
                            elif structure_no > 1:
                                configuration = system.create_configuration(
                                    name=confname
                                )
                    elif structure_no > 1:
                        if subsequent_as_configurations:
                            configuration = system.create_configuration()
                        else:
                            system = system_db.create_system()
                            configuration = system.create_configuration()

                    # Periodic
                    if "pbc" in header and header["pbc"] == "T T T":
                        configuration.periodicity = 3
                        if "Lattice" not in header:
                            raise ValueError(
                                f"record {record_no} ending at line {line_no} is "
                                "periodic, but no Lattice given"
                            )
                        tmp = [float(v) for v in header["Lattice"].split(" ")]
                        vectors = [tmp[i : i + 3] for i in range(0, 9, 3)]
                        configuration.cell.from_vectors(vectors)

                        # Stresses?
                        if save_properties and "stress" in header:
                            prop = (
                                "stress"
                                if "model" not in header
                                else f"stress#{header['model']}"
                            )
                            if not configuration.properties.exists(prop):
                                configuration.properties.add(
                                    prop,
                                    _type="json",
                                    units="MPa",
                                    description="Stress in Voigt notation 6-vector",
                                )
                            units = configuration.properties.units(prop)
                            factor = Q_("eV/Å^3").m_as(units)
                            stress = [
                                float(v) * factor for v in header["stress"].split(" ")
                            ]
                            configuration.properties.put(prop, stress)

                    configuration.atoms.append(symbol=data["species"])
                    configuration.atoms.set_coordinates(data["pos"], fractionals=False)
                    if save_properties and "forces" in data:
                        factor = Q_("eV/Å").m_as("kJ/mol/Å")
                        g = -factor * np.array(data["forces"])
                        configuration.atoms.set_gradients(g, fractionals=False)
                    if save_properties and "velocities" in data:
                        factor = Q_("Å*amu^0.5/eV^0.5").m_as("Å/fs")
                        velocities = factor * np.array(data["velocities"])
                        configuration.atoms.set_velocities(
                            velocities, fractionals=False
                        )

                    # Add other properties from the header
                    if save_properties and "energy" in header:
                        prop = (
                            "energy"
                            if "model" not in header
                            else f"energy#{header['model']}"
                        )
                        if not configuration.properties.exists(prop):
                            configuration.properties.add(
                                prop,
                                _type="float",
                                units="kJ/mol",
                                description="The energy -- usually potential energy",
                            )
                        units = configuration.properties.units(prop)
                        factor = Q_("eV").m_as(units)
                        configuration.properties.put(
                            prop, float(header["energy"]) * factor
                        )

                    configurations.append(configuration)

                    logger.debug(f"   added system {system_db.n_systems}: {record_no}")

                    # Set the system name
                    if system_name is not None and system_name != "":
                        lower_name = system_name.lower()
                        if lower_name in ("keep current name", "title"):
                            if sysname is not None:
                                system.name = sysname
                            else:
                                system.name = f"{path.stem}_{record_no}"
                        elif "canonical smiles" in lower_name:
                            system.name = configuration.canonical_smiles
                        elif "isomeric smiles" in lower_name:
                            system.name = configuration.isomeric_smiles
                        elif "smiles" in lower_name:
                            system.name = configuration.smiles
                        elif "iupac" in lower_name:
                            system.name = configuration.PC_iupac_name
                        elif "inchikey" in lower_name:
                            system.name = configuration.inchikey
                        elif "inchi" in lower_name:
                            system.name = configuration.inchi
                        elif "file name" in lower_name:
                            system.name = path.stem
                        elif "formula" in lower_name:
                            system.name = configuration.formula()[0]
                        elif "empirical formula" in lower_name:
                            system.name = configuration.formula()[1]
                        else:
                            system.name = system_name

                    # And the configuration name
                    if configuration_name is not None and configuration_name != "":
                        lower_name = configuration_name.lower()
                        if lower_name in ("keep current name", "title"):
                            if confname:
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
                        elif "file name" in lower_name:
                            configuration.name = path.stem
                        elif "formula" in lower_name:
                            configuration.name = configuration.formula()[0]
                        elif "empirical formula" in lower_name:
                            configuration.name = configuration.formula()[1]
                        else:
                            configuration.name = configuration_name
                    logger.debug(f"   added system {system_db.n_systems}: {record_no}")

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

    if printer:
        t1 = time.time()
        rate = structure_no / (t1 - t0)
        printer(
            f"    Read {structure_no} structures in {t1 - t0:.1f} "
            f"seconds = {rate:.2f} per second."
        )

        ns = len({c.system.name for c in configurations})
        nc = len({c.name for c in configurations})
        printer(
            f"    {ns} systems and {nc} configurations store the trajectory. The last"
            " is"
        )
        printer(f"               system = {configurations[-1].system.name}")
        printer(f"        configuration = {configurations[-1].name}")

    return configurations


@register_writer(".extxyz -- ASE style extended XYZ file")
def write_extxyz(
    path,
    configurations,
    extension=None,
    remove_hydrogens="no",
    printer=None,
    references=None,
    bibliography=None,
    append=False,
    extra_attributes="",
    **kwargs,
):
    """Write to .extxyz files.

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

    append : bool
        Whether to append to the file.

    extra_attributes : str
        Extra attributes of the configuration as a string with key="value" ...
        Where quotes are needed if the value contains blanks.
    """

    if isinstance(path, str):
        path = Path(path)
    path = path.expanduser().resolve()

    n_structures = len(configurations)
    last_percent = 0
    last_t = t0 = time.time()
    structure_no = 0

    mode = "a" if append else "w"
    with (
        gzip.open(path, mode=mode + "t")
        if path.suffix == ".gz"
        else (
            bz2.open(path, mode=mode + "t")
            if path.suffix == ".bz2"
            else open(path, mode)
        )
    ) as fd:
        for configuration in configurations:
            natoms = configuration.n_atoms
            text = [str(natoms)]
            header = "Properties=species:S:1:pos:R:3"
            model = ""

            # Other columns for the atoms...
            have_gradients = "no"
            if configuration.atoms.have_gradients:
                # In the atom tables
                header += ":REF_forces:R:3"
                have_gradients = "atoms"
            else:
                available = configuration.properties.list("gradients*")
                if len(available) > 0:
                    # If there are more than one, take the latest as it is the most
                    # likely one the user wants.
                    key = available[-1]
                    if "#" in key:
                        model = key.split("#", maxsplit=1)[1]
                    header += ":REF_forces:R:3"
                    have_gradients = key

            have_velocities = "no"
            if configuration.atoms.have_velocities:
                # In the atom tables
                header += ":velocities:R:3"
                have_velocities = "atoms"
            else:
                available = configuration.properties.list("velocities*")
                if len(available) > 0:
                    # If there are more than one, take the latest as it is the most
                    # likely one the user wants.
                    header += ":velocities:R:3"
                    have_velocities = available[-1]

            # See if the energy exists as a property. May be "potential energy"
            for prop in (
                "DfE0*",
                "enthalpy of formation*",
                "potential energy*",
                "E#*",
                "total energy*",
                "energy*",
            ):
                available = configuration.properties.list(prop)
                if len(available) > 0:
                    # If there are more than one, take the latest as it is the most
                    # likely one the user wants.
                    key = available[-1]
                    if "#" in key:
                        model = key.split("#", maxsplit=1)[1]
                    E = configuration.properties.get(key)[key]["value"]
                    units = configuration.properties.units(available[-1])
                    E *= Q_(units).m_as("eV")
                    header += f" REF_energy={E:.5f}"
                    break

            if configuration.periodicity == 3:
                header += ' pbc="T T T" Lattice="'
                vectors = configuration.cell.vectors()
                vectors = [f"{v:.5f}" for row in vectors for v in row]
                header += " ".join(vectors)
                header += '"'

                # See if the stress exists as a property
                available = configuration.properties.list("stress*")
                if len(available) > 0:
                    # If there are more than one, take the latest as it is the most
                    # likely one the user wants.
                    key = available[-1]
                    if "#" in key:
                        model = key.split("#", maxsplit=1)[1]
                    stresses = configuration.properties.get(key)[key]["value"]
                    units = configuration.properties.units(key)
                    header += ' REF_stress="'
                    factor = Q_(units).m_as("eV/Å^3")
                    stresses = [f"{factor * v:.7f}" for v in stresses]
                    header += " ".join(stresses)
                    header += '"'

            if "model=" not in extra_attributes and model != "":
                header += f' model="{model}"'

            header += f' SEAMM/system_name="{configuration.system.name}"'
            header += f' SEAMM/configuration_name="{configuration.name}"'

            if extra_attributes != "":
                header += " "
                header += extra_attributes

            text.append(header)

            # Now for the atoms, potentially with forces and velocities
            symbols = configuration.atoms.symbols
            xyzs = configuration.atoms.get_coordinates(fractionals=False)

            if have_gradients == "no":
                forces = [[0.0] * 3 for _ in range(natoms)]
            elif have_gradients == "atoms":
                gradients = configuration.atoms.get_gradients(
                    fractionals=False, as_array=True
                )
                factor = Q_("kJ/mol/Å").m_as("eV/Å")
                forces = (-gradients * factor).tolist()
            else:
                gradients = configuration.properties.get(have_gradients)[
                    have_gradients
                ]["value"]
                units = configuration.properties.units(have_gradients)
                factor = Q_(units).m_as("eV/Å")
                forces = (-np.array(gradients) * factor).tolist()

            if have_velocities == "no":
                velocities = [[0.0] * 3 for _ in range(natoms)]
            elif have_velocities == "atoms":
                velocities = configuration.atoms.get_velocities(
                    fractionals=False, as_array=True
                )
                factor = Q_("Å/fs").m_as("eV^0.5/amu^0.5")
                velocities = (velocities * factor).tolist()
            else:
                velocities = configuration.properties.get(have_velocities)[
                    have_velocities
                ]["value"]
                units = configuration.properties.units(have_velocities)
                factor = Q_(units).m_as("eV^0.5/amu^0.5")
                velocities = (np.array(velocities) * factor).tolist()

            for symbol, xyz, force, velocity in zip(symbols, xyzs, forces, velocities):
                line = f"{symbol:<2} "
                line += " ".join([f"{v:14.8f}" for v in xyz])

                if have_gradients != "no":
                    line += " "
                    line += " ".join([f"{v:14.8f}" for v in force])

                if have_velocities != "no":
                    line += " "
                    line += " ".join([f"{v:14.8f}" for v in velocity])

                text.append(line)
            text = "\n".join(text)

            structure_no += 1

            fd.write(text)
            fd.write("\n")

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
