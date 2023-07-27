"""
Implementation of the reader for XYZ files using OpenBabel
"""

import bz2
import gzip
import logging
import os
from pathlib import Path
import shutil
import string
import subprocess
import sys
import threading
import time
import re

from openbabel import openbabel
from read_structure_step.formats.registries import register_reader

if "OpenBabel_version" not in globals():
    OpenBabel_version = None

logger = logging.getLogger("read_structure_step.read_structure")


class OutputGrabber(object):
    """Class used to grab standard output or another stream.

    see https://stackoverflow.com/questions/24277488/in-python-how-to-capture-the-stdout-from-a-c-shared-library-to-a-variable/29834357#29834357  # noqa: E501
    """

    escape_char = "\b"

    def __init__(self, stream=None, threaded=False):
        self.origstream = stream
        self.threaded = threaded
        if self.origstream is None:
            self.origstream = sys.stdout
        self.origstreamfd = self.origstream.fileno()
        self.capturedtext = ""
        # Create a pipe so the stream can be captured:
        self.pipe_out, self.pipe_in = os.pipe()

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, type, value, traceback):
        self.stop()

    def start(self):
        """
        Start capturing the stream data.
        """
        self.capturedtext = ""
        # Save a copy of the stream:
        self.streamfd = os.dup(self.origstreamfd)
        # Replace the original stream with our write pipe:
        os.dup2(self.pipe_in, self.origstreamfd)
        if self.threaded:
            # Start thread that will read the stream:
            self.workerThread = threading.Thread(target=self.readOutput)
            self.workerThread.start()
            # Make sure that the thread is running and os.read() has executed:
            time.sleep(0.01)

    def stop(self):
        """
        Stop capturing the stream data and save the text in `capturedtext`.
        """
        # Print the escape character to make the readOutput method stop:
        self.origstream.write(self.escape_char)
        # Flush the stream to make sure all our data goes in before
        # the escape character:
        self.origstream.flush()
        if self.threaded:
            # wait until the thread finishes so we are sure that
            # we have until the last character:
            self.workerThread.join()
        else:
            self.readOutput()
        # Close the pipe:
        os.close(self.pipe_in)
        os.close(self.pipe_out)
        # Restore the original stream:
        os.dup2(self.streamfd, self.origstreamfd)
        # Close the duplicate stream:
        os.close(self.streamfd)

    def readOutput(self):
        """
        Read the stream data (one byte at a time)
        and save the text in `capturedtext`.
        """
        while True:
            char = os.read(self.pipe_out, 1).decode(self.origstream.encoding)
            if not char or self.escape_char in char:
                break
            self.capturedtext += char


def _find_charge(regex, input_file):
    text = re.search(regex, input_file)
    if text is not None:
        return text.group(2)


@register_reader(".xyz")
def load_xyz(
    file_name,
    configuration,
    extension=".xyz",
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
    """Read an XYZ input file.

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


    We'll use OpenBabel to read the file; however, OpenBabel is somewhat limited, so
    we'll first preprocess the file to extract extra data and also to fit it to the
    format that OpenBabel can handle.

    A "standard" .xyz file the following structure:

        #. The number of atoms on the first line
        #. A comment on the second, often the structure name
        #. symbol, x, y, z
        #. ...

    Some XYZ files, for instance those from tmQM encode data in the comment line::

        CSD_code = WIXKOE | q = 0 | S = 0 | Stoichiometry = C47H65LaN2O | MND = 7

    We'll try to handle this type of comment.

    The Minnesota Solvation database uses a slightly modified form:

        #. A comment, often the structure name and provenance
        #. A blank line (maybe the number of atoms, but seems to be blank)
        #. <charge> <multiplicity>
        #. symbol, x, y, z
        #. ...

    so three header lines, and it includes the charge and multiplicity which is very
    useful.

    OpenBabel appears to only work with "standard" files which have the number of atoms
    so this method will transform the MN standard to that.
    """
    global OpenBabel_version

    # Get the text in the file
    if isinstance(file_name, str):
        path = Path(file_name)
    else:
        path = file_name
    path = path.expanduser().resolve()

    # Get the information for progress output, if requested.
    n_structures = 0
    last_line = 0
    with (
        gzip.open(path, mode="rt")
        if path.suffix == ".gz"
        else bz2.open(path, mode="rt")
        if path.suffix == ".bz2"
        else open(path, "r")
    ) as fd:
        for line in fd:
            last_line += 1
            if line.strip() == "":
                n_structures += 1
    # may not have blank line at end
    if line.strip() != "":
        n_structures += 1
    if printer is not None:
        printer("")
        printer(f"    The XYZ file contains {n_structures} structures.")
        last_percent = 0
        t0 = time.time()
        last_t = t0

    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat("xyz")

    configurations = []
    structure_no = 0
    n_errors = 0
    obMol = openbabel.OBMol()

    total_lines = 0
    with (
        gzip.open(path, mode="rt")
        if path.suffix == ".gz"
        else bz2.open(path, mode="rt")
        if path.suffix == ".bz2"
        else open(path, "r")
    ) as fd:
        lines = []
        line_no = 0
        for line in fd:
            # lines have \n at end. It is not stripped as is more normal.
            total_lines += 1
            line_no += 1
            lines.append(line)
            if total_lines == last_line or line_no > 3 and line.strip() == "":
                # End of block, so examine the first lines and see which format
                file_type = "unknown"
                n_lines = len(lines)
                line1 = lines[0].strip()
                fields1 = line1.split()
                n_fields1 = len(fields1)

                if n_lines <= 1:
                    line2 = None
                    fields2 = []
                    n_fields2 = 0
                else:
                    line2 = lines[1].strip()
                    fields2 = line2.split()
                    n_fields2 = len(fields2)

                if n_lines <= 2:
                    line3 = None
                    fields3 = []
                    n_fields3 = 0
                else:
                    line3 = lines[2].strip()
                    fields3 = line3.split()
                    n_fields3 = len(fields3)

                # Check for "standard" file
                if n_fields1 == 1:
                    try:
                        n_atoms = int(fields1[0])
                    except Exception:
                        pass
                    else:
                        # Might be traditional file. Check 3rd line for atom
                        if n_fields3 == 4:
                            file_type = "standard"
                elif n_fields1 == 0:
                    # Might be standard file without atom count.
                    if n_fields3 == 4:
                        file_type = "standard"
                        n_atoms = 0
                        for line in lines[2:]:
                            n_fields = len(line.split())
                            if n_fields == 0:
                                break
                            else:
                                n_atoms += 1
                        # Put the count in line 1
                        lines[0] = str(n_atoms)

                # And Minnesota variant with three headers.
                if n_lines > 3 and n_fields3 == 2:
                    try:
                        charge = int(fields3[0])
                        multiplicity = int(fields3[1])
                    except Exception:
                        pass
                    else:
                        file_type = "Minnesota"
                        if n_fields2 != 0:
                            logger.warning(
                                "Minnesota style XYZ file, 2nd line is not blank:"
                                f"\n\t{lines[1]}"
                            )
                        # Count atoms
                        n_atoms = 0
                        for line in lines[3:]:
                            n_fields = len(line.split())
                            if n_fields == 0:
                                break
                            else:
                                n_atoms += 1
                        # Move comment to 2nd line
                        lines[1] = lines[0]
                        # Put the count in line 1
                        lines[0] = str(n_atoms) + "\n"
                        # Remove 3rd line with charge and multiplicity
                        del lines[2]

                # Reassemble an input file.
                input_data = "".join(lines)

                logger.info(f"Input data:\n\n{input_data}\n")

                title = lines[1].strip()

                lines = []
                line_no = 0

                # Now try to convert using OpenBabel
                out = OutputGrabber(sys.stderr)
                with out:
                    success = obConversion.ReadString(obMol, input_data)
                    if not success:
                        raise RuntimeError("obConversion failed")

                    if add_hydrogens:
                        obMol.AddHydrogens()

                    structure_no += 1
                    if structure_no > 1:
                        if subsequent_as_configurations:
                            configuration = system.create_configuration()
                        else:
                            system = system_db.create_system()
                            configuration = system.create_configuration()

                    configuration.from_OBMol(obMol)
                    configurations.append(configuration)

                # Check any stderr information from obabel.
                if out.capturedtext != "":
                    tmp = out.capturedtext
                    if (
                        "Failed to kekulize aromatic bonds in OBMol::PerceiveBondOrders"
                        not in tmp
                    ):
                        logger.warning(f"{structure_no}: {tmp}")

                if file_type == "Minnesota":
                    # Record the charge, and the spin state
                    configuration.charge = charge
                    configuration.spin_multiplicity = multiplicity

                    logger.info(f"{charge=} {multiplicity=}")
                elif "|" in title:
                    for tmp in title.split("|"):
                        if "=" in tmp:
                            key, val = tmp.split(maxsplit=1)
                            key = key.strip()
                            val = val.strip()
                            if key == "q":
                                try:
                                    configuration.charge = float(val)
                                except Exception:
                                    pass
                            elif key == "S":
                                try:
                                    configuration.spin_multiplicity = int(val)
                                except Exception:
                                    pass
                            elif key == "CSD_code":
                                title = val

                # Set the system name
                if system_name is not None and system_name != "":
                    lower_name = system_name.lower()
                    if "from file" in lower_name:
                        system.name = str(path)
                    elif lower_name == "title":
                        if len(title) > 0:
                            system.name = title
                        else:
                            system.name = str(path)
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
                    elif lower_name == "title":
                        if len(title) > 0:
                            configuration.name = title
                        else:
                            configuration.name = str(path)
                    elif "canonical smiles" in lower_name:
                        configuration.name = configuration.canonical_smiles
                    elif "smiles" in lower_name:
                        configuration.name = configuration.smiles
                    elif lower_name == "sequential":
                        configuration.name = str(structure_no)
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
                                f"\t{structure_no:6} ({percent}%) structures read in "
                                f"{t} seconds. About {t_left} seconds remaining."
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
