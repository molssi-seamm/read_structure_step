# -*- coding: utf-8 -*-

"""Non-graphical part of the Write Structure step in a SEAMM flowchart

In addition to the normal logger, two logger-like printing facilities are
defined: 'job' and 'printer'. 'job' send output to the main job.out file for
the job, and should be used very sparingly, typically to echo what this step
will do in the initial summary of the job.

'printer' sends output to the file 'step.out' in this steps working
directory, and is used for all normal output from this step.
"""

import logging
from pathlib import PurePath
import textwrap

import read_structure_step
from .write import write
import seamm
from seamm import data  # noqa: F401
from seamm_util import ureg, Q_  # noqa: F401
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("Write Structure")


class WriteStructure(seamm.Node):
    def __init__(self, flowchart=None, title="Write Structure", extension=None):
        """A step for Write Structure in a SEAMM flowchart.

        You may wish to change the title above, which is the string displayed
        in the box representing the step in the flowchart.

        Parameters:
            flowchart: The flowchart that contains this step.
            title: The name displayed in the flowchart.

            extension: ??

        Returns:
            None
        """
        logger.debug("Creating Write Structure {}".format(self))

        # Set the logging level for this module if requested
        # if 'write_structure_step_log_level' in self.options:
        #     logger.setLevel(self.options.write_structure_step_log_level)

        super().__init__(
            flowchart=flowchart, title=title, extension=extension, logger=logger
        )  # yapf: disable

        self.parameters = read_structure_step.WriteStructureParameters()

    @property
    def version(self):
        """The semantic version of this module."""
        return read_structure_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return read_structure_step.__git_revision__

    def description_text(self, P=None):
        """Create the text description of what this step will do.
        The dictionary of control values is passed in as P so that
        the code can test values, etc.

        Keyword arguments:
            P: An optional dictionary of the current values of the control
               parameters.
        """

        if not P:
            P = self.parameters.values_to_dict()

        structures = P["structures"]
        configs = P["configurations"]
        ignore_missing = P["ignore missing"]
        if isinstance(ignore_missing, bool):
            if ignore_missing:
                ignore_missing = "yes"
            else:
                ignore_missing = "no"

        n_per_file = P["number per file"]

        if structures == "current configuration":
            text = f"The current configuration will be written to {P['file']}. "
        elif structures == "current system":
            if configs == "all":
                text = (
                    "All the configurations of the current system will be written to "
                    f"{P['file']}. "
                )
            else:
                text = (
                    f"Configuration '{configs}' of the current system will be written "
                    f"to {P['file']}. "
                )
            if n_per_file != "all":
                text += (
                    " The output will be broken into multiple files containing no "
                    f"more than {n_per_file} structures each."
                )
        elif structures == "all systems":
            if configs == "all":
                text = (
                    "All the configurations of all the systems will be written to "
                    f"{P['file']}. "
                )
            else:
                text = (
                    f"Configuration '{configs}' of all systems will be written "
                    f"to {P['file']}. "
                )
                if ignore_missing == "yes":
                    pass
                elif ignore_missing == "no":
                    text += (
                        "It will be an error if a system does not have the named "
                        "configuration."
                    )
                else:
                    text += (
                        f"The value of {ignore_missing} will determine whether it is"
                        " an error if a system does not have the named configuration."
                    )
            if n_per_file != "all":
                text += (
                    " The output will be broken into multiple files containing no "
                    f"more than {n_per_file} structures each."
                )
        else:
            text = (
                f"The variable {structures} will determine which systems to write out, "
                f"and {configs} the configurations for those systems."
            )
            if n_per_file != "all":
                text += (
                    " The output will be broken into multiple files containing no "
                    f"more than {n_per_file} structures each."
                )

        text = textwrap.fill(text, initial_indent=4 * " ", subsequent_indent=4 * " ")
        return self.header + "\n" + text

    def run(self):
        """Run a Write Structure step."""
        next_node = super().run(printer)

        # Get the values of the parameters, dereferencing any variables
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        # What type of file?
        filename = P["file"].strip()
        path = PurePath(filename)
        file_type = P["file type"]

        if file_type != "from extension":
            extension = file_type.split()[0]
        else:
            extension = path.suffix
            if extension == ".gz":
                extension = path.with_suffix("").suffix

        if extension == "":
            raise RuntimeError(
                "Can't write the file without knowing the type (extension)"
            )

        # Print what we are doing
        printer.important(self.description_text(P))

        # Write the file into the system
        system_db = self.get_variable("_system_db")
        system, configuration = self.get_system_configuration(P)

        structures = P["structures"]
        configs = P["configurations"]
        errors = not P["ignore missing"]
        configurations = []
        if structures == "current configuration":
            n_systems = 1
            configurations.append(configuration)
        elif structures == "current system":
            n_systems = 1
            if configs == "all":
                for configuration in system.configurations:
                    configurations.append(configuration)
            else:
                cid = system.get_configuration_id(configs, errors=errors)
                if cid is not None:
                    configurations.append(system.get_configuration(cid))
        elif structures == "all systems":
            n_systems = system_db.n_systems
            if configs == "all":
                for system in system_db.systems:
                    for configuration in system.configurations:
                        configurations.append(configuration)
            else:
                for system in system_db.systems:
                    cid = system.get_configuration_id(configs, errors=errors)
                    if cid is not None:
                        configurations.append(system.get_configuration(cid))

        n_per_file = P["number per file"]
        n_configurations = len(configurations)
        if n_configurations > 1:
            printer.important(
                __(
                    f"\n    Writing out {n_configurations} structures from {n_systems} "
                    "systems.",
                    indent=4 * " ",
                )
            )
        if n_per_file == "all" or n_configurations <= n_per_file:
            write(
                filename,
                configurations,
                extension=extension,
                remove_hydrogens=P["remove hydrogens"],
                printer=printer.important,
                references=self.references,
                bibliography=self._bibliography,
            )
        else:
            n_per_file = int(n_per_file)
            if path.suffix == ".gz":
                base = path.with_suffix("")
                suffix = base.suffix + ".gz"
                stem = str(base.with_suffix(""))
            else:
                suffix = path.suffix
                stem = str(path.with_suffix(""))
            last = 1  # Note that counting from 1 for users.

            while last <= n_configurations:
                first = last
                last += n_per_file
                tmp_name = stem + f"_{first}" + suffix
                write(
                    tmp_name,
                    configurations[first - 1 : last - 1],
                    extension=extension,
                    remove_hydrogens=P["remove hydrogens"],
                    printer=printer.important,
                    references=self.references,
                    bibliography=self._bibliography,
                )

        # Finish the output
        if n_configurations == 1:
            printer.important(
                __(
                    f"\n    Wrote the structure with {configuration.n_atoms} "
                    "atoms."
                    f"\n           System name = {system.name}"
                    f"\n    Configuration name = {configuration.name}",
                    indent=4 * " ",
                )
            )
        printer.important("")

        return next_node
