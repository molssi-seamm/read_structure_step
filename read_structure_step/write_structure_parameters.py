# -*- coding: utf-8 -*-
"""
Control parameters for the Write Structure step in a SEAMM flowchart
"""

import logging

from . import formats
import seamm

logger = logging.getLogger(__name__)


_filetypes = sorted(formats.registries.REGISTERED_WRITERS.keys())


# The pre-2026.9.18 selection ('structures' + 'configurations'), translated on
# loading an older flowchart to the standard SEAMM structure-selection parameters.
def _translate_legacy_selection(data):
    """Translate the legacy 'structures'/'configurations' entries in ``data``."""
    if "structures" not in data:
        return data
    data = dict(data)
    structures = data.pop("structures")
    configs = data.pop("configurations", {"value": "all", "units": None})
    s_value = str(structures.get("value", "current configuration")).strip()
    c_value = str(configs.get("value", "all")).strip()

    def entry(value):
        return {"value": value, "units": None}

    if s_value == "current configuration":
        data["source systems"] = entry("current")
        data["source configurations"] = entry("current")
    elif s_value in ("current system", "all systems"):
        data["source systems"] = entry(
            "current" if s_value == "current system" else "all"
        )
        if c_value == "all":
            data["source configurations"] = entry("all")
        else:
            data["source configurations"] = entry("name is")
            data["source configuration name"] = entry(c_value)
    else:  # a $variable holding the configurations
        data["source systems"] = entry(s_value)
    return data


class WriteStructureParameters(seamm.Parameters):
    """The control parameters for Write Structure

    This is a dictionary of Parameters objects, which themselves are
    dictionaries.  You need to replace the 'time' example below with one or
    more definitions of the control parameters for your plugin and application.

    The fields of each Parameter are:

        default: the default value of the parameter, used to reset it
        kind: one of 'integer', 'float', 'string', 'boolean' or 'enum'
        default_units: the default units, used for reseting the value
        enumeration: a tuple of enumerated values. See below for more.
        format_string: a format string for 'pretty' output
        description: a short string used as a prompt in the GUI
        help_text: a longer string to display as help for the user

    While the 'kind' of a variable might be a numeric value, it may still have
    enumerated values, such as 'normal', 'precise', etc. In addition, any
    parameter can be set to a variable of expression, indicated by having '$'
    as the first character in the field.
    """

    parameters = {
        "file": {
            "default": "",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "File:",
            "help_text": "The file to write.",
        },
        "file type": {
            "default": "from extension",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("from extension", *_filetypes),
            "format_string": "s",
            "description": "Type of file:",
            "help_text": "The type of file, overrides the extension",
        },
        "append": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Append to the file:",
            "help_text": "Whether to append to the file rather than overwrite.",
        },
        "extra attributes": {
            "default": "",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "Extra attributes:",
            "help_text": (
                "Other attributes to add to the 'comment' line in the xyz file."
            ),
        },
        "remove hydrogens": {
            "default": "no",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("no", "nonpolar", "all"),
            "format_string": "s",
            "description": "Remove hydrogens:",
            "help_text": (
                "Whether to remove hydrogen atoms before writing, and if so, just the "
                "nonpolar ones or all."
            ),
        },
        # Which structures to write: the standard SEAMM structure selection
        # (default: the current configuration)
        **seamm.standard_parameters.structure_selection_parameters,
        "number per file": {
            "default": "all",
            "kind": "integer",
            "default_units": "",
            "enumeration": ("all",),
            "format_string": "",
            "description": "# structures per file:",
            "help_text": "The number of structures to write per file.",
        },
        "ignore missing": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": (
                "yes",
                "no",
            ),
            "format_string": "s",
            "description": "Ignore missing structures:",
            "help_text": (
                "Whether a selection that matches no structures is silently ignored "
                "(nothing is written) or is an error."
            ),
        },
    }

    def __init__(self, defaults={}, data=None):
        """
        Initialize the parameters, by default with the parameters defined above

        Args:
            defaults: A dictionary of parameters to initialize. The parameters
                above are used first and any given will override/add to them.
            data: A dictionary of keys and a subdictionary with value and units
                for updating the current, default values.
        """

        logger.debug("WriteStructureParameters.__init__")

        super().__init__(
            defaults={
                **WriteStructureParameters.parameters,
                **defaults,
            },
            data=data,
        )

    def update(self, data):
        """Update from a dictionary, translating the legacy 'structures' and
        'configurations' selection from flowcharts saved before 2026.9.18."""
        super().update(_translate_legacy_selection(data))
