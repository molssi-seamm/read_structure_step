from pathlib import Path
from . import formats
import re

from seamm_util.list_definition import parse_list


def guess_extension(file_name, use_file_name=False):
    """
    Returns the file format. It can either use the file name extension or
    guess based on signatures found in the file.

    Correctly handles .gz and .bz2 files.

    Parameters
    ----------
    file_name: str
        Name of the file

    use_file_name: bool, optional, default: False
        If set to True, uses the file name extension to identify the
        file format.

    Returns
    -------
    extension: str
        The file format.
    """

    if use_file_name is True:
        path = Path(file_name)
        suffixes = path.suffixes
        ext = ""
        if len(suffixes) > 0:
            ext = suffixes[-1]
            if ext in (".gz", ".bz2") and len(suffixes) > 1:
                ext = suffixes[-2]
        if ext == "":
            return None

        return ext.lower()

    available_extensions = formats.registries.REGISTERED_FORMAT_CHECKERS.keys()

    for extension in available_extensions:
        extension_checker = formats.registries.REGISTERED_FORMAT_CHECKERS[extension]

        if extension_checker(file_name) is True:
            return extension


def sanitize_file_format(file_format):
    """
    Returns a uniform file format string.

    Parameters
    ----------
    file_format: str
        Extension of the file.

    Returns
    -------
    file_format: str
        The sanitized file format.
    """

    if re.match(r"^\.?[a-zA-Z\d]+$", file_format) is None:
        raise NameError(
            "read_structure_step: the file format %s could not be validated"
            % file_format
        )

    file_format = file_format.lower()

    if file_format.startswith(".") is False:
        file_format = "." + file_format

    return file_format


def parse_indices(text, maximum):
    """Return the sorted list of 1-based indices in a SEAMM list expression.

    The syntax is that of ``seamm_util.parse_list``: values and ranges
    ``start:stop[:step]`` separated by commas, with the stop value included, e.g.
    ``"1:10:2, 20:end"`` giving 1, 3, 5, 7, 9, 20, ..., ``maximum``. The words
    ``end`` and ``last`` stand for ``maximum``.

    Parameters
    ----------
    text : str
        The list expression.
    maximum : int
        The number of items available (the value of ``end``).

    Returns
    -------
    [int]
        The sorted, unique indices, each between 1 and ``maximum``.
    """
    text = re.sub(r"\b(end|last)\b", str(maximum), str(text).strip(), flags=re.I)
    if text == "":
        raise ValueError("The list of structures to read is empty.")
    values = parse_list(text)
    result = set()
    for value in values:
        if isinstance(value, float):
            if not value.is_integer():
                raise ValueError(f"Structure indices must be integers, not {value}.")
            value = int(value)
        if value < 1 or value > maximum:
            raise ValueError(
                f"Structure index {value} is out of range: there are {maximum} "
                "structures in the file (use 'end' for the last)."
            )
        result.add(value)
    return sorted(result)
