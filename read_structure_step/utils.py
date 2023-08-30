from pathlib import Path
from . import formats
import re


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
    """Return a list of values in the given index expression.

    Handles expressions like "1-10 by 2, 20-end" which would result in
    1,3,5,7,9,20,21,22,23,24,25 if there were 25 items in the list.
    """
    result = set()
    for indices in text.split(","):
        increment = 1
        if "to" in indices:
            tmp = indices.split("to")
        else:
            if ":" in indices:
                tmp = indices.split(":")
                increment = 0
            else:
                tmp = indices.split("-")
        if len(tmp) == 1:
            if tmp[0].strip() == "end":
                result.add(maximum)
            else:
                result.add(int(tmp[0].strip()))
        else:
            start = int(tmp[0].strip())
            end = tmp[1]
            if "by" in end:
                end, by = end.split("by")
                by = int(by.strip())
            else:
                by = 1
            end = end.strip()
            if end == "end":
                end = maximum
                increment = 1
            else:
                end = int(end)
            result.update(range(start, end + increment, by))
    return sorted(result)
