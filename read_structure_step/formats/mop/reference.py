# -*- coding: utf-8 -*-

"""The reference data of a MOPAC (.mop) file -- the experimental or calculated
values on its reference line (H=, D=, I=, S=, CP=, GR=, with their sources, HR=,
DR=, IR=, ...) and the MOPAC heat of formation in its comment lines -- as SEAMM
properties, named as in the PM-8 reference set's SDF records (MOPAC session and
Paul, 2026-11-19).

The sources are keys into a BibTeX bibliography, references.bib, in the file's
directory or a parent directory. When there is one, each value also gets where in
the source it is (``,locator``) and how it was derived (``,note``), the cited
entries are kept with the structure (``bibliography``), and a source the
bibliography marks as calculated (``mopac_kind = {calculated}``) tags the values
with its key instead of ``experiment``.
"""

import logging
from pathlib import Path
import re

logger = logging.getLogger(__name__)

#: kJ/mol per kcal/mol, and J/mol/K per cal/mol/K
KCAL = 4.184

#: Reference-line fields -> (SEAMM property, type, units, factor to those units).
#: "{tag}" is "experiment", "estimate", or the key of a calculated source.
PROPERTY = {
    "H": ("enthalpy of formation#{tag}", "float", "kJ/mol", KCAL),
    "CP": ("constant pressure heat capacity#{tag}", "float", "J/mol/K", KCAL),
    "S": ("entropy#{tag}", "float", "J/mol/K", KCAL),
    "D": ("dipole moment#{tag}", "float", "D", 1.0),
    "I": ("ionization energy#{tag}", "float", "eV", 1.0),
}
#: The field giving each value's source
REFERENCE = {"H": "HR", "D": "DR", "I": "IR", "CP": "HR", "S": "HR", "G": "GR"}
#: The field giving each value's weight in the parameterization
WEIGHT = {"H": "HWT", "D": "DWT", "I": "IWT", "G": "GWT"}

BIBLIOGRAPHY_NAME = "references.bib"

_cache = {}  # path -> (mtime, {key: fields}, {key: BibTeX text})


def find_bibliography(path):
    """references.bib in the directory of `path` or the nearest parent, or None."""
    directory = Path(path).expanduser().resolve().parent
    for candidate in (directory, *directory.parents):
        bib = candidate / BIBLIOGRAPHY_NAME
        if bib.is_file():
            return bib
    return None


def load_bibliography(path):
    """The entries of a BibTeX file: ({key: fields}, {key: the entry's text}).

    Read once per file (and again if it changes)."""
    path = Path(path)
    mtime = path.stat().st_mtime
    cached = _cache.get(path)
    if cached is not None and cached[0] == mtime:
        return cached[1], cached[2]

    import bibtexparser

    text = path.read_text(errors="replace")
    parser = bibtexparser.bparser.BibTexParser(common_strings=True)
    parser.ignore_nonstandard_types = False
    entries = {}
    for entry in bibtexparser.loads(text, parser=parser).entries:
        entries[entry["ID"]] = entry
    # Each entry's own text, so a structure carries it exactly as written
    texts = {}
    for match in re.finditer(r"^@\w+\{([^,\s]+),.*?^\}", text, re.S | re.M):
        texts[match.group(1)] = match.group(0)
    _cache[path] = (mtime, entries, texts)
    return entries, texts


def parse_mop(text):
    """The parts of a .mop file: comment lines, keyword lines, title, reference
    line, geometry lines and anything after the geometry."""
    lines = text.splitlines()
    i = 0
    comments = []
    while i < len(lines) and lines[i].lstrip().startswith("*"):
        comments.append(lines[i])
        i += 1
    keyword_lines = [lines[i]] if i < len(lines) else [""]
    n_description = 2
    while True:
        tokens = keyword_lines[-1].split()
        if "&" in tokens:
            n_description -= 1
        elif "+" not in tokens:
            break
        i += 1
        if i >= len(lines):
            break
        keyword_lines.append(lines[i])
        if n_description == 0:
            break
    i += 1
    description = lines[i : i + n_description]
    i += n_description
    geometry = []
    while i < len(lines) and lines[i].strip() != "":
        geometry.append(lines[i])
        i += 1
    trailing = lines[i:]
    while trailing and trailing[-1].strip() == "":
        trailing.pop()
    mopac_hof = None
    if len(comments) >= 2:
        try:
            mopac_hof = float(comments[1].lstrip("* ").split()[0])
        except (ValueError, IndexError):
            mopac_hof = None
    return {
        "comments": comments,
        "keywords": keyword_lines,
        "title": description[0].strip() if description else "",
        "reference line": description[1] if len(description) > 1 else "",
        "geometry": geometry,
        "trailing": trailing,
        "MOPAC heat of formation": mopac_hof,
    }


def parse_reference_line(line):
    """The fields of a reference line, {field: raw value}, and its <...> labels."""
    fields = {}
    for match in re.finditer(r"(?:^|\s)([A-Z][A-Z0-9_]*)=(\S+)", line):
        fields[match.group(1)] = match.group(2)
    labels = re.findall(r"<[^>]*>", line)
    return fields, labels


def number_and_uncertainty(raw):
    """'-73.6,4.8' -> (-73.6, 4.8, []); '-13.9+"FORMIC ACID"' -> (-13.9, None,
    ['FORMIC ACID']): the value, its uncertainty and the species it is relative
    to."""
    match = re.match(r"([-+]?[0-9.]+(?:[EeDd][-+]?\d+)?)(?:,([0-9.]+))?", raw)
    if match is None:
        raise ValueError(f"Cannot read the value '{raw}'")
    value = float(match.group(1).replace("D", "E").replace("d", "e"))
    uncertainty = float(match.group(2)) if match.group(2) else None
    rest = raw[match.end() :]
    relative = [a or b for a, b in re.findall(r'\+"([^"]+)"|\+([^+"\s]+)', rest)]
    return value, uncertainty, relative


def geometry_is_plain_cartesian(parts, n_atoms):
    """Whether Cartesian coordinates carry all that the geometry block does: no
    internal coordinates, dummy atoms, fixed coordinates or trailing sections."""
    if parts["trailing"] or len(parts["geometry"]) != n_atoms:
        return False
    for line in parts["geometry"]:
        tokens = line.split()
        if len(tokens) != 7 or tokens[2] != "+1" or tokens[4] != "+1":
            return False
        if tokens[6] != "+1" or tokens[0].upper() in ("XX", "X", "99", "TV"):
            return False
    return True


def reference_properties(path, n_atoms):
    """The reference data of a .mop file as SEAMM properties.

    Parameters
    ----------
    path : str or pathlib.Path
        The .mop file.
    n_atoms : int
        The number of (real) atoms read from it.

    Returns
    -------
    ([(str, str, str, object)], [str])
        The properties as (name, type, units, value), and notes on what was not
        stored or not found.
    """
    path = Path(path)
    parts = parse_mop(path.read_text(errors="replace"))
    fields, labels = parse_reference_line(parts["reference line"])

    bibliography = find_bibliography(path)
    entries = texts = {}
    if bibliography is not None:
        try:
            entries, texts = load_bibliography(bibliography)
        except Exception as e:
            logger.warning(f"Could not read the bibliography {bibliography}: {e}")
            entries = texts = {}

    def tag_for(key):
        kind = entries.get(key, {}).get("mopac_kind", "") if key else ""
        if kind == "calculated":
            return key
        if kind == "estimate":
            return "estimate"
        return "experiment"

    properties = []
    notes = []
    cited = set()
    if parts["MOPAC heat of formation"] is not None:
        properties.append(
            (
                "enthalpy of formation#MOPAC#PM7",
                "float",
                "kJ/mol",
                round(parts["MOPAC heat of formation"] * KCAL, 6),
            )
        )
    for field, (name, _type, units, factor) in PROPERTY.items():
        # I=, IA=, IE=, IP= and I1= are all ionization energies
        if field == "I":
            keys = [k for k in ("I", "IA", "IE", "IP", "I1") if k in fields]
        else:
            keys = [field] if field in fields else []
        for key in keys:
            ref = fields.get(REFERENCE[field])
            try:
                value, uncertainty, relative = number_and_uncertainty(fields[key])
            except ValueError as e:
                notes.append(str(e))
                continue
            full_name = name.format(tag=tag_for(ref))
            base, method = full_name.split("#", 1)
            properties.append((full_name, _type, units, round(value * factor, 6)))
            if uncertainty is not None:
                properties.append(
                    (
                        f"{base},uncertainty#{method}",
                        _type,
                        units,
                        round(uncertainty * factor, 6),
                    )
                )
            if ref is not None:
                properties.append((f"{base},reference#{method}", "str", "", ref))
                entry = entries.get(ref)
                if entry is not None:
                    cited.add(ref)
                    if entry.get("mopac_locator"):
                        properties.append(
                            (
                                f"{base},locator#{method}",
                                "str",
                                "",
                                entry["mopac_locator"],
                            )
                        )
                    if entry.get("mopac_derivation"):
                        properties.append(
                            (
                                f"{base},note#{method}",
                                "str",
                                "",
                                entry["mopac_derivation"],
                            )
                        )
                elif entries:
                    notes.append(f"The source {ref} is not in {bibliography}.")
            if relative:
                properties.append(
                    (f"{base},relative to#{method}", "json", "", relative)
                )
            weight = WEIGHT.get(field)
            if weight and weight in fields:
                try:
                    properties.append(
                        (f"{base},weight#{method}", "float", "", float(fields[weight]))
                    )
                except ValueError:
                    notes.append(f"Cannot read the weight {weight}={fields[weight]}")
    if "GR" in fields or "GEOREF" in parts["reference line"]:
        properties.append(
            ("geometry,reference#experiment", "str", "", fields.get("GR", ""))
        )
        properties.append(("geometry,labels#experiment", "json", "", labels))
        if "GWT" in fields:
            try:
                properties.append(
                    ("geometry,weight#experiment", "float", "", float(fields["GWT"]))
                )
            except ValueError:
                notes.append(f"Cannot read the weight GWT={fields['GWT']}")
    if cited:
        properties.append(
            (
                "bibliography",
                "str",
                "",
                "\n".join(texts[k] for k in sorted(cited) if k in texts),
            )
        )
    properties.append(
        ("MOPAC keywords", "str", "", "\n".join(k.rstrip() for k in parts["keywords"]))
    )
    properties.append(
        ("MOPAC reference line", "str", "", parts["reference line"].rstrip())
    )
    properties.append(("MOPAC file", "str", "", path.name))
    if not geometry_is_plain_cartesian(parts, n_atoms):
        geometry = parts["geometry"]
        if parts["trailing"]:
            geometry = geometry + [""] + parts["trailing"]
        properties.append(("MOPAC geometry", "str", "", "\n".join(geometry)))
    return properties, notes
