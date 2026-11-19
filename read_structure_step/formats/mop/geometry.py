# -*- coding: utf-8 -*-

"""The geometry of a MOPAC input file as Cartesian coordinates, without running
MOPAC: Cartesians, internal coordinates (a Z-matrix, with dummy atoms) or a mix of
the two, as MOPAC reads them."""

import math
import re

import numpy as np

from molsystem import elements

#: The names MOPAC uses for dummy atoms
DUMMY = ("X", "XX", "99")


class UnsupportedGeometry(ValueError):
    """A geometry this module does not convert (e.g. translation vectors)."""


def _element(token):
    """The element symbol of an atom label, e.g. 'C', 'C(HETATM 1 ...)', '6',
    'Cl', or 'X' for a dummy atom."""
    name = re.sub(r"\(.*\)", "", token).strip()
    if name.upper() in DUMMY:
        return "X"
    if name.isdigit():
        return elements.to_symbols([int(name)])[0]
    match = re.match(r"[A-Za-z]+", name)
    if match is None:
        raise UnsupportedGeometry(f"Cannot read the atom '{token}'")
    symbol = match.group(0).capitalize()
    if symbol == "Tv":
        raise UnsupportedGeometry("Translation vectors (Tv) are not handled")
    if symbol not in elements.symbol_to_atno:
        raise UnsupportedGeometry(f"Unknown element '{token}'")
    return symbol


def _split(line):
    """The atom's label and its numbers. The label may contain blanks inside
    parentheses, e.g. 'N(HETATM    2  N1  HET B   2)'."""
    match = re.match(r"\s*([^\s(]+(?:\([^)]*\))?)(.*)", line)
    return match.group(1), match.group(2).split()


def _number(text):
    return float(text.replace("D", "E").replace("d", "e"))


def _place(a, b, c, r, theta, phi):
    """The position at distance r from a, angle theta (degrees) with a-b, and
    dihedral phi (degrees) with a-b-c."""
    theta = math.radians(theta)
    phi = math.radians(phi)
    bc = a - b
    bc /= np.linalg.norm(bc)
    n = np.cross(b - c, bc)
    if np.linalg.norm(n) < 1.0e-8:
        # a, b and c are collinear: any perpendicular will do
        trial = np.array([0.0, 0.0, 1.0]) if abs(bc[2]) < 0.9 else np.array([1.0, 0, 0])
        n = np.cross(trial, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    d = np.array(
        [
            -r * math.cos(theta),
            r * math.sin(theta) * math.cos(phi),
            r * math.sin(theta) * math.sin(phi),
        ]
    )
    return a + d[0] * bc + d[1] * m + d[2] * n


def mopac_geometry(lines):
    """The atoms and Cartesian coordinates of a MOPAC geometry block.

    Each line is an atom: its label, then three values each followed by an
    optimization flag, then for internal coordinates the atoms it refers to (na,
    nb, nc), and optionally a charge. Internal coordinates are the distance to
    atom na, the angle with na-nb and the dihedral with na-nb-nc. A line without
    those atoms, or with na = 0, has Cartesian coordinates, as MOPAC reads it
    (checked with MOPAC 23): so atom 1, and any atom in a mixed geometry. As in
    MOPAC, an atom with only a distance lies along +x from atom na, and one with
    no dihedral lies in the xy plane on the +y side. Dummy atoms (X, XX, 99) are
    used for the geometry and then dropped.

    Parameters
    ----------
    lines : [str]
        The lines of the geometry block, up to the blank line that ends it.

    Returns
    -------
    ([str], numpy.ndarray)
        The element symbols and the coordinates in Å (n, 3) of the real atoms.

    Raises
    ------
    UnsupportedGeometry
        For what is not converted, such as translation vectors.
    """
    symbols = []
    values = []
    refs = []
    for line in lines:
        if line.strip() == "":
            break
        label, rest = _split(line)
        symbols.append(_element(label))
        # value flag value flag value flag [na nb nc] [charge], or just the
        # three values
        if len(rest) == 3:
            numbers, tail = rest, []
        elif len(rest) >= 6:
            numbers, tail = rest[0:6:2], rest[6:]
        else:
            raise UnsupportedGeometry(f"Cannot read the line '{line}'")
        try:
            values.append([_number(x) for x in numbers])
        except ValueError:
            raise UnsupportedGeometry(f"Cannot read the line '{line}'")
        connect = [0, 0, 0]
        if len(tail) >= 3 and all(re.fullmatch(r"[+-]?\d+", x) for x in tail[:3]):
            connect = [int(x) for x in tail[:3]]
        refs.append(connect)

    n = len(symbols)
    if n == 0:
        raise UnsupportedGeometry("No atoms")

    xyz = np.zeros((n, 3))
    for k in range(n):
        na, nb, nc = refs[k]
        r, theta, phi = values[k]
        if na == 0:
            xyz[k] = values[k]
        elif nb == 0:  # along x from atom na
            xyz[k] = xyz[na - 1] + [r, 0.0, 0.0]
        elif nc == 0:  # in the xy plane, on the +y side
            a = xyz[na - 1]
            b = xyz[nb - 1]
            u = (b - a) / np.linalg.norm(b - a)
            v = np.array([0.0, 1.0, 0.0]) - u[1] * u
            if np.linalg.norm(v) < 1.0e-8:
                v = np.array([0.0, 0.0, 1.0]) - u[2] * u
            v /= np.linalg.norm(v)
            t = math.radians(theta)
            xyz[k] = a + r * (math.cos(t) * u + math.sin(t) * v)
        else:
            for ref in (na, nb, nc):
                if ref < 1 or ref > k:
                    raise UnsupportedGeometry(
                        f"Atom {k + 1} refers to atom {ref}, which is not before it"
                    )
            xyz[k] = _place(xyz[na - 1], xyz[nb - 1], xyz[nc - 1], r, theta, phi)

    keep = [k for k, s in enumerate(symbols) if s != "X"]
    return [symbols[k] for k in keep], xyz[keep]
