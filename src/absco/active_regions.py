"""The spectral regions in which each molecule is processed (``data/active_regions.csv``).

A band that does not overlap any region of a molecule is skipped for that molecule
(``preprocess.configure.molProcess``), and the regions are what ``absco-generate`` suggests
when a configuration names no molecules.  The file records, for each region, why it is there
(``basis``: the line list has lines, an LBLRTM continuum applies, or cross sections drive it)
and where it came from (``source``: the original 2018 table, or set from the line list).
``absco-active-regions`` checks the file against the AER line list.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass

from absco import paths

FILE = "active_regions.csv"
BASES = ("lines", "continuum", "cross_section")
SOURCES = ("2018", "derived")


@dataclass(frozen=True)
class Region:
    molecule: str
    wn1: float
    wn2: float
    basis: tuple
    source: str
    reference: str = ""
    note: str = ""


def load(path: str | None = None) -> list[Region]:
    """Read the region file (the packaged one by default) and check its structure."""
    path = path or paths.data_file(FILE)
    with open(path, newline="") as fh:
        rows = [r for r in csv.reader(fh) if r and not r[0].startswith("#")]
    head, rows = rows[0], rows[1:]
    expected = ["molecule", "wn1", "wn2", "basis", "source", "reference", "note"]
    if head != expected:
        raise ValueError("%s: header is %s, expected %s" % (path, head, expected))
    out = []
    for n, r in enumerate(rows, start=2):
        mol, wn1, wn2, basis, source, ref, note = (r + [""] * 7)[:7]
        tags = tuple(basis.split("+"))
        if not set(tags) <= set(BASES):
            raise ValueError("%s line %d: unknown basis %r (one or more of %s joined by +)" % (path, n, basis, ", ".join(BASES)))
        if source not in SOURCES:
            raise ValueError("%s line %d: unknown source %r (%s)" % (path, n, source, ", ".join(SOURCES)))
        a, b = float(wn1), float(wn2)
        if not a < b:
            raise ValueError("%s line %d: wn1 %g is not below wn2 %g" % (path, n, a, b))
        out.append(Region(mol.upper(), a, b, tags, source, ref, note))
    return out


def by_molecule(path: str | None = None) -> dict[str, list[list[float]]]:
    """``{molecule: [[wn1, wn2], ...]}``, the form ``preprocess.configure`` uses."""
    regions: dict[str, list[list[float]]] = {}
    for r in load(path):
        regions.setdefault(r.molecule, []).append([r.wn1, r.wn2])
    for mol, rs in regions.items():
        rs.sort()
        for (a1, b1), (a2, b2) in zip(rs, rs[1:]):
            if a2 < b1:
                raise ValueError("%s: regions %g-%g and %g-%g overlap" % (mol, a1, b1, a2, b2))
    return regions
