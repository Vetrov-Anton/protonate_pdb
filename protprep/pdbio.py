"""Minimal PDB input/output with no external dependencies.

We work at the level of ATOM/HETATM records: that is enough, because all the
"chemistry" is done by pdb2pqr and we only need to slice, splice and rename.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Iterable, List, Tuple

# the 20 standard residues plus every protonated variant pdb2pqr can emit
STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "ASH", "GLH", "LYN", "CYM", "CYX", "HID", "HIE", "HIP", "TYM", "AR0",
    "HSD", "HSE", "HSP", "MSE", "ACE", "NME", "NHE", "NH2",
}

WATER = {"HOH", "WAT", "SOL", "TIP3", "DOD"}


@dataclass
class Atom:
    record: str          # ATOM / HETATM
    serial: int
    name: str            # atom name (unpadded)
    altloc: str
    resname: str
    chain: str
    resseq: int
    icode: str
    x: float
    y: float
    z: float
    occ: float
    bfac: float
    element: str
    charge: str = ""

    # --- convenience accessors ---
    @property
    def res_key(self) -> Tuple[str, int, str]:
        return (self.chain, self.resseq, self.icode)

    @property
    def is_hydrogen(self) -> bool:
        if self.element:
            return self.element.upper() == "H"
        return self.name.lstrip("0123456789").startswith("H")

    def moved(self, xyz) -> "Atom":
        return replace(self, x=float(xyz[0]), y=float(xyz[1]), z=float(xyz[2]))


def _fmt_name(name: str, element: str) -> str:
    """Lay out the atom name in columns 13-16 following the PDB rules."""
    name = name.strip()
    if len(name) >= 4:
        return name[:4]
    el = (element or "").strip()
    # single-letter element -> the name is shifted one position to the right
    if len(el) == 1 and len(name) < 4:
        return f" {name:<3}"
    if not el and not name[0].isdigit() and len(name) < 4:
        return f" {name:<3}"
    return f"{name:<4}"


def parse_pdb(path: str) -> Tuple[List[Atom], List[str]]:
    """Return (atoms, remaining header lines)."""
    atoms: List[Atom] = []
    header: List[str] = []
    for line in open(path, "r", errors="replace"):
        rec = line[:6].strip()
        if rec in ("ATOM", "HETATM"):
            try:
                atoms.append(
                    Atom(
                        record=rec,
                        serial=int(line[6:11]),
                        name=line[12:16].strip(),
                        altloc=line[16],
                        resname=line[17:21].strip(),
                        chain=line[21] if line[21] != " " else "",
                        resseq=int(line[22:26]),
                        icode=line[26] if line[26] != " " else "",
                        x=float(line[30:38]),
                        y=float(line[38:46]),
                        z=float(line[46:54]),
                        occ=float(line[54:60] or 1.0),
                        bfac=float(line[60:66] or 0.0),
                        element=line[76:78].strip(),
                        charge=line[78:80].strip(),
                    )
                )
            except ValueError as err:  # a broken line is worth failing loudly
                raise ValueError(f"Cannot parse PDB line:\n{line}\n{err}")
        elif rec in ("HEADER", "TITLE", "CRYST1", "SSBOND", "LINK", "REMARK"):
            header.append(line.rstrip("\n"))
    return atoms, header


def guess_element(atom: Atom) -> str:
    if atom.element:
        return atom.element
    name = atom.name.strip().lstrip("0123456789")
    for two in ("CL", "BR", "ZN", "MG", "NA", "FE", "MN", "CA", "SE"):
        if name.upper().startswith(two) and len(atom.name.strip()) > 2:
            return two.capitalize()
    return name[0] if name else "X"


def atom_line(atom: Atom, serial: int) -> str:
    el = guess_element(atom)
    return (
        f"{atom.record:<6}{serial % 100000:>5} {_fmt_name(atom.name, el)}"
        f"{atom.altloc or ' '}{atom.resname:<4}{(atom.chain or ' '):>1}"
        f"{atom.resseq % 10000:>4}{atom.icode or ' '}   "
        f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
        f"{atom.occ:6.2f}{atom.bfac:6.2f}          {el.upper():>2}"
        f"{atom.charge:>2}"
    )


def write_pdb(path: str, atoms: Iterable[Atom], title: Iterable[str] = (),
              ter_between_chains: bool = True) -> None:
    serial = 0
    prev = None
    with open(path, "w") as out:
        for line in title:
            out.write(line.rstrip("\n") + "\n")
        for atom in atoms:
            if (
                ter_between_chains
                and prev is not None
                and (atom.chain != prev.chain or prev.record != atom.record)
            ):
                serial += 1
                out.write(
                    f"TER   {serial % 100000:>5}      {prev.resname:<4}"
                    f"{(prev.chain or ' '):>1}{prev.resseq % 10000:>4}"
                    f"{prev.icode or ' '}\n"
                )
            serial += 1
            out.write(atom_line(atom, serial) + "\n")
            prev = atom
        if prev is not None:
            serial += 1
            out.write(
                f"TER   {serial % 100000:>5}      {prev.resname:<4}"
                f"{(prev.chain or ' '):>1}{prev.resseq % 10000:>4}"
                f"{prev.icode or ' '}\n"
            )
        out.write("END\n")


def group_residues(atoms: Iterable[Atom]) -> List[Tuple[Tuple[str, int, str], List[Atom]]]:
    """Group atoms by residue, preserving file order."""
    out: List[Tuple[Tuple[str, int, str], List[Atom]]] = []
    for atom in atoms:
        if out and out[-1][0] == atom.res_key:
            out[-1][1].append(atom)
        else:
            out.append((atom.res_key, [atom]))
    return out
