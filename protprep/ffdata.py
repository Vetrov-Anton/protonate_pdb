"""Reading GROMACS force field data (*.rtp, *.r2b) and matching atom names.

The idea: pdb2pqr emits structures in AMBER nomenclature (HB2/HB3, CD1 in ILE,
OXT...) while the GROMACS amber99sb-ildn port expects its own (HB1/HB2, CD,
OC1/OC2...). Instead of a hard-coded table we read the rtp files from the force
field directory and derive the mapping automatically - that way the script
survives any local edit of the force field.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

SECTION_KEYWORDS = {
    "atoms", "bonds", "angles", "dihedrals", "impropers", "exclusions",
    "cmap", "bondedtypes",
}

# Explicit correspondences that cannot be derived automatically (symmetric pairs).
EXPLICIT_RENAME: Dict[str, Dict[str, str]] = {
    "*": {"OXT": "OC1", "O''": "OC1"},
}


@dataclass
class RtpEntry:
    name: str
    # name, type, charge, cgnr
    atoms: List[Tuple[str, str, float, int]] = field(default_factory=list)
    bonds: List[Tuple[str, str]] = field(default_factory=list)
    impropers: List[Tuple[str, ...]] = field(default_factory=list)

    @property
    def charge(self) -> float:
        return sum(a[2] for a in self.atoms)

    @property
    def atom_names(self) -> List[str]:
        return [a[0] for a in self.atoms]

    def parents(self) -> Dict[str, str]:
        """For every hydrogen, the heavy atom it is bonded to."""
        heavy = {n for n in self.atom_names if not _is_h(n)}
        out: Dict[str, str] = {}
        for a, b in self.bonds:
            for x, y in ((a, b), (b, a)):
                if _is_h(x) and y in heavy:
                    out[x] = y
        return out


def _is_h(name: str) -> bool:
    return name.lstrip("0123456789+-").upper().startswith("H")


def gromacs_top_dirs() -> List[str]:
    """Top directories of the installed GROMACS where *.ff live."""
    out = []
    gmxdata = os.environ.get("GMXDATA")
    if gmxdata:
        out.append(os.path.join(gmxdata, "top"))
    out += [
        "/usr/local/gromacs/share/gromacs/top",
        "/usr/share/gromacs/top",
        "/usr/local/share/gromacs/top",
    ]
    return [d for d in out if os.path.isdir(d)]


def resolve_ff_path(name: str) -> str:
    """Locate a force field directory: as given, nearby, then inside GROMACS.

    Accepts both a path and a bare name, with or without the .ff suffix.
    """
    variants = [name, name + ".ff"] if not name.endswith(".ff") else [name]
    for variant in variants:
        if os.path.isdir(variant):
            return os.path.abspath(variant)
    for top in gromacs_top_dirs():
        for variant in variants:
            cand = os.path.join(top, os.path.basename(variant))
            if os.path.isdir(cand):
                return cand
    available = sorted(
        {
            entry
            for top in gromacs_top_dirs()
            for entry in os.listdir(top)
            if entry.endswith(".ff")
        }
    )
    hint = ""
    if available:
        hint = ("\nThe installed GROMACS provides: " + ", ".join(available)
                + "\nPick one with --ff (a bare name is fine).")
    raise FileNotFoundError(f"Force field directory not found: {name}{hint}")


class ForceField:
    """Lazy wrapper around a <name>.ff directory."""

    def __init__(self, path: str):
        self.path = resolve_ff_path(path)
        self.rtp: Dict[str, RtpEntry] = {}
        for fname in sorted(os.listdir(self.path)):
            if fname.endswith(".rtp"):
                self._read_rtp(os.path.join(self.path, fname))
        self.r2b: Dict[str, Tuple[str, str, str]] = {}
        self.ff2gmx: Dict[str, str] = {}
        r2b_path = os.path.join(self.path, "aminoacids.r2b")
        if os.path.exists(r2b_path):
            self._read_r2b(r2b_path)

    # ------------------------------------------------------------------ rtp
    def _read_rtp(self, path: str) -> None:
        cur: Optional[RtpEntry] = None
        sec = ""
        for raw in open(path, errors="replace"):
            line = raw.split(";")[0].strip()
            if not line:
                continue
            m = re.fullmatch(r"\[\s*(\S+)\s*\]", line)
            if m:
                key = m.group(1)
                if key in SECTION_KEYWORDS:
                    sec = key
                else:
                    cur = RtpEntry(name=key)
                    self.rtp[key] = cur
                    sec = ""
                continue
            if cur is None or sec == "":
                continue
            parts = line.split()
            if sec == "atoms" and len(parts) >= 3:
                try:
                    cgnr = int(parts[3]) if len(parts) > 3 else len(cur.atoms) + 1
                    cur.atoms.append((parts[0], parts[1], float(parts[2]), cgnr))
                except ValueError:
                    pass
            elif sec == "bonds" and len(parts) >= 2:
                cur.bonds.append((parts[0], parts[1]))
            elif sec == "impropers" and len(parts) >= 4:
                cur.impropers.append(tuple(parts[:4]))

    def _read_r2b(self, path: str) -> None:
        for raw in open(path, errors="replace"):
            line = raw.split(";")[0].strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            gmx, main, nter, cter = parts[0], parts[1], parts[2], parts[3]
            self.r2b[gmx] = (main, nter, cter)
            if main != "-":
                self.ff2gmx.setdefault(main, gmx)

    # --------------------------------------------------------------- lookup
    def block_for(self, resname: str, position: str = "middle") -> Optional[RtpEntry]:
        """position: 'middle' | 'nter' | 'cter'."""
        gmx = self.ff2gmx.get(resname, resname)
        row = self.r2b.get(gmx)
        candidates: List[str] = []
        if position == "nter":
            if row and row[1] != "-":
                candidates.append(row[1])
            candidates += ["N" + resname, resname]
        elif position == "cter":
            if row and row[2] != "-":
                candidates.append(row[2])
            candidates += ["C" + resname, resname]
        else:
            if row and row[0] != "-":
                candidates.append(row[0])
            candidates.append(resname)
        for name in candidates:
            if name in self.rtp:
                return self.rtp[name]
        return None


# Pairs "protonated form -> deprotonated form": the charge of the second one
# must be exactly one lower. Some amber ports get this wrong (TYN, for example).
PROTONATION_PAIRS = [
    ("ASH", "ASP"), ("GLH", "GLU"), ("LYS", "LYN"), ("CYS", "CYM"),
    ("TYR", "TYN"), ("ARG", "ARN"), ("HIP", "HID"), ("HIP", "HIE"),
]


def check_pair_charges(ff: "ForceField", used: set) -> List[Tuple[str, float, str]]:
    """Find blocks with a broken (fractional) charge among those actually used.

    Fixes nothing: the force field directory is opened read-only.

    :returns: list of (block, discrepancy, reference block)
    """
    out = []
    for prot, deprot in PROTONATION_PAIRS:
        if deprot not in used:
            continue
        a, b = ff.rtp.get(prot), ff.rtp.get(deprot)
        if a is None or b is None:
            continue
        delta = (a.charge - 1.0) - b.charge
        if abs(delta) > 1e-3:
            out.append((deprot, delta, prot))
    return out


def _num_key(name: str):
    m = re.search(r"(\d+)$", name)
    return (0, int(m.group(1)), name) if m else (1, 0, name)


def reconcile_residue(
    struct_names: List[str],
    parents: Dict[str, str],
    block: RtpEntry,
) -> Tuple[Dict[str, str], List[str], List[str]]:
    """Match the atom names of a residue against the names in the rtp block.

    :param struct_names: atom names in the structure (in file order)
    :param parents: hydrogen -> heavy atom (taken from the structure geometry)
    :returns: (renames, extra atoms, missing rtp atoms)
    """
    rtp_names = block.atom_names
    rename: Dict[str, str] = {}
    remaining_struct = list(struct_names)
    remaining_rtp = [n for n in rtp_names if n not in ("-C", "+N")]

    def take(src: str, dst: str) -> None:
        rename[src] = dst
        if src in remaining_struct:
            remaining_struct.remove(src)
        if dst in remaining_rtp:
            remaining_rtp.remove(dst)

    # 1. explicit rules (OXT -> OC1 and friends)
    for src, dst in EXPLICIT_RENAME["*"].items():
        if src in remaining_struct and dst in remaining_rtp:
            take(src, dst)

    # 2. names that already match
    for name in list(remaining_struct):
        if name in remaining_rtp:
            take(name, name)

    # 3. heavy atoms: if exactly one candidate is left on each side with the
    #    same leading character (element), pair them up (ILE CD1 -> CD,
    #    C-terminal O -> OC2)
    left_heavy = [n for n in remaining_struct if not _is_h(n)]
    right_heavy = [n for n in remaining_rtp if not _is_h(n)]
    for elem in {n[0] for n in left_heavy}:
        lhs = sorted([n for n in left_heavy if n[0] == elem], key=_num_key)
        rhs = sorted([n for n in right_heavy if n[0] == elem], key=_num_key)
        if lhs and len(lhs) == len(rhs):
            for src, dst in zip(lhs, rhs):
                take(src, dst)

    # 4. hydrogens: group by their heavy parent and hand out names in order
    rtp_parents = block.parents()
    by_parent_struct: Dict[str, List[str]] = {}
    for name in remaining_struct:
        if _is_h(name):
            par = parents.get(name, "?")
            by_parent_struct.setdefault(rename.get(par, par), []).append(name)
    by_parent_rtp: Dict[str, List[str]] = {}
    for name in remaining_rtp:
        if _is_h(name):
            by_parent_rtp.setdefault(rtp_parents.get(name, "?"), []).append(name)
    for par, lhs in by_parent_struct.items():
        rhs = by_parent_rtp.get(par, [])
        lhs = sorted(lhs, key=_num_key)
        rhs = sorted(rhs, key=_num_key)
        if len(lhs) == len(rhs):
            for src, dst in zip(lhs, rhs):
                take(src, dst)

    extra = [n for n in remaining_struct]
    missing = [n for n in remaining_rtp]
    return rename, extra, missing
