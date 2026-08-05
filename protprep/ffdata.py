"""Чтение данных силового поля GROMACS (*.rtp, *.r2b) и приведение имён атомов.

Идея: pdb2pqr выдаёт структуру в номенклатуре AMBER (HB2/HB3, CD1 у ILE, OXT...),
а gromacs-порт amber99sb-ildn ждёт свою (HB1/HB2, CD, OC1/OC2...). Вместо
захардкоженной таблицы мы читаем сам rtp из папки силового поля и сопоставляем
имена автоматически - тогда скрипт переживёт любую правку ff.
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

# Явные соответствия, которые нельзя вывести автоматически (симметричные пары).
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
        """Для каждого водорода - тяжёлый атом, к которому он привязан."""
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
    """Каталоги top установленного GROMACS, где лежат *.ff."""
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
    """Ищет каталог силового поля: как указано, рядом, затем в GROMACS.

    Принимает и путь, и просто имя - с суффиксом .ff и без него.
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
        hint = ("\nВ установленном GROMACS есть: " + ", ".join(available)
                + "\nУкажите нужное через --ff (можно просто именем).")
    raise FileNotFoundError(f"Каталог силового поля не найден: {name}{hint}")


class ForceField:
    """Ленивая обёртка над каталогом <name>.ff."""

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

    # ------------------------------------------------------------- запросы
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


# Пары "протонированная форма -> депротонированная": заряд второй обязан быть
# ровно на 1 меньше. В некоторых портах amber это нарушено (например TYN).
PROTONATION_PAIRS = [
    ("ASH", "ASP"), ("GLH", "GLU"), ("LYS", "LYN"), ("CYS", "CYM"),
    ("TYR", "TYN"), ("ARG", "ARN"), ("HIP", "HID"), ("HIP", "HIE"),
]


def check_pair_charges(ff: "ForceField", used: set) -> List[Tuple[str, float, str]]:
    """Ищет блоки с неверным (дробным) зарядом среди реально использованных.

    Ничего не чинит: каталог силового поля открывается только на чтение.

    :returns: список (блок, расхождение, эталонный блок)
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
    """Сопоставляет имена атомов остатка с именами из rtp.

    :param struct_names: имена атомов в структуре (в порядке файла)
    :param parents: водород -> тяжёлый атом (из геометрии структуры)
    :returns: (переименования, лишние атомы, недостающие атомы rtp)
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

    # 1. явные правила (OXT -> OC1 и т.п.)
    for src, dst in EXPLICIT_RENAME["*"].items():
        if src in remaining_struct and dst in remaining_rtp:
            take(src, dst)

    # 2. совпадающие имена
    for name in list(remaining_struct):
        if name in remaining_rtp:
            take(name, name)

    # 3. тяжёлые атомы: если осталось ровно по одному кандидату с одинаковым
    #    первым символом (элементом) - связываем (ILE CD1 -> CD, C-конец O -> OC2)
    left_heavy = [n for n in remaining_struct if not _is_h(n)]
    right_heavy = [n for n in remaining_rtp if not _is_h(n)]
    for elem in {n[0] for n in left_heavy}:
        lhs = sorted([n for n in left_heavy if n[0] == elem], key=_num_key)
        rhs = sorted([n for n in right_heavy if n[0] == elem], key=_num_key)
        if lhs and len(lhs) == len(rhs):
            for src, dst in zip(lhs, rhs):
                take(src, dst)

    # 4. водороды: группируем по тяжёлому "родителю" и раздаём имена по порядку
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
