"""Локальные операции с протонами: постановка/снятие титруемых H и кэпирование.

Всё, что pdb2pqr отказывается делать для силового поля AMBER (например,
протонировать ASP на конце цепи или депротонировать TYR), доделываем здесь
геометрически.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .geometry import best_torsion, dihedral, place_atom
from .pdbio import Atom

# ---------------------------------------------------------------- титруемые H
# имя H -> (ref, middle, root, длина связи, угол, торсион | None = подобрать)
PROTON_GEOMETRY: Dict[str, Dict[str, tuple]] = {
    "ASP": {"HD2": ("OD1", "CG", "OD2", 0.98, 113.0, 0.0)},
    "GLU": {"HE2": ("OE1", "CD", "OE2", 0.98, 113.0, 0.0)},
    "LYS": {"HZ3": ("CD", "CE", "NZ", 1.01, 109.5, None)},
    "CYS": {"HG": ("CA", "CB", "SG", 1.34, 96.0, None)},
    "TYR": {"HH": ("CE1", "CZ", "OH", 0.97, 108.0, 180.0)},
    "ARG": {"HH22": ("NE", "CZ", "NH2", 1.01, 120.0, 180.0)},
    "HIS": {
        "HD1": ("CD2", "CG", "ND1", 1.01, 125.6, 180.0),
        "HE2": ("CG", "CD2", "NE2", 1.01, 125.6, 180.0),
    },
}

# семейство -> {имя состояния: набор титруемых протонов, которые должны быть}
STATE_PROTONS: Dict[str, Dict[str, frozenset]] = {
    "ASP": {"ASP": frozenset(), "ASH": frozenset({"HD2"})},
    "GLU": {"GLU": frozenset(), "GLH": frozenset({"HE2"})},
    "LYS": {"LYS": frozenset({"HZ3"}), "LYN": frozenset()},
    "CYS": {"CYS": frozenset({"HG"}), "CYM": frozenset(), "CYX": frozenset()},
    "TYR": {"TYR": frozenset({"HH"}), "TYN": frozenset()},
    "ARG": {"ARG": frozenset({"HH22"}), "ARN": frozenset()},
    "HIS": {
        "HIP": frozenset({"HD1", "HE2"}),
        "HID": frozenset({"HD1"}),
        "HIE": frozenset({"HE2"}),
    },
}

FAMILY_OF: Dict[str, str] = {
    name: family for family, table in STATE_PROTONS.items() for name in table
}
FAMILY_OF.update({"HIS": "HIS", "HSD": "HIS", "HSE": "HIS", "HSP": "HIS"})


def family(resname: str) -> Optional[str]:
    return FAMILY_OF.get(resname.upper())


def _by_name(atoms: Sequence[Atom]) -> Dict[str, Atom]:
    return {a.name: a for a in atoms}


def _template(res: Sequence[Atom], name: str, xyz) -> Atom:
    """Новый атом-водород по образцу остатка."""
    ref = res[0]
    return Atom(
        record=ref.record, serial=0, name=name, altloc=" ", resname=ref.resname,
        chain=ref.chain, resseq=ref.resseq, icode=ref.icode,
        x=float(xyz[0]), y=float(xyz[1]), z=float(xyz[2]),
        occ=1.0, bfac=0.0, element="H",
    )


def add_proton(res: List[Atom], hname: str, fam: str,
               cloud: Optional[np.ndarray] = None) -> Tuple[bool, str]:
    """Ставит титруемый протон hname. Возвращает (успех, сообщение)."""
    geom = PROTON_GEOMETRY[fam].get(hname)
    if geom is None:
        return False, f"нет геометрии для {hname}"
    ref, mid, root, bond, ang, tors = geom
    idx = _by_name(res)
    for need in (ref, mid, root):
        if need not in idx:
            return False, f"в структуре нет тяжёлого атома {need}"
    a, b, c = idx[ref], idx[mid], idx[root]

    if tors is None:
        # ищем свободный поворот с учётом уже стоящих на root водородов
        siblings = [
            at for at in res
            if at.element == "H" and at.name != hname
            and np.linalg.norm(np.array([at.x - c.x, at.y - c.y, at.z - c.z])) < 1.4
        ]
        if siblings:
            base = dihedral(a, b, c, siblings[0])
            best, best_gap = base + 120.0, -1.0
            for cand in (base + 120.0, base + 240.0):
                gap = min(
                    abs(((cand - dihedral(a, b, c, s) + 180) % 360) - 180)
                    for s in siblings
                )
                if gap > best_gap:
                    best, best_gap = cand, gap
            tors = best
        else:
            tors, _ = best_torsion(a, b, c, bond, ang,
                                   cloud if cloud is not None else np.empty((0, 3)),
                                   preferred=180.0)
    pos = place_atom(a, b, c, bond, ang, tors)
    res.append(_template(res, hname, pos))
    return True, ""


def set_state(res: List[Atom], target: str,
              cloud: Optional[np.ndarray] = None) -> Tuple[List[Atom], List[str]]:
    """Приводит остаток к состоянию target (ASH/HID/LYN/...). Возвращает
    (новый список атомов, список замечаний)."""
    notes: List[str] = []
    fam = family(target)
    if fam is None:
        return res, [f"состояние {target} не поддерживается"]
    wanted = STATE_PROTONS[fam][target]
    titratable = set(PROTON_GEOMETRY[fam])

    out = [a for a in res if not (a.name in titratable and a.name not in wanted)]
    removed = {a.name for a in res} - {a.name for a in out}
    for name in sorted(removed):
        notes.append(f"снят протон {name}")

    have = {a.name for a in out}
    for name in sorted(wanted - have):
        ok, msg = add_proton(out, name, fam, cloud)
        notes.append(f"добавлен протон {name}" if ok else f"НЕ добавлен {name}: {msg}")

    out = [replace(a, resname=target) for a in out]
    return out, notes


# ------------------------------------------------------------------ кэпы
def _mk(res_ref: Atom, name: str, resname: str, resseq: int, xyz,
        element: str) -> Atom:
    return Atom(
        record="ATOM", serial=0, name=name, altloc=" ", resname=resname,
        chain=res_ref.chain, resseq=resseq, icode="", x=float(xyz[0]),
        y=float(xyz[1]), z=float(xyz[2]), occ=1.0, bfac=0.0, element=element,
    )


def _methyl(root, neighbor, ref, res_ref, resname, resseq, names, cloud):
    """Три водорода метила root; поворот подбираем по минимуму контактов."""
    tors, _ = best_torsion(ref, neighbor, root, 1.09, 109.5, cloud, preferred=60.0,
                           step=20.0)
    out = []
    for i, nm in enumerate(names):
        pos = place_atom(ref, neighbor, root, 1.09, 109.5, tors + 120.0 * i)
        out.append(_mk(res_ref, nm, resname, resseq, pos, "H"))
    return out


def cap_nterm_ace(res: List[Atom], cloud: np.ndarray) -> Tuple[List[Atom], List[Atom], List[str]]:
    """ACE-кэп перед первым остатком. Возвращает (атомы ACE, изменённый остаток, заметки)."""
    idx = _by_name(res)
    notes: List[str] = []
    for need in ("N", "CA", "C"):
        if need not in idx:
            return [], res, [f"нет атома {need}, ACE не поставлен"]
    n, ca, c = idx["N"], idx["CA"], idx["C"]

    # убираем "лишние" протоны NH3+
    body = [a for a in res if a.name not in ("H1", "H2", "H3", "HN1", "HN2", "HN3")]
    resseq = res[0].resseq - 1

    phi, _ = best_torsion(c, ca, n, 1.335, 121.7, cloud, preferred=-60.0, step=15.0)
    c_ace = place_atom(c, ca, n, 1.335, 121.7, phi)
    ace_c = _mk(res[0], "C", "ACE", resseq, c_ace, "C")
    o_pos = place_atom(ca, n, ace_c, 1.229, 122.9, 0.0)
    ace_o = _mk(res[0], "O", "ACE", resseq, o_pos, "O")
    ch3_pos = place_atom(ca, n, ace_c, 1.508, 116.6, 180.0)
    ace_ch3 = _mk(res[0], "CH3", "ACE", resseq, ch3_pos, "C")
    hs = _methyl(ace_ch3, ace_c, ace_o, res[0], "ACE", resseq,
                 ["HH31", "HH32", "HH33"], cloud)

    # амидный H на N (для PRO его нет)
    if res[0].resname != "PRO":
        h_pos = place_atom(ace_o, ace_c, n, 1.01, 119.8, 180.0)
        body.append(_template(body, "H", h_pos))
        notes.append("N-конец: NH3+ -> амидный H + ACE")
    else:
        notes.append("N-конец: PRO кэпирован ACE (протоны на N убраны)")

    ace = [hs[0], ace_ch3, hs[1], hs[2], ace_c, ace_o]  # порядок как в rtp
    return ace, body, notes


def cap_cterm(res: List[Atom], cloud: np.ndarray,
              kind: str = "NME") -> Tuple[List[Atom], List[Atom], List[str]]:
    """NME/NHE-кэп после последнего остатка."""
    idx = _by_name(res)
    notes: List[str] = []
    for need in ("N", "CA", "C", "O"):
        if need not in idx:
            return [], res, [f"нет атома {need}, {kind} не поставлен"]
    n, ca, c, o = idx["N"], idx["CA"], idx["C"], idx["O"]

    body = [a for a in res if a.name not in ("OXT", "HXT", "OT2", "HO")]
    resseq = res[0].resseq + 1

    psi = dihedral(n, ca, c, o) + 180.0
    n_pos = place_atom(n, ca, c, 1.335, 116.6, psi)
    cap_n = _mk(res[0], "N", kind, resseq, n_pos, "N")

    if kind == "NME":
        h_pos = place_atom(o, c, cap_n, 1.01, 119.8, 180.0)
        cap_h = _mk(res[0], "H", kind, resseq, h_pos, "H")
        ch3_pos = place_atom(o, c, cap_n, 1.449, 121.9, 0.0)
        cap_ch3 = _mk(res[0], "CH3", kind, resseq, ch3_pos, "C")
        hs = _methyl(cap_ch3, cap_n, c, res[0], kind, resseq,
                     ["HH31", "HH32", "HH33"], cloud)
        cap = [cap_n, cap_h, cap_ch3] + hs
        notes.append("C-конец: COO- -> NME")
    else:  # NHE - амид NH2
        h1 = place_atom(o, c, cap_n, 1.01, 119.8, 180.0)
        h2 = place_atom(o, c, cap_n, 1.01, 119.8, 0.0)
        cap = [
            cap_n,
            _mk(res[0], "H1", kind, resseq, h1, "H"),
            _mk(res[0], "H2", kind, resseq, h2, "H"),
        ]
        notes.append("C-конец: COO- -> NHE (амид)")
    return cap, body, notes


def coords(atoms: Sequence[Atom]) -> np.ndarray:
    if not atoms:
        return np.empty((0, 3))
    return np.array([[a.x, a.y, a.z] for a in atoms])
