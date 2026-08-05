"""Разбор пользовательского задания: какие остатки зафиксировать и как."""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# --- целевые имена остатков в номенклатуре amber/gromacs-порта ------------
PROTONATED = "protonated"
DEPROTONATED = "deprotonated"

# состояние -> (имя остатка, ветка титрования)
# ветка нужна, чтобы подсунуть pdb2pqr фиктивное pKa: pH < pKa  => протонированная
# форма, pH >= pKa => депротонированная.
STATE_TABLE: Dict[str, Dict[str, str]] = {
    "ASP": {"ASP": DEPROTONATED, "ASH": PROTONATED},
    "GLU": {"GLU": DEPROTONATED, "GLH": PROTONATED},
    "LYS": {"LYS": PROTONATED, "LYN": DEPROTONATED},
    "CYS": {"CYS": PROTONATED, "CYM": DEPROTONATED, "CYX": PROTONATED},
    "TYR": {"TYR": PROTONATED, "TYN": DEPROTONATED},
    "ARG": {"ARG": PROTONATED, "ARN": DEPROTONATED},
    "HIS": {"HIP": PROTONATED, "HID": DEPROTONATED, "HIE": DEPROTONATED},
}

# синонимы -> целевое имя, отдельно для каждого типа остатка.
# Однобуквенные сокращения (p/d/n/c) добавляются ниже автоматически.
ALIASES: Dict[str, Dict[str, str]] = {
    "ASP": {"protonated": "ASH", "neutral": "ASH", "deprotonated": "ASP",
            "charged": "ASP", "ash": "ASH", "asph": "ASH", "asp": "ASP"},
    "GLU": {"protonated": "GLH", "neutral": "GLH", "deprotonated": "GLU",
            "charged": "GLU", "glh": "GLH", "gluh": "GLH", "glu": "GLU"},
    "LYS": {"protonated": "LYS", "charged": "LYS", "deprotonated": "LYN",
            "neutral": "LYN", "lyn": "LYN", "lysn": "LYN", "lys": "LYS"},
    "CYS": {"protonated": "CYS", "neutral": "CYS", "deprotonated": "CYM",
            "charged": "CYM", "cym": "CYM", "cys": "CYS", "cyx": "CYX",
            "ss": "CYX", "disulfide": "CYX"},
    "TYR": {"protonated": "TYR", "neutral": "TYR", "deprotonated": "TYN",
            "charged": "TYN", "tyn": "TYN", "tym": "TYN", "tyr": "TYR"},
    "ARG": {"protonated": "ARG", "charged": "ARG", "deprotonated": "ARN",
            "neutral": "ARN", "arn": "ARN", "ar0": "ARN", "arg": "ARG"},
    "HIS": {"protonated": "HIP", "charged": "HIP", "deprotonated": "HIE",
            "neutral": "HIE", "hip": "HIP", "hid": "HID", "hie": "HIE",
            "hish": "HIP", "hisd": "HID", "hise": "HIE", "d": "HID",
            "e": "HIE", "his": "HIE"},
}

# p/d/n/c - короткая запись для protonated/deprotonated/neutral/charged
SHORT = {"p": "protonated", "d": "deprotonated", "n": "neutral", "c": "charged"}
for _table in ALIASES.values():
    for _short, _long in SHORT.items():
        _table[_short] = _table[_long]
# у гистидина 'd'/'e' уже заняты тавтомерами HID/HIE - это важнее
ALIASES["HIS"]["d"] = "HID"
ALIASES["HIS"]["e"] = "HIE"

TITRATABLE = set(STATE_TABLE)

# --- термини --------------------------------------------------------------
NTER_STATES = {
    "NH3+": "NH3+", "CHARGED": "NH3+", "NH3": "NH3+", "PROTONATED": "NH3+",
    "P": "NH3+", "C": "NH3+",
    "NH2": "NH2", "NEUTRAL": "NH2", "DEPROTONATED": "NH2",
    "N": "NH2", "D": "NH2",
    "ACE": "ACE", "CAP": "ACE", "ACETYL": "ACE",
    "KEEP": "KEEP", "NONE": "KEEP",
}
CTER_STATES = {
    "COO-": "COO-", "COO": "COO-", "CHARGED": "COO-", "DEPROTONATED": "COO-",
    "D": "COO-", "C": "COO-",
    "COOH": "COOH", "NEUTRAL": "COOH", "PROTONATED": "COOH",
    "P": "COOH", "N": "COOH",
    "NME": "NME", "CAP": "NME", "NMETHYL": "NME",
    "NHE": "NHE", "NH2": "NHE", "AMIDE": "NHE",
    "KEEP": "KEEP", "NONE": "KEEP",
}


@dataclass(frozen=True)
class ResidueSpec:
    chain: str
    resid: int
    icode: str
    state: str            # уже нормализованное имя остатка (ASH/HID/...)
    raw: str = ""

    @property
    def key(self) -> Tuple[str, int, str]:
        return (self.chain, self.resid, self.icode)


@dataclass(frozen=True)
class TerminusSpec:
    chain: str
    end: str              # 'N' или 'C'
    state: str


@dataclass
class Spec:
    ph: float = 7.0
    # 'propka' - считать локальные сдвиги pKa; 'standard' - взять табличные
    pka_source: str = "propka"
    residues: List[ResidueSpec] = field(default_factory=list)
    termini: List[TerminusSpec] = field(default_factory=list)

    def residue_map(self) -> Dict[Tuple[str, int, str], ResidueSpec]:
        return {r.key: r for r in self.residues}

    def terminus(self, chain: str, end: str) -> Optional[str]:
        for t in self.termini:
            if t.chain in (chain, "*") and t.end == end:
                return t.state
        return None


class SpecError(ValueError):
    pass


_RES_RE = re.compile(
    r"^\s*(?P<chain>[A-Za-z0-9])?\s*[:/ ]?\s*(?P<resid>-?\d+)(?P<icode>[A-Za-z])?"
    r"\s*[:= ]\s*(?P<state>\S+)\s*$"
)


def _normalize_state(parent: str, state: str) -> str:
    parent = parent.upper()
    key = state.strip().lower()
    table = ALIASES.get(parent)
    if table is None:
        raise SpecError(
            f"Остаток {parent} не титруется этим скриптом "
            f"(поддерживаются: {', '.join(sorted(TITRATABLE))})"
        )
    if key not in table:
        raise SpecError(
            f"Неизвестное состояние '{state}' для {parent}. "
            f"Доступно: {', '.join(sorted(set(table.values()) | {'protonated', 'deprotonated'}))}"
        )
    return table[key]


def normalize_residue_state(parent_resname: str, state: str) -> str:
    return _normalize_state(parent_resname, state)


def parse_terminus_token(end: str, value: str) -> str:
    table = NTER_STATES if end == "N" else CTER_STATES
    key = value.strip().upper()
    if key not in table:
        raise SpecError(
            f"Неизвестное состояние {end}-конца: '{value}'. "
            f"Доступно: {', '.join(sorted(set(table.values())))}"
        )
    return table[key]


def parse_pka_source(value: str) -> str:
    """'propka' | 'standard' (синонимы: model, table, none, no)."""
    key = value.strip().lower()
    if key in ("propka", "local", "yes"):
        return "propka"
    if key in ("standard", "model", "table", "none", "no", "off"):
        return "standard"
    raise SpecError(
        f"Неизвестный источник pKa: '{value}'. Доступно: propka | standard"
    )


def parse_fix_token(token: str) -> Tuple[str, int, str, str]:
    """'A:145:ASH' | 'A/145/protonated' | 'A 145 ASH' -> (chain, resid, icode, state)"""
    parts = [p for p in re.split(r"[:/,]|\s+", token.strip()) if p]
    if len(parts) == 3:
        chain, resid, state = parts
    elif len(parts) == 2:
        chain, resid, state = "", parts[0], parts[1]
    else:
        raise SpecError(f"Не могу разобрать '{token}'. Формат: ЦЕПЬ:НОМЕР:СОСТОЯНИЕ")
    m = re.fullmatch(r"(-?\d+)([A-Za-z]?)", resid)
    if not m:
        raise SpecError(f"Плохой номер остатка в '{token}'")
    return chain, int(m.group(1)), m.group(2).upper(), state


def load_spec(path: str) -> Spec:
    """JSON или простой текстовый формат."""
    text = open(path, errors="replace").read()
    if os.path.splitext(path)[1].lower() == ".json" or text.lstrip().startswith("{"):
        return _load_json(text)
    return _load_text(text)


def _load_json(text: str) -> Spec:
    data = json.loads(text)
    spec = Spec(ph=float(data.get("ph", 7.0)),
                pka_source=parse_pka_source(str(data.get("pka", "propka"))))
    for item in data.get("residues", []):
        spec.residues.append(
            ResidueSpec(
                chain=str(item.get("chain", "")),
                resid=int(item["resid"]),
                icode=str(item.get("icode", "")).upper(),
                state=str(item["state"]),
                raw=str(item["state"]),
            )
        )
    for item in data.get("termini", []):
        end = str(item["end"]).upper()[0]
        spec.termini.append(
            TerminusSpec(
                chain=str(item.get("chain", "*")),
                end=end,
                state=parse_terminus_token(end, str(item["state"])),
            )
        )
    return spec


def _load_text(text: str) -> Spec:
    spec = Spec()
    for lineno, raw in enumerate(text.splitlines(), 1):
        line = raw.split("#")[0].split(";")[0].strip()
        if not line:
            continue
        parts = [p for p in re.split(r"[:=,]|\s+", line) if p]
        head = parts[0].lower()
        try:
            if head == "ph":
                spec.ph = float(parts[1])
            elif head == "pka":
                spec.pka_source = parse_pka_source(parts[1])
            elif head in ("nter", "n-term", "nterm", "cter", "c-term", "cterm"):
                end = "N" if head.startswith("n") else "C"
                chain = parts[1] if len(parts) == 3 else "*"
                spec.termini.append(
                    TerminusSpec(chain, end, parse_terminus_token(end, parts[-1]))
                )
            elif len(parts) >= 3 and parts[1].lstrip("-").isdigit():
                chain, resid, icode, state = parse_fix_token(" ".join(parts[:3]))
                spec.residues.append(ResidueSpec(chain, resid, icode, state, state))
            elif len(parts) >= 2 and parts[0].lstrip("-").isdigit():
                chain, resid, icode, state = parse_fix_token(" ".join(parts[:2]))
                spec.residues.append(ResidueSpec(chain, resid, icode, state, state))
            else:
                raise SpecError(f"не понял строку: '{raw.strip()}'")
        except (IndexError, ValueError) as err:
            raise SpecError(f"Строка {lineno}: {err}")
    return spec
