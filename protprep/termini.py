"""Термини: что доступно в силовом поле и как это называется.

Важно: скрипт **никогда** ничего не пишет в каталог силового поля. Если нужного
блока в нём нет - выдаётся ошибка, а не самодельный блок с придуманными
зарядами.

В amber-портах GROMACS файлы aminoacids.n.tdb / aminoacids.c.tdb пустые: тип
конца задаётся не модификацией, а отдельным блоком rtp (NALA/CALA). Нейтральных
вариантов (NH2 и COOH) там обычно просто нет. Если вы добавили их в своё
силовое поле сами, назовите блоки по этой схеме - и они будут подхвачены:

    ZXXX - нейтральный N-конец  (NH2),  например ZALA
    JXXX - нейтральный C-конец (COOH),  например JALA
"""

from __future__ import annotations

import os
from typing import List, Optional, Set, Tuple

from .ffdata import ForceField

NEUTRAL_NTER_PREFIX = "Z"
NEUTRAL_CTER_PREFIX = "J"


def neutral_nter_name(resname: str) -> str:
    return NEUTRAL_NTER_PREFIX + resname[:3]


def neutral_cter_name(resname: str) -> str:
    return NEUTRAL_CTER_PREFIX + resname[:3]


def system_residuetypes() -> Tuple[Optional[str], Set[str]]:
    """Путь к residuetypes.dat установленного GROMACS и известные ему имена."""
    candidates: List[str] = []
    gmxdata = os.environ.get("GMXDATA")
    if gmxdata:
        candidates.append(os.path.join(gmxdata, "top", "residuetypes.dat"))
    candidates += [
        "/usr/local/gromacs/share/gromacs/top/residuetypes.dat",
        "/usr/share/gromacs/top/residuetypes.dat",
        "/usr/local/share/gromacs/top/residuetypes.dat",
    ]
    for cand in candidates:
        if os.path.exists(cand):
            names = {
                line.split()[0]
                for line in open(cand, errors="replace")
                if line.split() and not line.startswith(";")
            }
            return cand, names
    return None, set()


def check_neutral_block(ff: ForceField, resname: str, end: str) -> Optional[str]:
    """Имя блока нейтрального конца, если он есть в силовом поле."""
    name = neutral_nter_name(resname) if end == "N" else neutral_cter_name(resname)
    return name if name in ff.rtp else None
