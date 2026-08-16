"""Termini: what the force field offers and how it is named.

Important: this tool **never** writes anything into the force field directory.
If a required building block is missing, it reports an error instead of
inventing a block with made-up charges.

In the GROMACS amber ports aminoacids.n.tdb / aminoacids.c.tdb are empty: the
terminus type comes not from a modification but from a separate rtp block
(NALA/CALA). Neutral variants (NH2 and COOH) are usually absent altogether.
If you added them to your own force field, name the blocks following this
scheme and they will be picked up:

    ZXXX - neutral N-terminus  (NH2),  e.g. ZALA
    JXXX - neutral C-terminus (COOH),  e.g. JALA
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
    """Path to the installed GROMACS residuetypes.dat and the names it knows."""
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
    """Name of the neutral terminus block, if the force field has one."""
    name = neutral_nter_name(resname) if end == "N" else neutral_cter_name(resname)
    return name if name in ff.rtp else None
