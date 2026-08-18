"""Running gmx pdb2gmx - shared by the CLI and the Python API."""

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from typing import List, Optional


@dataclass
class Pdb2gmxResult:
    """Outcome of a gmx pdb2gmx run."""

    returncode: int
    command: List[str]
    workdir: str
    log_path: Optional[str]
    stdout: str = ""
    stderr: str = ""
    gro: Optional[str] = None
    top: Optional[str] = None

    @property
    def ok(self) -> bool:
        return self.returncode == 0

    @property
    def total_charge(self) -> Optional[float]:
        """Total system charge as reported by pdb2gmx, if it said anything."""
        value = None
        for line in (self.stderr + "\n" + self.stdout).splitlines():
            if "Total charge" in line:
                for token in line.replace("e", " ").split():
                    try:
                        value = float(token)
                    except ValueError:
                        continue
        return value

    def charge_lines(self) -> List[str]:
        return [l.strip() for l in (self.stderr + "\n" + self.stdout).splitlines()
                if "Total charge" in l]

    def tail(self, n: int = 25) -> List[str]:
        text = (self.stderr or self.stdout).strip()
        return text.splitlines()[-n:] if text else []


def find_gmx(name: str = "gmx") -> Optional[str]:
    """gmx is often not on PATH (GMXRC needed) - look in the usual places."""
    found = shutil.which(name)
    if found:
        return found
    for cand in (
        "/usr/local/gromacs/bin/gmx", "/usr/local/gromacs/bin/gmx_mpi",
        "/usr/bin/gmx", "/opt/gromacs/bin/gmx",
    ):
        if os.path.exists(cand):
            return cand
    return None


def run_pdb2gmx(
    pdb: str,
    ff_dir: str,
    outdir: str,
    water: str = "tip3p",
    gmx: str = "gmx",
    extra_args: Optional[List[str]] = None,
) -> Pdb2gmxResult:
    """Build a topology from an already prepared PDB.

    pdb2gmx looks for the force field in the current directory, so the command
    runs from the directory that holds the *.ff folder.
    """
    exe = find_gmx(gmx)
    if exe is None:
        raise FileNotFoundError(
            f"cannot find the executable '{gmx}' (source GMXRC first)"
        )
    ff_dir = os.path.abspath(ff_dir)
    workdir = os.path.dirname(ff_dir)
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)
    gro = os.path.join(outdir, "conf.gro")
    top = os.path.join(outdir, "topol.top")
    cmd = [
        exe, "pdb2gmx", "-f", os.path.abspath(pdb),
        "-o", gro, "-p", top, "-i", os.path.join(outdir, "posre.itp"),
        "-ff", os.path.basename(ff_dir)[:-3], "-water", water,
    ] + list(extra_args or [])
    proc = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
    log_path = os.path.join(outdir, "pdb2gmx.log")
    with open(log_path, "w") as fh:
        fh.write(proc.stdout + "\n" + proc.stderr)
    return Pdb2gmxResult(
        returncode=proc.returncode, command=cmd, workdir=workdir,
        log_path=log_path, stdout=proc.stdout, stderr=proc.stderr,
        gro=gro if proc.returncode == 0 else None,
        top=top if proc.returncode == 0 else None,
    )
