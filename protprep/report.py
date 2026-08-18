"""Writing the protonation reports (TSV + JSON)."""

from __future__ import annotations

import json
import os
from typing import TYPE_CHECKING, List, Optional

if TYPE_CHECKING:  # pragma: no cover - typing only
    from .pipeline import Result
    from .spec import Spec

TSV_NAME = "protonation_report.tsv"
JSON_NAME = "protonation_report.json"

COLUMNS = ["chain", "resid", "icode", "input", "final", "pKa", "model_pKa",
           "buried", "fixed", "notes"]


def _interesting(r) -> bool:
    """Skip residues that were never titrated and never touched."""
    return not (r.pka is None and not r.forced and r.original == r.final)


def write_tsv(result: "Result", path: str) -> str:
    with open(path, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for r in result.residues:
            if not _interesting(r):
                continue
            pka = f"{r.pka:.2f}" if r.pka is not None else ""
            mpka = f"{r.model_pka:.2f}" if r.model_pka is not None else ""
            bur = f"{r.buried:.2f}" if isinstance(r.buried, (int, float)) else ""
            fh.write(
                f"{r.chain}\t{r.resid}\t{r.icode}\t{r.original}\t{r.final}\t"
                f"{pka}\t{mpka}\t{bur}\t{'yes' if r.forced else ''}\t"
                f"{'; '.join(r.notes)}\n"
            )
    return path


def write_json(result: "Result", path: str, spec: Optional["Spec"] = None) -> str:
    spec = spec or result.spec
    with open(path, "w") as fh:
        json.dump(
            {
                "ph": getattr(spec, "ph", None),
                "pka_source": getattr(spec, "pka_source", None),
                "hydrogens": getattr(spec, "hydrogens", None),
                "output_pdb": result.output_pdb,
                "force_field": result.ff_dir,
                "pdb2gmx": result.pdb2gmx_cmd,
                "termini": result.termini,
                "warnings": result.warnings,
                "residues": result.records(),
            },
            fh, ensure_ascii=False, indent=1,
        )
    return path


def write_reports(result: "Result", outdir: str,
                  spec: Optional["Spec"] = None) -> List[str]:
    """Write both reports into `outdir`, returning the paths."""
    os.makedirs(outdir, exist_ok=True)
    return [
        write_tsv(result, os.path.join(outdir, TSV_NAME)),
        write_json(result, os.path.join(outdir, JSON_NAME), spec),
    ]
