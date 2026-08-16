"""Command line interface."""

from __future__ import annotations

import argparse
import json
import logging
import os
import subprocess
import sys
from typing import List

from .pipeline import ForceFieldError, Result, prepare
from .spec import (Spec, SpecError, TerminusSpec, ResidueSpec, load_spec,
                   parse_fix_token, parse_terminus_token)

LOG = logging.getLogger("protprep")

EPILOG = """
examples:
  # Asp145 protonated, His264 as the HIE tautomer, everything else from PROPKA
  protonate -f 6CFO.pdb --ph 7.4 --fix A:145:p --fix A:264:HIE

  # cap both ends of chain A, leave chain B charged
  protonate -f 6CFO.pdb --nter A:ACE --cter A:NME

  # the same thing through a spec file
  protonate -f 6CFO.pdb --spec fixed.txt

spec file format (text or json):
  pH 7.4
  A 145 p               # ASP -> ASH
  A 264 HIE             # histidine tautomer
  B 87  d               # LYS -> LYN
  *  63 HIE             # in every chain at once
  nter A ACE
  cter A NME
  pka standard          # do not run PROPKA
  hydrogens fixed       # H only on the residues from --fix

pKa source:
  PROPKA (local shifts) by default. With --standard-pka PROPKA is not run at
  all: every unpinned group, N- and C-termini included, follows the standard
  pKa values of free amino acids at the requested pH.

hydrogens:
  the whole structure is protonated by default. With --only-fixed-h hydrogens
  are kept only on the residues from --fix (and on the caps); the rest leave
  without hydrogens and gmx pdb2gmx rebuilds them. Residue names, and thus the
  states, are preserved.

states:
  p = protonated, d = deprotonated, n = neutral, c = charged
  or an explicit name: ASP|ASH, GLU|GLH, LYS|LYN, CYS|CYM|CYX, TYR|TYN,
                       ARG|ARN, HID|HIE|HIP  (for HIS: d = HID, e = HIE, p = HIP)
termini:
  N -> NH3+ (p, c) | NH2 (n, d) | ACE ;  C -> COO- (d, c) | COOH (p, n) | NME | NHE

The force field is opened read-only. If a required block is missing, the tool
says so and writes nothing.
"""


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="protonate",
        description="Protonate a protein at a given pH (PROPKA) with the "
                    "states of selected residues pinned by hand.",
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-f", "--pdb", required=True, help="input PDB (without hydrogens)")
    p.add_argument("-o", "--outdir", default="prepared", help="output directory")
    p.add_argument("-n", "--name", default="protonated.pdb", help="output PDB name")
    p.add_argument("--ph", type=float, default=None,
                   help="pH of the medium (default 7.0)")
    p.add_argument("--standard-pka", "--no-propka", dest="standard_pka",
                   action="store_true",
                   help="do not run PROPKA: every unpinned group, termini "
                        "included, takes the state implied by the standard "
                        "(tabulated) pKa at the given pH")
    p.add_argument("--only-fixed-h", "--strip-other-h", dest="only_fixed_h",
                   action="store_true",
                   help="place hydrogens only on the residues from --fix (and "
                        "on the caps), leaving all the others without "
                        "hydrogens - gmx pdb2gmx will add them")
    p.add_argument("--ff", default="amber-99sb-ildn.ff",
                   help="GROMACS force field: a path or a bare name")
    p.add_argument("--spec", help="spec file (txt or json)")
    p.add_argument("--fix", action="append", default=[],
                   metavar="CHAIN:NUMBER:STATE",
                   help="pin the state of a residue (repeatable)")
    p.add_argument("--nter", action="append", default=[], metavar="CHAIN:STATE",
                   help="N-terminus state of a chain: NH3+|NH2|ACE")
    p.add_argument("--cter", action="append", default=[], metavar="CHAIN:STATE",
                   help="C-terminus state of a chain: COO-|COOH|NME|NHE")
    p.add_argument("--drop-water", action="store_true", help="discard water")
    p.add_argument("--drop-het", action="store_true",
                   help="discard ligands/ions (HETATM)")
    p.add_argument("--no-debump", action="store_true",
                   help="pdb2pqr without debumping")
    p.add_argument("--no-opt", action="store_true",
                   help="pdb2pqr without hydrogen bond optimisation")
    p.add_argument("--water", default="tip3p", help="water model for pdb2gmx")
    p.add_argument("--run-pdb2gmx", action="store_true",
                   help="check the result right away by running gmx pdb2gmx")
    p.add_argument("--gmx", default="gmx", help="gromacs executable")
    p.add_argument("-v", "--verbose", action="store_true")
    return p


def build_spec(args) -> Spec:
    spec = load_spec(args.spec) if args.spec else Spec()
    if args.ph is not None:
        spec.ph = args.ph
    if args.standard_pka:
        spec.pka_source = "standard"
    if args.only_fixed_h:
        spec.hydrogens = "fixed"
    for token in args.fix:
        chain, resid, icode, state = parse_fix_token(token)
        spec.residues.append(ResidueSpec(chain, resid, icode, state, state))
    for token in args.nter:
        chain, _, value = token.partition(":")
        if not value:
            chain, value = "*", chain
        spec.termini.append(TerminusSpec(chain, "N", parse_terminus_token("N", value)))
    for token in args.cter:
        chain, _, value = token.partition(":")
        if not value:
            chain, value = "*", chain
        spec.termini.append(TerminusSpec(chain, "C", parse_terminus_token("C", value)))
    return spec


def write_report(result: Result, spec: Spec, outdir: str) -> None:
    tsv = os.path.join(outdir, "protonation_report.tsv")
    with open(tsv, "w") as fh:
        fh.write("chain\tresid\ticode\tinput\tfinal\tpKa\tmodel_pKa\tburied\t"
                 "fixed\tnotes\n")
        for r in result.residues:
            if r.pka is None and not r.forced and r.original == r.final:
                continue
            pka = f"{r.pka:.2f}" if r.pka is not None else ""
            mpka = f"{r.model_pka:.2f}" if r.model_pka is not None else ""
            bur = f"{r.buried:.2f}" if isinstance(r.buried, (int, float)) else ""
            fh.write(
                f"{r.chain}\t{r.resid}\t{r.icode}\t{r.original}\t{r.final}\t"
                f"{pka}\t{mpka}\t{bur}\t{'yes' if r.forced else ''}\t"
                f"{'; '.join(r.notes)}\n"
            )

    js = os.path.join(outdir, "protonation_report.json")
    with open(js, "w") as fh:
        json.dump(
            {
                "ph": spec.ph,
                "pka_source": spec.pka_source,
                "hydrogens": spec.hydrogens,
                "output_pdb": result.output_pdb,
                "force_field": result.ff_dir,
                "pdb2gmx": result.pdb2gmx_cmd,
                "termini": result.termini,
                "warnings": result.warnings,
                "residues": [
                    {
                        "chain": r.chain, "resid": r.resid, "icode": r.icode,
                        "input": r.original, "final": r.final, "pKa": r.pka,
                        "model_pKa": r.model_pka, "buried": r.buried,
                        "fixed": r.forced, "notes": r.notes,
                    }
                    for r in result.residues
                ],
            },
            fh, ensure_ascii=False, indent=1,
        )


def summarize(result: Result, spec: Spec) -> None:
    changed = [r for r in result.residues if r.original != r.final]
    forced = [r for r in result.residues if r.forced]
    source = ("standard pKa (PROPKA was not run)"
              if spec.pka_source == "standard" else "local pKa from PROPKA")
    print(f"\npH = {spec.ph:.2f}; {source}")
    print(f"structure: {result.output_pdb}")
    print(f"residues with a shifted state: {len(changed)} "
          f"(pinned by hand: {len(forced)})")
    if changed:
        label = "pKa(table)" if spec.pka_source == "standard" else "pKa(PROPKA)"
        print(f"\n  chain residue    was -> now   {label}  source")
        for r in changed:
            pka = f"{r.pka:8.2f}" if r.pka is not None else "       -"
            src = ("PINNED" if r.forced
                   else ("table" if spec.pka_source == "standard" else "propka"))
            print(f"  {r.chain:>5} {r.resid:>7}   {r.original:>4} -> {r.final:<5} "
                  f"{pka}   {src}")
    missing_forced = [r for r in forced if r.original == r.final]
    if missing_forced:
        print("\n  pinned, state already matched the input name:")
        for r in missing_forced:
            print(f"  {r.chain}:{r.resid} {r.final}")
    if result.termini:
        print("\nchain termini:")
        for note in result.termini:
            print(f"  {note}")
    if result.warnings:
        print("\nwarnings:")
        for w in result.warnings:
            print(f"  ! {w}")
    print(f"\nnext:\n  {result.pdb2gmx_cmd}\n")


def find_gmx(name: str = "gmx") -> str | None:
    """gmx is often not on PATH (GMXRC needed) - look in the usual places."""
    import shutil as _sh
    found = _sh.which(name)
    if found:
        return found
    for cand in (
        "/usr/local/gromacs/bin/gmx", "/usr/local/gromacs/bin/gmx_mpi",
        "/usr/bin/gmx", "/opt/gromacs/bin/gmx",
    ):
        if os.path.exists(cand):
            return cand
    return None


def run_pdb2gmx(result: Result, args) -> int:
    # pdb2gmx looks for the force field in the current directory - work there
    workdir = os.path.dirname(os.path.abspath(result.ff_dir))
    ffname = os.path.basename(result.ff_dir)[:-3]
    gmx = find_gmx(args.gmx)
    if gmx is None:
        print(f"cannot find the executable '{args.gmx}' (source GMXRC first)",
              file=sys.stderr)
        return 3
    cmd = [
        gmx, "pdb2gmx", "-f", os.path.abspath(result.output_pdb),
        "-o", os.path.join(os.path.abspath(args.outdir), "conf.gro"),
        "-p", os.path.join(os.path.abspath(args.outdir), "topol.top"),
        "-i", os.path.join(os.path.abspath(args.outdir), "posre.itp"),
        "-ff", ffname, "-water", args.water,
    ]
    print(f"\n>>> {' '.join(cmd)}   (in {workdir})")
    proc = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
    log = os.path.join(args.outdir, "pdb2gmx.log")
    with open(log, "w") as fh:
        fh.write(proc.stdout + "\n" + proc.stderr)
    if proc.returncode == 0:
        print("pdb2gmx finished successfully, the topology is built.")
        for line in proc.stderr.splitlines():
            if "Total charge" in line or "charge" in line.lower():
                print("  " + line.strip())
    else:
        print(f"pdb2gmx failed (exit code {proc.returncode}), details in {log}")
        tail = (proc.stderr or proc.stdout).strip().splitlines()[-25:]
        for line in tail:
            print("  " + line)
    return proc.returncode


def main(argv: List[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s: %(message)s",
    )
    try:
        spec = build_spec(args)
        result = prepare(
            input_pdb=args.pdb,
            spec=spec,
            ff_dir=args.ff,
            outdir=args.outdir,
            output_name=args.name,
            keep_water=not args.drop_water,
            keep_het=not args.drop_het,
            debump=not args.no_debump,
            opt=not args.no_opt,
            water_model=args.water,
            verbose=args.verbose,
        )
    except ForceFieldError as err:
        print("\nStopping: the force field lacks what is needed, and I will "
              "not edit it.\n", file=sys.stderr)
        for problem in err.problems:
            print(f"  * {problem}", file=sys.stderr)
        print("\nNothing was written.", file=sys.stderr)
        return 3
    except (SpecError, ValueError, FileNotFoundError) as err:
        print(f"error: {err}", file=sys.stderr)
        return 2
    write_report(result, spec, args.outdir)
    summarize(result, spec)
    if args.run_pdb2gmx:
        return run_pdb2gmx(result, args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
