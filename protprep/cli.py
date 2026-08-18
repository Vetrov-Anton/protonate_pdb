"""Command line interface."""

from __future__ import annotations

import argparse
import logging
import sys
from typing import List

from . import report as report_mod
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
    report_mod.write_reports(result, outdir, spec)


def summarize(result: Result, spec: Spec) -> None:
    print()
    print(result.summary())
    print()


def run_pdb2gmx(result: Result, args) -> int:
    try:
        run = result.run_pdb2gmx(outdir=args.outdir, water=args.water,
                                 gmx=args.gmx)
    except FileNotFoundError as err:
        print(str(err), file=sys.stderr)
        return 3
    print(f"\n>>> {' '.join(run.command)}   (in {run.workdir})")
    if run.ok:
        print("pdb2gmx finished successfully, the topology is built.")
        for line in run.charge_lines():
            print("  " + line)
    else:
        print(f"pdb2gmx failed (exit code {run.returncode}), "
              f"details in {run.log_path}")
        for line in run.tail():
            print("  " + line)
    return run.returncode


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
