"""Командный интерфейс."""

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
примеры:
  # Asp145 протонирован, His264 - тавтомер HIE, остальное - по PROPKA
  protonate -f 6CFO.pdb --ph 7.4 --fix A:145:p --fix A:264:HIE

  # закрыть концы цепи A кэпами, цепь B оставить заряженной
  protonate -f 6CFO.pdb --nter A:ACE --cter A:NME

  # то же самое через файл задания
  protonate -f 6CFO.pdb --spec fixed.txt

формат файла задания (текст или json):
  pH 7.4
  A 145 p               # ASP -> ASH
  A 264 HIE             # тавтомер гистидина
  B 87  d               # LYS -> LYN
  *  63 HIE             # во всех цепях сразу
  nter A ACE
  cter A NME

источник pKa:
  по умолчанию - PROPKA (локальные сдвиги). С --standard-pka PROPKA не
  запускается вовсе: все незафиксированные группы, включая N- и C-концы,
  берутся по стандартным pKa свободных аминокислот при заданном pH.

состояния:
  p = protonated, d = deprotonated, n = neutral, c = charged
  или прямое имя: ASP|ASH, GLU|GLH, LYS|LYN, CYS|CYM|CYX, TYR|TYN,
                  ARG|ARN, HID|HIE|HIP   (у HIS: d = HID, e = HIE, p = HIP)
концы:
  N -> NH3+ (p, c) | NH2 (n, d) | ACE ;  C -> COO- (d, c) | COOH (p, n) | NME | NHE

Силовое поле используется только на чтение. Если в нём нет нужного блока,
скрипт сообщает об этом и не пишет ничего.
"""


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="protonate",
        description="Протонирование белка при заданном pH (PROPKA) с "
                    "принудительно заданными состояниями отдельных остатков.",
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-f", "--pdb", required=True, help="входной PDB (без водородов)")
    p.add_argument("-o", "--outdir", default="prepared", help="каталог результата")
    p.add_argument("-n", "--name", default="protonated.pdb", help="имя выходного PDB")
    p.add_argument("--ph", type=float, default=None, help="pH среды (по умолчанию 7.0)")
    p.add_argument("--standard-pka", "--no-propka", dest="standard_pka",
                   action="store_true",
                   help="не запускать PROPKA: все незафиксированные группы, "
                        "включая N- и C-концы, берутся в состоянии по "
                        "стандартным (табличным) pKa при заданном pH")
    p.add_argument("--ff", default="amber-99sb-ildn.ff",
                   help="каталог силового поля GROMACS (*.ff)")
    p.add_argument("--spec", help="файл задания (txt или json)")
    p.add_argument("--fix", action="append", default=[], metavar="ЦЕПЬ:НОМЕР:СОСТОЯНИЕ",
                   help="зафиксировать состояние остатка (можно повторять)")
    p.add_argument("--nter", action="append", default=[], metavar="ЦЕПЬ:СОСТОЯНИЕ",
                   help="состояние N-конца цепи: NH3+|NH2|ACE")
    p.add_argument("--cter", action="append", default=[], metavar="ЦЕПЬ:СОСТОЯНИЕ",
                   help="состояние C-конца цепи: COO-|COOH|NME|NHE")
    p.add_argument("--drop-water", action="store_true", help="выбросить воду")
    p.add_argument("--drop-het", action="store_true",
                   help="выбросить лиганды/ионы (HETATM)")
    p.add_argument("--no-debump", action="store_true", help="pdb2pqr без debump")
    p.add_argument("--no-opt", action="store_true",
                   help="pdb2pqr без оптимизации водородных связей")
    p.add_argument("--water", default="tip3p", help="модель воды для pdb2gmx")
    p.add_argument("--run-pdb2gmx", action="store_true",
                   help="сразу проверить результат запуском gmx pdb2gmx")
    p.add_argument("--gmx", default="gmx", help="исполняемый файл gromacs")
    p.add_argument("-v", "--verbose", action="store_true")
    return p


def build_spec(args) -> Spec:
    spec = load_spec(args.spec) if args.spec else Spec()
    if args.ph is not None:
        spec.ph = args.ph
    if args.standard_pka:
        spec.pka_source = "standard"
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
    source = ("стандартные pKa (PROPKA не запускалась)"
              if spec.pka_source == "standard" else "локальные pKa из PROPKA")
    print(f"\npH = {spec.ph:.2f}; {source}")
    print(f"структура: {result.output_pdb}")
    print(f"остатков со сдвинутым состоянием: {len(changed)} "
          f"(из них зафиксировано вручную: {len(forced)})")
    if changed:
        label = "pKa(табл.)" if spec.pka_source == "standard" else "pKa(PROPKA)"
        print(f"\n  цепь остаток   было -> стало   {label}  источник")
        for r in changed:
            pka = f"{r.pka:8.2f}" if r.pka is not None else "       -"
            src = ("ЗАДАНО" if r.forced
                   else ("таблица" if spec.pka_source == "standard" else "propka"))
            print(f"  {r.chain:>4} {r.resid:>7}   {r.original:>4} -> {r.final:<5} "
                  f"{pka}   {src}")
    missing_forced = [r for r in forced if r.original == r.final]
    if missing_forced:
        print("\n  зафиксированы, состояние совпало с исходным именем:")
        for r in missing_forced:
            print(f"  {r.chain}:{r.resid} {r.final}")
    if result.termini:
        print("\nконцы цепей:")
        for note in result.termini:
            print(f"  {note}")
    if result.warnings:
        print("\nпредупреждения:")
        for w in result.warnings:
            print(f"  ! {w}")
    print(f"\nдальше:\n  {result.pdb2gmx_cmd}\n")


def find_gmx(name: str = "gmx") -> str | None:
    """gmx часто не в PATH (нужен GMXRC) - поищем в типичных местах."""
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
    workdir = os.path.dirname(os.path.abspath(result.output_pdb))
    # pdb2gmx ищет силовое поле в текущем каталоге - работаем оттуда
    workdir = os.path.dirname(os.path.abspath(result.ff_dir))
    ffname = os.path.basename(result.ff_dir)[:-3]
    gmx = find_gmx(args.gmx)
    if gmx is None:
        print(f"не нашёл исполняемый файл '{args.gmx}' (нужен source GMXRC)",
              file=sys.stderr)
        return 3
    cmd = [
        gmx, "pdb2gmx", "-f", os.path.abspath(result.output_pdb),
        "-o", os.path.join(os.path.abspath(args.outdir), "conf.gro"),
        "-p", os.path.join(os.path.abspath(args.outdir), "topol.top"),
        "-i", os.path.join(os.path.abspath(args.outdir), "posre.itp"),
        "-ff", ffname, "-water", args.water,
    ]
    print(f"\n>>> {' '.join(cmd)}   (в каталоге {workdir})")
    proc = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
    log = os.path.join(args.outdir, "pdb2gmx.log")
    with open(log, "w") as fh:
        fh.write(proc.stdout + "\n" + proc.stderr)
    if proc.returncode == 0:
        print("pdb2gmx отработал успешно, топология собрана.")
        for line in proc.stderr.splitlines():
            if "Total charge" in line or "charge" in line.lower():
                print("  " + line.strip())
    else:
        print(f"pdb2gmx завершился с ошибкой (код {proc.returncode}), "
              f"подробности: {log}")
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
        print("\nОстанавливаюсь: силовое поле не содержит нужного, а править "
              "его я не буду.\n", file=sys.stderr)
        for problem in err.problems:
            print(f"  * {problem}", file=sys.stderr)
        print("\nНичего не записано.", file=sys.stderr)
        return 3
    except (SpecError, ValueError, FileNotFoundError) as err:
        print(f"ошибка: {err}", file=sys.stderr)
        return 2
    write_report(result, spec, args.outdir)
    summarize(result, spec)
    if args.run_pdb2gmx:
        return run_pdb2gmx(result, args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
