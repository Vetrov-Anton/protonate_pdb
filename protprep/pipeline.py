"""Main pipeline: a hydrogen-free PDB -> a PDB ready for pdb2gmx."""

from __future__ import annotations

import logging
import os
import shutil
import tempfile
from dataclasses import dataclass, field, replace
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from . import chemistry, ffdata, termini as termini_mod
from .pdbio import (STANDARD_AA, WATER, Atom, group_residues, parse_pdb,
                    write_pdb)
from .spec import (Spec, SpecError, TITRATABLE, normalize_residue_state)

LOG = logging.getLogger("protprep")

# states pdb2pqr can produce on its own
PDB2PQR_STATES = {"ASP", "ASH", "GLU", "GLH", "LYS", "LYN", "CYS", "CYM",
                  "CYX", "HID", "HIE", "HIP", "TYR", "ARG"}

FORCED_PKA_HIGH = 1000.0    # pH < pKa -> protonated form
FORCED_PKA_LOW = -1000.0    # pH >= pKa -> deprotonated form

# Standard (model) pKa values of free amino acids - the same numbers PROPKA
# uses as its reference. Needed for --standard-pka, where no local shifts are
# computed at all.
MODEL_PKA: Dict[str, float] = {
    "ASP": 3.80, "GLU": 4.50, "HIS": 6.50, "CYS": 9.00,
    "TYR": 10.00, "LYS": 10.50, "ARG": 12.50,
}
MODEL_PKA_NTERM = 8.00
MODEL_PKA_CTERM = 3.20


@dataclass
class ResidueReport:
    chain: str
    resid: int
    icode: str
    original: str
    final: str
    pka: Optional[float] = None
    model_pka: Optional[float] = None
    buried: Optional[float] = None
    forced: bool = False
    notes: List[str] = field(default_factory=list)


class ForceFieldError(RuntimeError):
    """The force field lacks something. We do not touch it - we just refuse."""

    def __init__(self, problems: List[str]):
        self.problems = problems
        super().__init__("; ".join(problems))


@dataclass
class Result:
    """Everything a run produced: the structure, the per-residue report and
    the ready-made pdb2gmx command."""

    output_pdb: str
    residues: List[ResidueReport]
    termini: List[str]
    warnings: List[str]
    ff_dir: str
    pdb2gmx_cmd: str
    spec: Optional[Spec] = None

    # ------------------------------------------------------------ queries
    @property
    def states(self) -> Dict[Tuple[str, int], str]:
        """{(chain, resid): final residue name} for every titratable residue."""
        return {(r.chain, r.resid): r.final for r in self.residues}

    def changed(self) -> List[ResidueReport]:
        """Residues whose state differs from the input file."""
        return [r for r in self.residues if r.original != r.final]

    def fixed(self) -> List[ResidueReport]:
        """Residues whose state was pinned by the user."""
        return [r for r in self.residues if r.forced]

    def records(self) -> List[dict]:
        """The report as plain dicts (handy for json/pandas)."""
        return [
            {
                "chain": r.chain, "resid": r.resid, "icode": r.icode,
                "input": r.original, "final": r.final, "pKa": r.pka,
                "model_pKa": r.model_pka, "buried": r.buried,
                "fixed": r.forced, "notes": r.notes,
            }
            for r in self.residues
        ]

    def to_dataframe(self, only_changed: bool = False):
        """pandas.DataFrame with the report (pandas must be installed)."""
        import pandas as pd

        rows = self.records()
        if only_changed:
            rows = [r for r in rows if r["input"] != r["final"]]
        return pd.DataFrame(rows)

    # ------------------------------------------------------------ actions
    def save_reports(self, outdir: Optional[str] = None) -> List[str]:
        """Write protonation_report.tsv / .json next to the structure."""
        from . import report as report_mod

        outdir = outdir or os.path.dirname(os.path.abspath(self.output_pdb))
        return report_mod.write_reports(self, outdir, self.spec)

    def run_pdb2gmx(self, outdir: Optional[str] = None, water: str = "tip3p",
                    gmx: str = "gmx", extra_args: Optional[List[str]] = None):
        """Build the topology with gmx pdb2gmx. Returns a Pdb2gmxResult."""
        from . import gmx as gmx_mod

        outdir = outdir or os.path.dirname(os.path.abspath(self.output_pdb))
        return gmx_mod.run_pdb2gmx(
            self.output_pdb, self.ff_dir, outdir, water=water, gmx=gmx,
            extra_args=extra_args,
        )

    def summary(self) -> str:
        """The same overview the command line prints, as a string."""
        standard = getattr(self.spec, "pka_source", "propka") == "standard"
        ph = getattr(self.spec, "ph", float("nan"))
        source = ("standard pKa (PROPKA was not run)" if standard
                  else "local pKa from PROPKA")
        changed, forced = self.changed(), self.fixed()
        out = [
            f"pH = {ph:.2f}; {source}",
            f"structure: {self.output_pdb}",
            f"residues with a shifted state: {len(changed)} "
            f"(pinned by hand: {len(forced)})",
        ]
        if changed:
            label = "pKa(table)" if standard else "pKa(PROPKA)"
            out.append("")
            out.append(f"  chain residue    was -> now   {label}  source")
            for r in changed:
                pka = f"{r.pka:8.2f}" if r.pka is not None else "       -"
                src = ("PINNED" if r.forced else ("table" if standard else "propka"))
                out.append(
                    f"  {r.chain:>5} {r.resid:>7}   {r.original:>4} -> "
                    f"{r.final:<5} {pka}   {src}"
                )
        already = [r for r in forced if r.original == r.final]
        if already:
            out.append("")
            out.append("  pinned, state already matched the input name:")
            out += [f"  {r.chain}:{r.resid} {r.final}" for r in already]
        if self.termini:
            out.append("")
            out.append("chain termini:")
            out += [f"  {n}" for n in self.termini]
        if self.warnings:
            out.append("")
            out.append("warnings:")
            out += [f"  ! {w}" for w in self.warnings]
        out.append("")
        out.append(f"next:\n  {self.pdb2gmx_cmd}")
        return "\n".join(out)

    def __str__(self) -> str:  # pragma: no cover - convenience only
        return self.summary()


# ---------------------------------------------------------------- pdb2pqr
def _model_pka_rows(biomolecule) -> Tuple[List[dict], str]:
    """Stand-in for the PROPKA calculation: tabulated pKa for every group."""
    from collections import OrderedDict

    rows = []
    for residue in biomolecule.residues:
        name = str(getattr(residue, "name", "")).upper()
        pka = MODEL_PKA.get(name)
        if pka is None:
            continue
        rows.append(
            OrderedDict(
                res_num=residue.res_seq,
                ins_code=(getattr(residue, "ins_code", "") or "").strip(),
                res_name=name,
                chain_id=residue.chain_id,
                group_label=name,
                group_type=None,
                pKa=pka,
                model_pKa=pka,
                buried=None,
                coupled_group=None,
            )
        )
    return rows, "PROPKA was not run: standard (model) pKa values were used"


def run_pdb2pqr(
    in_pdb: str,
    out_pdb: str,
    ph: float,
    forced: Dict[str, float],
    debump: bool = True,
    opt: bool = True,
    log_level: str = "ERROR",
    standard_pka: bool = False,
) -> List[dict]:
    """Run pdb2pqr with the pKa values of the pinned residues overridden.

    :param standard_pka: do not run PROPKA, use tabulated pKa instead - then
        every unpinned group, termini included, ends up in the state the
        standard values imply at this pH.
    """
    from pdb2pqr import biomolecule as pqr_biomolecule
    from pdb2pqr import main as pqr_main

    original = pqr_biomolecule.Biomolecule.apply_pka_values
    original_propka = pqr_main.run_propka

    def patched(self, force_field, ph_, pkadic):
        for key, value in forced.items():
            pkadic[key] = value
        return original(self, force_field, ph_, pkadic)

    def patched_propka(_args, biomolecule):
        return _model_pka_rows(biomolecule)

    parser = pqr_main.build_main_parser()
    argv = [
        "--ff=AMBER", "--ffout=AMBER", "--keep-chain",
        "--titration-state-method", "propka", "--with-ph", str(ph),
        "--pdb-output", out_pdb, "--log-level", log_level,
    ]
    if not debump:
        argv.append("--nodebump")
    if not opt:
        argv.append("--noopt")
    argv += [in_pdb, out_pdb + ".pqr"]
    args = parser.parse_args(argv)

    pqr_biomolecule.Biomolecule.apply_pka_values = patched
    if standard_pka:
        pqr_main.run_propka = patched_propka
    # main_driver() does not apply --log-level itself, so silence the chatty
    # pdb2pqr/propka loggers here (they are back to normal with -v/verbose)
    level = getattr(logging, log_level.upper(), logging.ERROR)
    noisy = [
        logging.getLogger(name)
        for name in list(logging.root.manager.loggerDict)
        if name.lower().startswith(("pdb2pqr", "propka"))
    ]
    saved = [lg.level for lg in noisy]
    for lg in noisy:
        lg.setLevel(level)
    try:
        _missed, pka_rows, _bio = pqr_main.main_driver(args)
    finally:
        pqr_biomolecule.Biomolecule.apply_pka_values = original
        pqr_main.run_propka = original_propka
        for lg, lvl in zip(noisy, saved):
            lg.setLevel(lvl)
    return pka_rows or []


# -------------------------------------------------------------- splitting
def _split_input(atoms: Sequence[Atom], keep_water: bool, keep_het: bool):
    """(residues for pdb2pqr, everything else), filtering altlocs."""
    titratable_like, other = [], []
    for atom in atoms:
        if atom.altloc not in (" ", "", "A"):
            continue
        name = atom.resname.upper()
        if name in STANDARD_AA and atom.record == "ATOM":
            if atom.element.upper() == "H" or (
                not atom.element and atom.name.lstrip("0123456789").startswith("H")
            ):
                continue  # drop input hydrogens, pdb2pqr will place them
            titratable_like.append(atom)
        elif name in WATER:
            if keep_water:
                if atom.element.upper() != "H":
                    titratable_like.append(replace(atom, record="ATOM"))
        elif keep_het:
            other.append(atom)
    return titratable_like, other


# ------------------------------------------------------------------ main
def prepare(
    input_pdb: str,
    spec: Spec,
    ff_dir: str,
    outdir: str,
    output_name: str = "protonated.pdb",
    keep_water: bool = True,
    keep_het: bool = True,
    debump: bool = True,
    opt: bool = True,
    water_model: str = "tip3p",
    verbose: bool = False,
) -> Result:
    warnings: List[str] = []
    ff = ffdata.ForceField(ff_dir)

    atoms, _header = parse_pdb(input_pdb)
    protein, other = _split_input(atoms, keep_water, keep_het)
    if not protein:
        raise ValueError("No standard residue found in the structure")

    orig_names = {
        key: res[0].resname.upper() for key, res in group_residues(protein)
    }

    # --- 1. what the user asked for --------------------------------------
    targets: Dict[Tuple[str, int, str], str] = {}
    forced_pka: Dict[str, float] = {}
    for item in spec.residues:
        wildcard = item.chain in ("", "*")
        matched = 0
        for key in _resolve_keys(item, orig_names, os.path.basename(input_pdb)):
            parent = orig_names[key]
            parent_family = {"HID": "HIS", "HIE": "HIS",
                             "HIP": "HIS"}.get(parent, parent)
            try:
                if parent_family not in TITRATABLE:
                    raise SpecError(
                        f"{key[0]}:{key[1]} is a {parent}, which is not titratable"
                    )
                target = normalize_residue_state(parent_family, item.state)
            except SpecError:
                if not wildcard:
                    raise
                # chain '*': other chains may hold anything at that number -
                # such residues are simply skipped
                warnings.append(
                    f"{key[0]}:{key[1]} ({parent}) skipped: state "
                    f"'{item.state}' does not apply to it"
                )
                continue
            matched += 1
            targets[key] = target
            pqr_key = f"{parent} {key[1]} {key[0]}".strip()
            branch_prot = target in ("ASH", "GLH", "HIP", "LYS", "CYS", "CYX",
                                     "TYR", "ARG")
            forced_pka[pqr_key] = (FORCED_PKA_HIGH if branch_prot
                                   else FORCED_PKA_LOW)
        if wildcard and matched == 0:
            raise SpecError(
                f"Residue {item.resid}{item.icode}: state '{item.state}' "
                "does not apply in any chain"
            )

    # --- 2. pdb2pqr + propka ---------------------------------------------
    tmpdir = tempfile.mkdtemp(prefix="protprep_")
    try:
        stage_in = os.path.join(tmpdir, "input.pdb")
        stage_out = os.path.join(tmpdir, "pqr.pdb")
        write_pdb(stage_in, protein)
        standard = spec.pka_source == "standard"
        LOG.info(
            "Running pdb2pqr at pH %.2f (%s) ...", spec.ph,
            "standard pKa" if standard else "local pKa from PROPKA",
        )
        pka_rows = run_pdb2pqr(
            stage_in, stage_out, spec.ph, forced_pka, debump=debump, opt=opt,
            log_level="INFO" if verbose else "ERROR", standard_pka=standard,
        )
        prot_atoms, _ = parse_pdb(stage_out)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    pka_map: Dict[Tuple[str, int, str], dict] = {}
    for row in pka_rows:
        label = str(row.get("group_label", ""))
        if not label.startswith(str(row.get("res_name", ""))):
            continue
        pka_map[(str(row["chain_id"]).strip(), int(row["res_num"]),
                 str(row.get("ins_code", "")).strip())] = row

    # --- 3. finishing the states by hand ---------------------------------
    residues = group_residues(prot_atoms)
    reports: List[ResidueReport] = []
    cloud = chemistry.coords(prot_atoms)
    fixed_residues: List[Tuple[Tuple[str, int, str], List[Atom]]] = []

    for key, res in residues:
        current = res[0].resname.upper()
        notes: List[str] = []
        target = targets.get(key)
        if target and target != current:
            res, notes = chemistry.set_state(list(res), target, cloud)
            for note in notes:
                if note.startswith("FAILED"):
                    warnings.append(f"{key[0]}:{key[1]} {target}: {note}")
        elif target:
            notes.append("state came straight out of pdb2pqr")
        fixed_residues.append((key, list(res)))
        if current in WATER or orig_names.get(key, "") in WATER:
            continue          # water does not belong in the protonation report
        row = pka_map.get(key)
        reports.append(
            ResidueReport(
                chain=key[0], resid=key[1], icode=key[2],
                original=orig_names.get(key, current),
                final=res[0].resname.upper(),
                pka=row["pKa"] if row else None,
                model_pka=row["model_pKa"] if row else None,
                buried=row.get("buried") if row else None,
                forced=key in targets,
                notes=notes,
            )
        )

    for key, target in targets.items():
        final = next((r.final for r in reports
                      if (r.chain, r.resid, r.icode) == key), None)
        if final != target:
            warnings.append(
                f"{key[0]}:{key[1]} asked for {target}, got {final}"
            )

    # --- 4. termini -------------------------------------------------------
    if spec.pka_source == "standard":
        # in the amber port both termini exist only in their charged form
        if spec.ph > MODEL_PKA_NTERM:
            warnings.append(
                f"at pH {spec.ph:.2f} a standard N-terminus (pKa "
                f"{MODEL_PKA_NTERM:.1f}) would be neutral, but the force field "
                "has no neutral N-terminus - it stays NH3+ (or use "
                "--nter CHAIN:ACE)"
            )
        if spec.ph < MODEL_PKA_CTERM:
            warnings.append(
                f"at pH {spec.ph:.2f} a standard C-terminus (pKa "
                f"{MODEL_PKA_CTERM:.1f}) would be protonated, but the force "
                "field has no COOH terminus - it stays COO- (or use "
                "--cter CHAIN:NME)"
            )
    chain_order: List[str] = []
    by_chain: Dict[str, List[int]] = {}
    for idx, (key, res) in enumerate(fixed_residues):
        if res[0].resname.upper() in WATER:
            continue
        by_chain.setdefault(key[0], []).append(idx)
        if key[0] not in chain_order:
            chain_order.append(key[0])

    term_notes: List[str] = []
    problems: List[str] = []                       # blocking problems
    caps: Dict[int, Tuple[List[Atom], str]] = {}   # index -> (cap atoms, 'before'|'after')
    positions: Dict[int, str] = {}                 # index -> nter/cter/middle

    for chain in chain_order:
        idxs = by_chain[chain]
        first, last = idxs[0], idxs[-1]
        positions.setdefault(first, "nter")
        positions.setdefault(last, "cter")

        n_state = spec.terminus(chain, "N") or "NH3+"
        c_state = spec.terminus(chain, "C") or "COO-"

        # --- N-terminus
        res = fixed_residues[first][1]
        if n_state == "ACE":
            cap, body, notes = chemistry.cap_nterm_ace(list(res), cloud)
            if cap:
                fixed_residues[first] = (fixed_residues[first][0], body)
                caps[first] = (cap, "before")
                positions[first] = "middle"
            term_notes += [f"chain {chain}: {n}" for n in notes]
        elif n_state == "NH2":
            body = [a for a in res if a.name != "H3"]
            resname = body[0].resname.upper()
            new_name = termini_mod.check_neutral_block(ff, resname, "N")
            if new_name is None:
                problems.append(
                    f"chain {chain}: a neutral N-terminus (NH2) for {resname} "
                    f"is impossible - force field "
                    f"{os.path.basename(ff.path)} has no "
                    f"[ {termini_mod.neutral_nter_name(resname)} ] block. "
                    f"Use a cap (--nter {chain}:ACE) or add such a block to "
                    "the force field yourself."
                )
            else:
                fixed_residues[first] = (
                    fixed_residues[first][0],
                    [replace(a, resname=new_name) for a in body],
                )
                positions[first] = "middle"   # the block is found by name
                term_notes.append(
                    f"chain {chain}: neutral N-terminus (NH2, {new_name})"
                )
        else:
            term_notes.append(f"chain {chain}: charged N-terminus (NH3+)")

        # --- C-terminus
        res = fixed_residues[last][1]
        if c_state in ("NME", "NHE"):
            cap, body, notes = chemistry.cap_cterm(list(res), cloud, kind=c_state)
            if cap:
                fixed_residues[last] = (fixed_residues[last][0], body)
                caps[last] = (cap, "after")
                positions[last] = "middle"
            term_notes += [f"chain {chain}: {n}" for n in notes]
        elif c_state == "COOH":
            body = list(res)
            resname = body[0].resname.upper()
            new_name = termini_mod.check_neutral_block(ff, resname, "C")
            if new_name is None:
                problems.append(
                    f"chain {chain}: a neutral C-terminus (COOH) for {resname} "
                    f"is impossible - force field "
                    f"{os.path.basename(ff.path)} has no "
                    f"[ {termini_mod.neutral_cter_name(resname)} ] block. "
                    f"Use a cap (--cter {chain}:NME or {chain}:NHE) or add "
                    "such a block to the force field yourself."
                )
            else:
                names = {a.name for a in body}
                if "OXT" in names and "HO" not in names:
                    idx = {a.name: a for a in body}
                    pos = chemistry.place_atom(
                        idx["O"], idx["C"], idx["OXT"], 0.97, 113.0, 0.0
                    )
                    body.append(chemistry._template(body, "HO", pos))
                fixed_residues[last] = (
                    fixed_residues[last][0],
                    [replace(a, resname=new_name) for a in body],
                )
                positions[last] = "middle"
                term_notes.append(
                    f"chain {chain}: neutral C-terminus (COOH, {new_name})"
                )
        else:
            term_notes.append(f"chain {chain}: charged C-terminus (COO-)")

    # --- 5. matching atom names to the force field nomenclature ----------
    final_atoms: List[Atom] = []
    stripped_res = stripped_h = 0
    for idx, (key, res) in enumerate(fixed_residues):
        if idx in caps and caps[idx][1] == "before":
            final_atoms.extend(caps[idx][0])
        resname = res[0].resname.upper()
        if resname in WATER:
            # pdb2pqr emits water as WAT/OW/HW/HW (both H share a name);
            # bring it to what the rtp expects: HOH / OW HW1 HW2
            hcount = 0
            for atom in res:
                if atom.element.upper() == "H":
                    hcount += 1
                    final_atoms.append(
                        replace(atom, resname="HOH", name=f"HW{hcount}",
                                record="HETATM")
                    )
                else:
                    final_atoms.append(
                        replace(atom, resname="HOH", name="OW", record="HETATM")
                    )
            continue
        position = positions.get(idx, "middle")
        block = ff.block_for(resname, position)
        where = f"{key[0]}:{key[1]}{key[2]} {resname}"
        if block is None:
            problems.append(
                f"{where}: force field {os.path.basename(ff.path)} has no rtp "
                f"block for this residue ({position})"
            )
            final_atoms.extend(res)
            continue
        if not _terminal_block_ok(block, position, resname):
            end = "N" if position == "nter" else "C"
            flag = "--nter" if position == "nter" else "--cter"
            cap = "ACE" if position == "nter" else "NME"
            problems.append(
                f"{where}: the force field has no {end}-terminal block for "
                f"this state (only [ {block.name} ] was found). Options: cap "
                f"the terminus ({flag} {key[0]}:{cap}) or leave the residue in "
                "its standard state."
            )
        parents = _hydrogen_parents(res)
        rename, extra, missing = ffdata.reconcile_residue(
            [a.name for a in res], parents, block
        )
        # "hydrogens only on pinned residues" mode: the rest lose their H,
        # pdb2gmx will rebuild them from the rtp block named by the residue
        drop_h = spec.hydrogens == "fixed" and key not in targets
        for atom in res:
            if drop_h and atom.element.upper() == "H":
                stripped_h += 1
                continue
            final_atoms.append(replace(atom, name=rename.get(atom.name, atom.name)))
        if drop_h:
            stripped_res += 1
        if extra:
            problems.append(
                f"{where} (block {block.name}): atoms {', '.join(extra)} have "
                "no counterpart in the rtp - pdb2gmx will complain about them"
            )
        if missing:
            warnings.append(
                f"{where} (block {block.name}): atoms {', '.join(missing)} are "
                "absent (pdb2gmx will add the hydrogens itself)"
            )
        if idx in caps and caps[idx][1] == "after":
            final_atoms.extend(caps[idx][0])

    # --- 6. checks that pdb2gmx would otherwise trip over -----------------
    protein_used = {a.resname.upper() for a in final_atoms}
    final_atoms.extend(other)
    het_used = {a.resname.upper() for a in other}
    used = protein_used | het_used

    rt_path, known = termini_mod.system_residuetypes()
    if known:
        # an unknown *protein* name is fatal: GROMACS would treat the residue as
        # non-protein and break the chain in two
        unknown = sorted(n for n in protein_used if n not in known)
        if unknown:
            problems.append(
                "residue names " + ", ".join(unknown) + " are absent from "
                f"residuetypes.dat ({rt_path}) - GROMACS will treat them as "
                "non-protein and break the chain. Either avoid these states or "
                "add the names to residuetypes.dat by hand."
            )
        # ligands and ions are carried over untouched, so an unknown name there
        # is the user's business - just say it out loud
        unknown_het = sorted(n for n in het_used - protein_used if n not in known)
        if unknown_het:
            warnings.append(
                "ligands/ions " + ", ".join(unknown_het) + " are unknown to "
                f"residuetypes.dat ({rt_path}) and carry no hydrogens: they "
                "need their own rtp/hdb entries, or drop them (keep_het=False "
                "/ --drop-het)"
            )

    if spec.hydrogens == "fixed":
        term_notes.append(
            f"hydrogens kept only on the residues from --fix and on the caps; "
            f"{stripped_h} H removed from the other {stripped_res} residues - "
            "pdb2gmx will rebuild them (residue names, and therefore states, "
            "are preserved)"
        )

    for block, delta, ref in ffdata.check_pair_charges(ff, used):
        problems.append(
            f"block [ {block} ] in the force field carries a non-integer "
            f"charge: it must equal the charge of {ref} minus 1, discrepancy "
            f"{delta:+.4f} e. This is a bug in the force field itself - fix it "
            "yourself or do not use this state."
        )

    if problems:
        raise ForceFieldError(problems)

    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, output_name)
    title = [
        "REMARK   1 Prepared by protprep (pdb2pqr/propka + fixed states)",
        f"REMARK   1 pH = {spec.ph:.2f}; force field: {os.path.basename(ff.path)}",
    ]
    write_pdb(out_path, final_atoms, title)

    cmd = _pdb2gmx_command(out_path, ff, water_model)
    return Result(
        output_pdb=out_path, residues=reports, termini=term_notes,
        warnings=warnings, ff_dir=ff.path, pdb2gmx_cmd=cmd, spec=spec,
    )


def _resolve_keys(item, orig_names: Dict[Tuple[str, int, str], str],
                  fname: str) -> List[Tuple[str, int, str]]:
    """Residue keys a spec line refers to.

    Chain '*' (or an empty chain) means "in every chain" - handy for
    homo-oligomers.
    """
    exact = (item.chain, item.resid, item.icode)
    if item.chain not in ("", "*"):
        if exact not in orig_names:
            raise SpecError(
                f"Residue {item.chain}:{item.resid}{item.icode} not found in {fname}"
            )
        return [exact]
    keys = [k for k in orig_names
            if k[1] == item.resid and k[2] == item.icode]
    if not keys:
        raise SpecError(f"Residue {item.resid}{item.icode} not found in {fname}")
    return keys


def _terminal_block_ok(block: ffdata.RtpEntry, position: str, resname: str) -> bool:
    """Does the block we found look like a real terminal one (not an internal)?"""
    names = set(block.atom_names)
    if position == "nter":
        return "H1" in names and ("H2" in names or resname == "PRO")
    if position == "cter":
        return bool(names & {"OC1", "OC2", "OXT", "HO"})
    return True


def _hydrogen_parents(res: Sequence[Atom]) -> Dict[str, str]:
    """H -> nearest heavy atom of the same residue."""
    heavy = [a for a in res if a.element.upper() != "H"]
    if not heavy:
        return {}
    hxyz = np.array([[a.x, a.y, a.z] for a in heavy])
    out: Dict[str, str] = {}
    for atom in res:
        if atom.element.upper() != "H":
            continue
        d = np.linalg.norm(hxyz - np.array([atom.x, atom.y, atom.z]), axis=1)
        out[atom.name] = heavy[int(np.argmin(d))].name
    return out


def _pdb2gmx_command(out_pdb: str, ff: ffdata.ForceField, water: str) -> str:
    ffname = os.path.basename(ff.path)[:-3]
    where = os.path.dirname(os.path.abspath(ff.path))
    return (
        f"cd {where} && gmx pdb2gmx -f {os.path.abspath(out_pdb)} "
        f"-o conf.gro -p topol.top -ff {ffname} -water {water}"
    )
