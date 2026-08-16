"""Основной конвейер: PDB без водородов -> PDB, готовый для pdb2gmx."""

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

# внутренние состояния, которые pdb2pqr сам умеет выдавать
PDB2PQR_STATES = {"ASP", "ASH", "GLU", "GLH", "LYS", "LYN", "CYS", "CYM",
                  "CYX", "HID", "HIE", "HIP", "TYR", "ARG"}

FORCED_PKA_HIGH = 1000.0    # pH < pKa -> протонированная форма
FORCED_PKA_LOW = -1000.0    # pH >= pKa -> депротонированная форма

# Стандартные (модельные) pKa свободных аминокислот - те же значения, которые
# PROPKA использует как точку отсчёта. Нужны для режима --standard-pka, когда
# расчёт локальных сдвигов не проводится вовсе.
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
    """Силовому полю чего-то не хватает. Мы его не трогаем - просто отказ."""

    def __init__(self, problems: List[str]):
        self.problems = problems
        super().__init__("; ".join(problems))


@dataclass
class Result:
    output_pdb: str
    residues: List[ResidueReport]
    termini: List[str]
    warnings: List[str]
    ff_dir: str
    pdb2gmx_cmd: str


# ---------------------------------------------------------------- pdb2pqr
def _model_pka_rows(biomolecule) -> Tuple[List[dict], str]:
    """Замена расчёту PROPKA: табличные (модельные) pKa для всех групп."""
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
    return rows, "PROPKA не запускался: взяты стандартные (модельные) pKa"


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
    """Запускает pdb2pqr, подменив pKa у зафиксированных остатков.

    :param standard_pka: не запускать PROPKA, а взять табличные pKa - тогда
        все незафиксированные группы (включая термини) получают состояние,
        которое им положено при данном pH по стандартным значениям.
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
    try:
        _missed, pka_rows, _bio = pqr_main.main_driver(args)
    finally:
        pqr_biomolecule.Biomolecule.apply_pka_values = original
        pqr_main.run_propka = original_propka
    return pka_rows or []


# ------------------------------------------------------------- сортировка
def _split_input(atoms: Sequence[Atom], keep_water: bool, keep_het: bool):
    """(остатки для pdb2pqr, прочее) с фильтром altloc."""
    titratable_like, other = [], []
    for atom in atoms:
        if atom.altloc not in (" ", "", "A"):
            continue
        name = atom.resname.upper()
        if name in STANDARD_AA and atom.record == "ATOM":
            if atom.element.upper() == "H" or (
                not atom.element and atom.name.lstrip("0123456789").startswith("H")
            ):
                continue  # входные водороды выбрасываем, их поставит pdb2pqr
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
        raise ValueError("В структуре не найдено ни одного стандартного остатка")

    orig_names = {
        key: res[0].resname.upper() for key, res in group_residues(protein)
    }

    # --- 1. чего хочет пользователь -------------------------------------
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
                        f"{key[0]}:{key[1]} - это {parent}, он не титруется"
                    )
                target = normalize_residue_state(parent_family, item.state)
            except SpecError:
                if not wildcard:
                    raise
                # цепь '*': в других цепях под этим номером может стоять
                # что угодно - такие остатки просто пропускаем
                warnings.append(
                    f"{key[0]}:{key[1]} ({parent}) пропущен: состояние "
                    f"'{item.state}' к нему неприменимо"
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
                f"Остаток {item.resid}{item.icode}: ни в одной цепи состояние "
                f"'{item.state}' неприменимо"
            )

    # --- 2. pdb2pqr + propka --------------------------------------------
    tmpdir = tempfile.mkdtemp(prefix="protprep_")
    try:
        stage_in = os.path.join(tmpdir, "input.pdb")
        stage_out = os.path.join(tmpdir, "pqr.pdb")
        write_pdb(stage_in, protein)
        standard = spec.pka_source == "standard"
        LOG.info(
            "Запускаю pdb2pqr при pH %.2f (%s) ...", spec.ph,
            "стандартные pKa" if standard else "локальные pKa из PROPKA",
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

    # --- 3. доводка состояний руками ------------------------------------
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
                if note.startswith("НЕ добавлен"):
                    warnings.append(f"{key[0]}:{key[1]} {target}: {note}")
        elif target:
            notes.append("состояние получено сразу из pdb2pqr")
        fixed_residues.append((key, list(res)))
        if current in WATER or orig_names.get(key, "") in WATER:
            continue          # вода в отчёт о протонировании не идёт
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
                f"{key[0]}:{key[1]} запрошено {target}, получено {final}"
            )

    # --- 4. термини ------------------------------------------------------
    if spec.pka_source == "standard":
        # в amber-порте оба конца существуют только в заряженном виде
        if spec.ph > MODEL_PKA_NTERM:
            warnings.append(
                f"при pH {spec.ph:.2f} стандартный N-конец (pKa "
                f"{MODEL_PKA_NTERM:.1f}) был бы нейтрален, но в силовом поле "
                "нейтрального N-конца нет - остаётся NH3+ (или задайте "
                "--nter ЦЕПЬ:ACE)"
            )
        if spec.ph < MODEL_PKA_CTERM:
            warnings.append(
                f"при pH {spec.ph:.2f} стандартный C-конец (pKa "
                f"{MODEL_PKA_CTERM:.1f}) был бы протонирован, но в силовом поле "
                "COOH-конца нет - остаётся COO- (или задайте --cter ЦЕПЬ:NME)"
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
    problems: List[str] = []                       # блокирующие проблемы
    caps: Dict[int, Tuple[List[Atom], str]] = {}   # индекс -> (атомы кэпа, 'before'|'after')
    positions: Dict[int, str] = {}                 # индекс -> nter/cter/middle

    for chain in chain_order:
        idxs = by_chain[chain]
        first, last = idxs[0], idxs[-1]
        positions.setdefault(first, "nter")
        positions.setdefault(last, "cter")

        n_state = spec.terminus(chain, "N") or "NH3+"
        c_state = spec.terminus(chain, "C") or "COO-"

        # --- N-конец
        res = fixed_residues[first][1]
        if n_state == "ACE":
            cap, body, notes = chemistry.cap_nterm_ace(list(res), cloud)
            if cap:
                fixed_residues[first] = (fixed_residues[first][0], body)
                caps[first] = (cap, "before")
                positions[first] = "middle"
            term_notes += [f"цепь {chain}: {n}" for n in notes]
        elif n_state == "NH2":
            body = [a for a in res if a.name != "H3"]
            resname = body[0].resname.upper()
            new_name = termini_mod.check_neutral_block(ff, resname, "N")
            if new_name is None:
                problems.append(
                    f"цепь {chain}: нейтральный N-конец (NH2) для {resname} "
                    f"невозможен - в силовом поле "
                    f"{os.path.basename(ff.path)} нет блока "
                    f"[ {termini_mod.neutral_nter_name(resname)} ]. "
                    "Возьмите кэп (--nter "
                    f"{chain}:ACE) либо добавьте такой блок в ff сами."
                )
            else:
                fixed_residues[first] = (
                    fixed_residues[first][0],
                    [replace(a, resname=new_name) for a in body],
                )
                positions[first] = "middle"   # блок ищем по имени напрямую
                term_notes.append(
                    f"цепь {chain}: N-конец нейтральный (NH2, {new_name})"
                )
        else:
            term_notes.append(f"цепь {chain}: N-конец заряжен (NH3+)")

        # --- C-конец
        res = fixed_residues[last][1]
        if c_state in ("NME", "NHE"):
            cap, body, notes = chemistry.cap_cterm(list(res), cloud, kind=c_state)
            if cap:
                fixed_residues[last] = (fixed_residues[last][0], body)
                caps[last] = (cap, "after")
                positions[last] = "middle"
            term_notes += [f"цепь {chain}: {n}" for n in notes]
        elif c_state == "COOH":
            body = list(res)
            resname = body[0].resname.upper()
            new_name = termini_mod.check_neutral_block(ff, resname, "C")
            if new_name is None:
                problems.append(
                    f"цепь {chain}: нейтральный C-конец (COOH) для {resname} "
                    f"невозможен - в силовом поле "
                    f"{os.path.basename(ff.path)} нет блока "
                    f"[ {termini_mod.neutral_cter_name(resname)} ]. "
                    f"Возьмите кэп (--cter {chain}:NME или {chain}:NHE) либо "
                    "добавьте такой блок в ff сами."
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
                    f"цепь {chain}: C-конец нейтральный (COOH, {new_name})"
                )
        else:
            term_notes.append(f"цепь {chain}: C-конец заряжен (COO-)")

    # --- 5. приведение имён атомов к номенклатуре силового поля ----------
    final_atoms: List[Atom] = []
    stripped_res = stripped_h = 0
    for idx, (key, res) in enumerate(fixed_residues):
        if idx in caps and caps[idx][1] == "before":
            final_atoms.extend(caps[idx][0])
        resname = res[0].resname.upper()
        if resname in WATER:
            # pdb2pqr выдаёт воду как WAT/OW/HW/HW (имена H совпадают);
            # приводим к виду, который ждёт rtp: HOH / OW HW1 HW2
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
                f"{where}: в силовом поле {os.path.basename(ff.path)} нет "
                f"блока rtp для этого остатка ({position})"
            )
            final_atoms.extend(res)
            continue
        if not _terminal_block_ok(block, position, resname):
            end = "N" if position == "nter" else "C"
            flag = "--nter" if position == "nter" else "--cter"
            cap = "ACE" if position == "nter" else "NME"
            problems.append(
                f"{where}: в силовом поле нет {end}-концевого блока для этого "
                f"состояния (нашёлся только [ {block.name} ]). Варианты: "
                f"закрыть конец кэпом ({flag} {key[0]}:{cap}) либо оставить "
                f"остаток в стандартном состоянии."
            )
        parents = _hydrogen_parents(res)
        rename, extra, missing = ffdata.reconcile_residue(
            [a.name for a in res], parents, block
        )
        # режим "водороды только у зафиксированных": у остальных остатков H
        # снимаются, их достроит pdb2gmx по имени остатка из rtp
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
                f"{where} (блок {block.name}): атомам {', '.join(extra)} нет "
                f"соответствия в rtp - pdb2gmx на них ругнётся"
            )
        if missing:
            warnings.append(
                f"{where} (блок {block.name}): нет атомов "
                f"{', '.join(missing)} (pdb2gmx достроит водороды сам)"
            )
        if idx in caps and caps[idx][1] == "after":
            final_atoms.extend(caps[idx][0])

    # --- 6. проверки, после которых pdb2gmx точно споткнётся -------------
    final_atoms.extend(other)
    used = {a.resname.upper() for a in final_atoms}

    rt_path, known = termini_mod.system_residuetypes()
    unknown = sorted(n for n in used if known and n not in known)
    if unknown:
        problems.append(
            "имена остатков " + ", ".join(unknown) + " отсутствуют в "
            f"residuetypes.dat ({rt_path}) - GROMACS сочтёт их не белком и "
            "разорвёт цепь. Либо не используйте эти состояния, либо добавьте "
            "имена в residuetypes.dat вручную."
        )

    if spec.hydrogens == "fixed":
        term_notes.append(
            f"водороды оставлены только у остатков из --fix и у кэпов; "
            f"с {stripped_res} остальных остатков снято {stripped_h} H - "
            "их достроит pdb2gmx (имена остатков, а значит и состояния, "
            "сохранены)"
        )

    for block, delta, ref in ffdata.check_pair_charges(ff, used):
        problems.append(
            f"блок [ {block} ] в силовом поле имеет нецелый заряд: он должен "
            f"быть равен заряду {ref} минус 1, расхождение {delta:+.4f} e. "
            "Это ошибка самого силового поля - исправьте её сами или не "
            "используйте это состояние."
        )

    if problems:
        raise ForceFieldError(problems)

    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, output_name)
    title = [
        "REMARK   1 Подготовлено protprep (pdb2pqr/propka + фиксированные состояния)",
        f"REMARK   1 pH = {spec.ph:.2f}; силовое поле: {os.path.basename(ff.path)}",
    ]
    write_pdb(out_path, final_atoms, title)

    cmd = _pdb2gmx_command(out_path, ff, water_model)
    return Result(
        output_pdb=out_path, residues=reports, termini=term_notes,
        warnings=warnings, ff_dir=ff.path, pdb2gmx_cmd=cmd,
    )


def _resolve_keys(item, orig_names: Dict[Tuple[str, int, str], str],
                  fname: str) -> List[Tuple[str, int, str]]:
    """Ключи остатков, к которым относится строка задания.

    Цепь '*' (или пустая) означает "во всех цепях" - удобно для гомоолигомеров.
    """
    exact = (item.chain, item.resid, item.icode)
    if item.chain not in ("", "*"):
        if exact not in orig_names:
            raise SpecError(
                f"Остаток {item.chain}:{item.resid}{item.icode} не найден в {fname}"
            )
        return [exact]
    keys = [k for k in orig_names
            if k[1] == item.resid and k[2] == item.icode]
    if not keys:
        raise SpecError(f"Остаток {item.resid}{item.icode} не найден в {fname}")
    return keys


def _terminal_block_ok(block: ffdata.RtpEntry, position: str, resname: str) -> bool:
    """Похож ли найденный блок на настоящий концевой (а не на внутренний)?"""
    names = set(block.atom_names)
    if position == "nter":
        return "H1" in names and ("H2" in names or resname == "PRO")
    if position == "cter":
        return bool(names & {"OC1", "OC2", "OXT", "HO"})
    return True


def _hydrogen_parents(res: Sequence[Atom]) -> Dict[str, str]:
    """H -> ближайший тяжёлый атом того же остатка."""
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
