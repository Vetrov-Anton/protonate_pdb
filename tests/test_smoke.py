"""Дымовой тест: прогон фрагмента через весь конвейер.

Запуск: pytest -q   (нужен установленный GROMACS - берём из него силовое поле)
"""

from __future__ import annotations

import os

import pytest

from protprep.ffdata import gromacs_top_dirs, resolve_ff_path
from protprep.pipeline import prepare
from protprep.spec import ResidueSpec, Spec, TerminusSpec

HERE = os.path.dirname(os.path.abspath(__file__))
FRAGMENT = os.path.join(HERE, os.pardir, "examples", "fragment.pdb")


def _ff():
    if not gromacs_top_dirs():
        pytest.skip("GROMACS не установлен - неоткуда взять силовое поле")
    try:
        return resolve_ff_path("amber99sb-ildn")
    except FileNotFoundError:
        pytest.skip("нет amber99sb-ildn.ff")


def _states(result):
    return {(r.chain, r.resid): r.final for r in result.residues}


def test_propka_run_and_fixed_states(tmp_path):
    spec = Spec(ph=7.0)
    spec.residues = [
        ResidueSpec("A", 31, "", "p"),        # ASP -> ASH
        ResidueSpec("A", 63, "", "HID"),      # тавтомер гистидина
        ResidueSpec("A", 54, "", "d"),        # LYS -> LYN
        ResidueSpec("A", 62, "", "CYM"),
    ]
    result = prepare(FRAGMENT, spec, _ff(), str(tmp_path))
    states = _states(result)
    assert states[("A", 31)] == "ASH"
    assert states[("A", 63)] == "HID"
    assert states[("A", 54)] == "LYN"
    assert states[("A", 62)] == "CYM"
    assert os.path.exists(result.output_pdb)
    with open(result.output_pdb) as fh:
        atoms = [l for l in fh if l.startswith(("ATOM", "HETATM"))]
    assert any(l[12:16].strip() == "HD2" and l[17:20] == "ASH" for l in atoms)


def test_standard_pka_leaves_everything_default(tmp_path):
    spec = Spec(ph=7.0, pka_source="standard")
    spec.residues = [ResidueSpec("A", 31, "", "p")]
    result = prepare(FRAGMENT, spec, _ff(), str(tmp_path))
    states = _states(result)
    assert states[("A", 31)] == "ASH"
    # при стандартных pKa и pH 7 сдвигаются только гистидины (в нейтральные)
    for (chain, resid), final in states.items():
        assert final not in ("GLH", "LYN", "CYM", "TYN", "ARN") or resid == 31


def test_caps_are_built(tmp_path):
    spec = Spec(ph=7.0, pka_source="standard")
    spec.termini = [TerminusSpec("A", "N", "ACE"), TerminusSpec("A", "C", "NME")]
    result = prepare(FRAGMENT, spec, _ff(), str(tmp_path))
    with open(result.output_pdb) as fh:
        names = {l[17:20] for l in fh if l.startswith("ATOM")}
    assert "ACE" in names and "NME" in names


def test_only_fixed_hydrogens(tmp_path):
    """Водороды остаются только у остатков из --fix; состояния сохраняются."""
    spec = Spec(ph=7.0, pka_source="standard", hydrogens="fixed")
    spec.residues = [ResidueSpec("A", 31, "", "p"), ResidueSpec("A", 63, "", "HID")]
    result = prepare(FRAGMENT, spec, _ff(), str(tmp_path))

    per_residue = {}
    for line in open(result.output_pdb):
        if not line.startswith("ATOM"):
            continue
        key = (line[21], int(line[22:26]))
        is_h = line[76:78].strip() == "H"
        counts = per_residue.setdefault(key, [0, line[17:20].strip()])
        counts[0] += int(is_h)

    assert per_residue[("A", 31)][0] > 0 and per_residue[("A", 31)][1] == "ASH"
    assert per_residue[("A", 63)][0] > 0 and per_residue[("A", 63)][1] == "HID"
    bare = [k for k, (n_h, _) in per_residue.items() if n_h == 0]
    assert len(bare) == len(per_residue) - 2
    # имена остатков (то есть состояния) остаются на месте
    assert all(name for _, (_, name) in per_residue.items())


def _ff_fingerprint(path):
    out = {}
    for root, _dirs, files in os.walk(path):
        for name in files:
            full = os.path.join(root, name)
            out[os.path.relpath(full, path)] = os.path.getsize(full)
    return out


def test_force_field_is_never_touched(tmp_path):
    """TYR- в amber-портах не собирается; чем бы ни кончилось - ff не тронут."""
    from protprep.pipeline import ForceFieldError

    ff = _ff()
    before = _ff_fingerprint(ff)
    spec = Spec(ph=7.0, pka_source="standard")
    spec.residues = [ResidueSpec("A", 35, "", "d")]   # TYR -> TYN
    outdir = tmp_path / "out"
    try:
        prepare(FRAGMENT, spec, ff, str(outdir))
    except ForceFieldError:
        # штатный отказ: каталог результата даже не создаётся
        assert not os.path.exists(outdir)
    assert _ff_fingerprint(ff) == before
