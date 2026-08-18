"""The Python API: Protonator, protonate() and the Result helpers."""

from __future__ import annotations

import os

import pytest

from protprep import ForceFieldError, Protonator, protonate
from protprep.ffdata import gromacs_top_dirs, resolve_ff_path

HERE = os.path.dirname(os.path.abspath(__file__))
FRAGMENT = os.path.join(HERE, os.pardir, "examples", "fragment.pdb")


def _ff():
    if not gromacs_top_dirs():
        pytest.skip("GROMACS is not installed - no force field to borrow")
    try:
        return resolve_ff_path("amber99sb-ildn")
    except FileNotFoundError:
        pytest.skip("no amber99sb-ildn.ff")


def test_chaining_and_states(tmp_path):
    prot = (
        Protonator(FRAGMENT, ff=_ff(), ph=7.0)
        .fix("A", 31, "p")          # three arguments
        .fix("A:63:HID")            # single token
        .fix_many([("A", 54, "d"), "A:62:CYM"])
    )
    result = prot.run(str(tmp_path))

    assert result.states[("A", 31)] == "ASH"
    assert result.states[("A", 63)] == "HID"
    assert result.states[("A", 54)] == "LYN"
    assert result.states[("A", 62)] == "CYM"
    assert {(r.chain, r.resid) for r in result.fixed()} == {
        ("A", 31), ("A", 63), ("A", 54), ("A", 62)
    }
    assert os.path.exists(result.output_pdb)
    # reports are written by default
    assert os.path.exists(tmp_path / "protonation_report.tsv")
    assert os.path.exists(tmp_path / "protonation_report.json")


def test_summary_and_records(tmp_path):
    result = Protonator(FRAGMENT, ff=_ff(), pka="standard").fix("A:31:p").run(
        str(tmp_path)
    )
    text = result.summary()
    assert "ASP -> ASH" in text and "PINNED" in text
    assert "standard pKa" in text

    records = result.records()
    assert records and set(records[0]) >= {"chain", "resid", "input", "final"}
    assert all(r["fixed"] is False or r["resid"] == 31 for r in records)


def test_caps_and_only_fixed_hydrogens(tmp_path):
    result = (
        Protonator(FRAGMENT, ff=_ff())
        .standard_pka()
        .only_fixed_hydrogens()
        .fix("A", 31, "p")
        .cap("A", n="ACE", c="NME")
        .run(str(tmp_path))
    )
    per_residue = {}
    for line in open(result.output_pdb):
        if not line.startswith("ATOM"):
            continue
        key = (line[21], int(line[22:26]))
        counts = per_residue.setdefault(key, [0, line[17:20].strip()])
        counts[0] += int(line[76:78].strip() == "H")

    names = {name for _, name in per_residue.values()}
    assert {"ACE", "NME"} <= names
    with_h = {k for k, (n_h, _) in per_residue.items() if n_h}
    # the pinned residue plus the two caps keep their hydrogens
    assert ("A", 31) in with_h and len(with_h) == 3


def test_one_shot_helper(tmp_path):
    result = protonate(
        FRAGMENT, str(tmp_path), fix=["A:31:p", ("A", 63, "HIE")],
        nter=[("A", "ACE")], ff=_ff(), ph=7.0, pka="standard",
    )
    assert result.states[("A", 31)] == "ASH"
    assert result.states[("A", 63)] == "HIE"
    assert any("ACE" in note for note in result.termini)


def test_spec_file_roundtrip(tmp_path):
    spec_file = tmp_path / "spec.txt"
    spec_file.write_text(
        "pH 7.0\npka standard\nhydrogens fixed\nA 31 p\nnter A ACE\n"
    )
    prot = Protonator.from_spec_file(FRAGMENT, str(spec_file), ff=_ff())
    assert prot.spec.pka_source == "standard"
    assert prot.spec.hydrogens == "fixed"
    result = prot.run(str(tmp_path / "out"))
    assert result.states[("A", 31)] == "ASH"


def test_force_field_error_is_raised(tmp_path):
    prot = Protonator(FRAGMENT, ff=_ff()).fix("A", 35, "d")   # TYR -> TYN
    with pytest.raises(ForceFieldError) as excinfo:
        prot.run(str(tmp_path / "out"))
    assert excinfo.value.problems          # the reasons are listed
    assert not os.path.exists(tmp_path / "out")
