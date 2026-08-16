# protprep — protonation with pinned residue states

[![tests](https://github.com/Vetrov-Anton/protonate_pdb/actions/workflows/ci.yml/badge.svg)](https://github.com/Vetrov-Anton/protonate_pdb/actions/workflows/ci.yml)

Prepares a hydrogen-free PDB for `gmx pdb2gmx`: local pKa values come from
**PROPKA/PDB2PQR** at the requested pH, but the states of the residues you list
are set **by hand** — something neither PROPKA nor PDB2PQR lets you do.
On top of that it does what PDB2PQR refuses to do for AMBER: histidine
tautomers, deprotonated TYR/ARG, and ACE/NME/NHE caps.

> **The force field is opened read-only.** The tool never writes or copies
> anything into a `*.ff` directory and never touches `residuetypes.dat`. If a
> required block is missing, it stops with an explanation and writes nothing.

## Installation

One command, no cloning required:

```bash
pipx install git+https://github.com/Vetrov-Anton/protonate_pdb.git    # recommended
# or, if you have no pipx:
pip install --user git+https://github.com/Vetrov-Anton/protonate_pdb.git
```

This puts a `protonate` command on your PATH. The dependencies (`pdb2pqr`,
`propka`, `numpy`) are installed automatically; the only separate requirement
is GROMACS, which is where the force fields come from.

Upgrade with `pipx upgrade protprep`, remove with `pipx uninstall protprep`.

<details>
<summary>If you want to hack on the code</summary>

```bash
git clone https://github.com/Vetrov-Anton/protonate_pdb.git
cd protonate_pdb
./install.sh          # venv next to the sources (editable), ./protonate command
./install.sh --pipx   # isolated, via pipx
./install.sh --user   # pip install --user
python -m pytest -q   # tests (need an installed GROMACS)
```

The `./protonate` script also works with no installation at all: it picks up a
local `.venv`, falling back to the system python.
</details>

## Quick start

A check on the small example from the repository (~1 s; if you installed via
pipx, use your own PDB or grab `examples/fragment.pdb`):

```bash
protonate -f examples/fragment.pdb -o prepared --ph 7.0 \
    --ff amber99sb-ildn --fix A:31:p --fix A:63:HIE
```

What it looks like on a real structure:

```bash
protonate -f 6CFO.pdb -o prepared --ph 7.4 --ff amber99sb-ildn \
    --fix A:167:p --fix A:63:HIE --fix B:99:d \
    --nter A:ACE --cter A:NME \
    --run-pdb2gmx
```

The same thing can live in a spec file (see `examples/example_spec.txt`):

```bash
protonate -f 6CFO.pdb -o prepared --spec examples/example_spec.txt --run-pdb2gmx
```

The force field is looked up in three places in order: the path you gave, the
current directory, and the installed GROMACS (`$GMXDATA/top`,
`/usr/local/gromacs/...`). So `--ff amber99sb-ildn` works from anywhere, while
`--ff ./my-custom.ff` picks up your own. If nothing matches, the tool prints the
list of force fields it can see.

## Pinning states

`--fix CHAIN:NUMBER:STATE`, repeatable. Instead of a chain you may write `*` —
the rule then applies to that residue number in every chain (handy for
homo-oligomers; residues of a different type are skipped with a warning).

Shorthand: **`p`** = protonated, **`d`** = deprotonated, **`n`** = neutral,
**`c`** = charged. The full words work too.

| residue | `p` (protonated) | `d` (deprotonated) | `n` | `c` |
|---------|------------------|--------------------|-----|-----|
| ASP     | `ASH` (COOH)     | `ASP` (COO−)       | `ASH` | `ASP` |
| GLU     | `GLH`            | `GLU`              | `GLH` | `GLU` |
| LYS     | `LYS` (NH3+)     | `LYN` (NH2)        | `LYN` | `LYS` |
| CYS     | `CYS` (SH)       | `CYM` (S−)         | `CYS` | `CYM` |
| TYR     | `TYR` (OH)       | `TYN` (O−)         | `TYR` | `TYN` |
| ARG     | `ARG` (+)        | `ARN` (neutral)    | `ARN` | `ARG` |
| HIS     | `HIP` (+)        | `HID`              | `HIE` | `HIP` |

For histidine the letters denote tautomers: **`d`** = `HID` (proton on ND1),
**`e`** = `HIE` (on NE2), **`p`** = `HIP` (both, charged), `n` = `HIE`.
Explicit names (`ASH`, `LYN`, `CYM`, `HID`, …) are accepted as well.
`CYX` means a cysteine in a disulfide (no HG).

## Only the selected side chains, everything else from standard pKa

`--standard-pka` (a.k.a. `--no-propka`) switches PROPKA off completely: you pin
the side chains you care about with `--fix`, and **every other group, N- and
C-termini included**, takes the state the standard pKa values of free amino
acids imply at the requested pH — with no environment-induced shifts at all.

```bash
protonate -f 6CFO.pdb -o prepared --ph 7.0 --standard-pka \
    --fix A:167:p --fix A:63:HIE
```

The table (the same values PROPKA uses as its reference) and the outcome at
pH 7:

| group      | standard pKa | state at pH 7 |
|------------|-------------:|---------------|
| ASP        | 3.80  | `ASP` (COO−) |
| GLU        | 4.50  | `GLU` (COO−) |
| HIS        | 6.50  | neutral; pdb2pqr picks the `HID`/`HIE` tautomer from the H-bond network |
| CYS        | 9.00  | `CYS` (SH) |
| TYR        | 10.00 | `TYR` (OH) |
| LYS        | 10.50 | `LYS` (NH3+) |
| ARG        | 12.50 | `ARG` (+) |
| N-terminus | 8.00  | `NH3+` |
| C-terminus | 3.20  | `COO−` |

It is also faster (~1 s instead of ~11 s on 6CFO) — PROPKA is not run at all.
The flag works at any pH, not only 7; if at your pH a standard terminus ought to
change state but the force field has no such variant, you get a warning (see the
termini section). In a spec file this is the line `pka standard`.

## Hydrogens only on the selected residues

`--only-fixed-h` (a.k.a. `--strip-other-h`) keeps hydrogens **only on the
residues from `--fix`** (and on the ACE/NME/NHE caps, if you asked for them).
Every other residue goes into the file bare, and `gmx pdb2gmx` builds its
hydrogens:

```bash
protonate -f 6CFO.pdb -o prepared --ph 7.0 --standard-pka --only-fixed-h \
    --fix A:167:p --fix A:63:HIE --ff amber99sb-ildn
```

Note that **residue names are preserved**, so the states do not go anywhere: if
PROPKA decided that Glu263 is `GLH`, it leaves as `GLH`, just without
hydrogens, and pdb2gmx builds them from the `GLH` block. What is stripped are
the hydrogen coordinates, not the protonation decision. If you want everything
except `--fix` in its standard state, add `--standard-pka`.

Crystal waters keep their hydrogens (PDB2PQR placed and optimised them).

Verified on 6CFO: with and without the flag pdb2gmx produces **identical
topologies** (only the comment lines with paths and dates differ) and the same
total charge; the hydrogen coordinates of the pinned residues match
character for character — pdb2gmx leaves them alone and only adds what is
missing.

## Chain termini

```
--nter A:NH3+ | NH2 | ACE          shorthand: p, c = NH3+ ; n, d = NH2
--cter A:COO- | COOH | NME | NHE   shorthand: d, c = COO- ; p, n = COOH
```

* `NH3+` / `COO-` — charged termini, the default behaviour;
* `ACE` / `NME` / `NHE` — caps built geometrically (trans peptide bond,
  standard lengths and angles, methyl rotation chosen to minimise contacts);
  their hydrogens are placed too;
* `NH2` / `COOH` — neutral uncapped termini. **The GROMACS amber ports do not
  have them**: `aminoacids.n.tdb` / `aminoacids.c.tdb` are empty and the
  terminus type comes from separate rtp blocks (`NALA`, `CALA`), among which no
  neutral variants exist. Such a request therefore fails with a suggestion to
  use a cap. If you added the blocks to your own force field, name them `ZXXX`
  (neutral N-terminus) and `JXXX` (neutral C-terminus, COOH) — e.g. `ZALA`,
  `JALA` — and they will be picked up automatically.

For non-standard states **at the very end of a chain** (a protonated ASP/GLU and
the like) the force field has neither `NASH` nor `CASH` — the tool catches this
and stops, suggesting a cap.

## What you get

In the `-o` directory:

| file | what it is |
|------|------------|
| `protonated.pdb` | structure with hydrogens, named per your force field |
| `protonation_report.tsv` | table: residue, PROPKA pKa, final state, pinned or not |
| `protonation_report.json` | the same plus warnings and the pdb2gmx command |
| `pdb2gmx.log`, `conf.gro`, `topol.top` | if you ran with `--run-pdb2gmx` |

Then it is plain `gmx pdb2gmx`; the ready command is printed at the end and
stored in the JSON report. The `-ignh` flag is **not needed** and not wanted:
atom names already match the rtp, so pdb2gmx takes our hydrogens as they are
(visible in the line "Now there are N residues with M atoms" — the atom count
does not change).

## How it works

1. Standard residues and water are taken from the input PDB; hydrogens and
   altlocs other than A are dropped; ligands/ions (`HETATM`) are set aside and
   appended back untouched (`--drop-het`, `--drop-water` if you do not want
   them).
2. PDB2PQR runs (`--ff=AMBER --ffout=AMBER`) with PROPKA at your pH. Right
   before the pKa values are applied the dictionary is patched: pinned residues
   get ±1000, so the desired form is built by PDB2PQR itself — with its own
   hydrogen bond optimisation and debumping. With `--standard-pka` it is not the
   dictionary but the calculation that is replaced: PROPKA gives way to a table
   of model pKa values, and nothing else in the pipeline changes.
3. Whatever PDB2PQR refuses to do for AMBER (HIS tautomer, `TYN`, `ARN`,
   `ASH`/`GLH` at the termini) is finished off geometrically: the proton is
   placed from internal coordinates (bond/angle/torsion), or removed.
4. Caps are built and atom names are matched against the `*.rtp` files of your
   force field. The mapping is derived from the rtp itself (heavy atoms by
   element, hydrogens by their parent from the bond list), so local edits of the
   force field do not break the tool.
5. Final checks — and if something does not add up, it stops without writing:
   * no rtp block for a residue or for its position in the chain;
   * atoms left over after matching that the rtp does not know;
   * a residue name unknown to `residuetypes.dat` (GROMACS would treat it as
     non-protein);
   * a non-integer block charge. `TYN` in this amber port, for example, carries
     −0.399 instead of −1 — a bug in the force field, and yours to fix.

## Limitations worth knowing

* **Disulfides** are still detected by pdb2gmx (by distance). To pin one, set
  `CYX` explicitly.
* **Ligands.** `A5X`, `MG` and other HETATM records are carried over as they
  are, without hydrogens: they need their own rtp/hdb entries (your force field
  already ships `TPP.rtp`, `API.rtp`, `unl_GMX.rtp` — the tool does not check
  whether they suffice).
* **PROPKA is not symmetric.** On the 6CFO homotetramer chains B and D come out
  with different states (B206 → HIE, D206 → HID, for instance). If you need
  symmetry, pin those residues explicitly with `*`.
* **`TYN` and `ARN` cannot be used with this force field**: they are absent from
  `residuetypes.dat`, and `TYN` also has a broken charge. The tool says so and
  stops.
* Water is protonated by PDB2PQR and renamed to `HOH`/`OW`/`HW1`/`HW2`.

## All the options

```
-f/--pdb            input PDB
-o/--outdir         output directory (default: prepared)
-n/--name           output PDB name
--ph                pH of the medium
--standard-pka      skip PROPKA: unpinned groups (and termini) follow the
                    standard tabulated pKa (= --no-propka)
--only-fixed-h      hydrogens only on the residues from --fix and on the caps,
                    pdb2gmx builds the rest (= --strip-other-h)
--ff                force field: a path or a bare name; looked up nearby and
                    inside GROMACS
--spec              spec file (txt or json)
--fix               CHAIN:NUMBER:STATE, repeatable
--nter/--cter       CHAIN:STATE, repeatable
--drop-water        discard water
--drop-het          discard ligands and ions
--no-debump/--no-opt  turn off the corresponding PDB2PQR stages
--water             water model for pdb2gmx (tip3p)
--run-pdb2gmx       check the result by building the topology right away
--gmx               path to gmx (also searched in /usr/local/gromacs/bin)
```

Exit codes: `0` — success, `2` — bad request or input, `3` — the force field
lacks something (nothing was written).
