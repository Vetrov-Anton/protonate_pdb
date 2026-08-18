# Examples

* `fragment.pdb` — residues 30–70 of chain A from 6CFO plus the crystal waters
  (~15 KB). Enough to check in a second that everything works.
* `protprep_demo.ipynb` — a notebook walking through the Python API: pinning
  states, reading the report, standard pKa, hydrogens only on pinned residues,
  caps, building the topology, and what an unsupported state looks like.
* `example_spec.txt` — a spec file for the full 6CFO homotetramer: pinned side
  chains in both alpha subunits, caps on chains A and C.

Quick check (the force field is taken from the installed GROMACS by name):

```bash
protonate -f examples/fragment.pdb -o /tmp/prot_demo --ph 7.0 \
    --ff amber99sb-ildn --fix A:31:p --fix A:63:HIE
```

Expected: Asp31 becomes `ASH`, His63 becomes `HIE`, all other groups follow
PROPKA, and the chain termini stay charged. Add `--run-pdb2gmx` and the
topology is built right away.

The full 6CFO (downloaded from RCSB, not part of the repository):

```bash
wget https://files.rcsb.org/download/6CFO.pdb
protonate -f 6CFO.pdb -o prepared --spec examples/example_spec.txt \
    --ff amber99sb-ildn --drop-het --run-pdb2gmx
```
