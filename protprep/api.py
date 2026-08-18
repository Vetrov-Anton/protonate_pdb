"""Python API: the same pipeline as the CLI, driven from code.

    from protprep import Protonator

    result = (
        Protonator("6CFO.pdb", ff="amber99sb-ildn", ph=7.4)
        .fix("A", 167, "p")           # ASP -> ASH
        .fix("*", 63, "HIE")          # the same tautomer in every chain
        .cap("A", n="ACE", c="NME")
        .run("prepared")
    )
    print(result.summary())
    topology = result.run_pdb2gmx()
"""

from __future__ import annotations

import os
from typing import Iterable, Optional, Tuple, Union

from .pipeline import ForceFieldError, Result, prepare
from .spec import (ResidueSpec, Spec, SpecError, TerminusSpec, load_spec,
                   parse_fix_token, parse_hydrogens, parse_pka_source,
                   parse_terminus_token)

__all__ = ["Protonator", "protonate", "ForceFieldError", "SpecError", "Result"]

FixItem = Union[str, Tuple]


class Protonator:
    """Fluent wrapper around the protonation pipeline.

    Every configuration method returns ``self``, so calls can be chained; the
    work happens only in :meth:`run`, which returns a :class:`~protprep.pipeline.Result`.

    :param pdb: input structure without hydrogens
    :param ff: force field - a path or a bare name looked up inside GROMACS
    :param ph: pH of the medium
    :param pka: ``"propka"`` (local shifts) or ``"standard"`` (tabulated values)
    :param hydrogens: ``"all"`` or ``"fixed"`` (H only on the pinned residues)
    :param keep_water: carry crystal waters over into the result
    :param keep_het: carry ligands and ions (HETATM) over into the result
    :param debump: let pdb2pqr debump the structure
    :param opt: let pdb2pqr optimise the hydrogen bond network
    :param water_model: water model used in the suggested pdb2gmx command
    """

    def __init__(
        self,
        pdb: str,
        ff: str = "amber-99sb-ildn.ff",
        ph: float = 7.0,
        pka: str = "propka",
        hydrogens: str = "all",
        spec: Optional[Spec] = None,
        keep_water: bool = True,
        keep_het: bool = True,
        debump: bool = True,
        opt: bool = True,
        water_model: str = "tip3p",
        verbose: bool = False,
    ):
        self.pdb = pdb
        self.ff = ff
        self.keep_water = keep_water
        self.keep_het = keep_het
        self.debump = debump
        self.opt = opt
        self.water_model = water_model
        self.verbose = verbose
        self.spec = spec or Spec()
        self.spec.ph = ph
        self.spec.pka_source = parse_pka_source(pka)
        self.spec.hydrogens = parse_hydrogens(hydrogens)

    # ------------------------------------------------------------ loading
    @classmethod
    def from_spec_file(cls, pdb: str, spec_path: str, **kwargs) -> "Protonator":
        """Build from a spec file (the same format the CLI's --spec takes)."""
        spec = load_spec(spec_path)
        kwargs.setdefault("ph", spec.ph)
        kwargs.setdefault("pka", spec.pka_source)
        kwargs.setdefault("hydrogens", spec.hydrogens)
        return cls(pdb, spec=spec, **kwargs)

    # ------------------------------------------------------ configuration
    def fix(self, chain: str, resid: Optional[int] = None,
            state: Optional[str] = None, icode: str = "") -> "Protonator":
        """Pin the protonation state of one residue.

        Accepts either three arguments or a single token::

            .fix("A", 145, "p")
            .fix("A:145:p")
            .fix("*", 63, "HIE")      # every chain at once
        """
        if resid is None and state is None:
            chain, resid, icode, state = parse_fix_token(chain)
        elif state is None:
            raise SpecError("fix() needs a state, e.g. fix('A', 145, 'p')")
        self.spec.residues.append(
            ResidueSpec(str(chain), int(resid), icode.upper(), str(state), str(state))
        )
        return self

    def fix_many(self, items: Iterable[FixItem]) -> "Protonator":
        """Pin several residues: strings ``"A:145:p"`` or tuples ``("A", 145, "p")``."""
        for item in items:
            if isinstance(item, str):
                self.fix(item)
            else:
                self.fix(*item)
        return self

    def nter(self, chain: str, state: str) -> "Protonator":
        """N-terminus of a chain: NH3+ | NH2 | ACE (``*`` = every chain)."""
        self.spec.termini.append(
            TerminusSpec(str(chain), "N", parse_terminus_token("N", state))
        )
        return self

    def cter(self, chain: str, state: str) -> "Protonator":
        """C-terminus of a chain: COO- | COOH | NME | NHE (``*`` = every chain)."""
        self.spec.termini.append(
            TerminusSpec(str(chain), "C", parse_terminus_token("C", state))
        )
        return self

    def cap(self, chain: str = "*", n: Optional[str] = "ACE",
            c: Optional[str] = "NME") -> "Protonator":
        """Cap both ends of a chain in one call."""
        if n:
            self.nter(chain, n)
        if c:
            self.cter(chain, c)
        return self

    def standard_pka(self, on: bool = True) -> "Protonator":
        """Skip PROPKA and use the standard (tabulated) pKa values."""
        self.spec.pka_source = "standard" if on else "propka"
        return self

    def only_fixed_hydrogens(self, on: bool = True) -> "Protonator":
        """Keep hydrogens only on the pinned residues and on the caps."""
        self.spec.hydrogens = "fixed" if on else "all"
        return self

    def at_ph(self, ph: float) -> "Protonator":
        self.spec.ph = float(ph)
        return self

    # -------------------------------------------------------------- work
    def run(self, outdir: str = "prepared", name: str = "protonated.pdb",
            save_reports: bool = True) -> Result:
        """Run the pipeline and write the structure into ``outdir``.

        :raises ForceFieldError: the force field lacks a required block; in
            that case nothing is written at all.
        """
        result = prepare(
            input_pdb=self.pdb,
            spec=self.spec,
            ff_dir=self.ff,
            outdir=outdir,
            output_name=name,
            keep_water=self.keep_water,
            keep_het=self.keep_het,
            debump=self.debump,
            opt=self.opt,
            water_model=self.water_model,
            verbose=self.verbose,
        )
        if save_reports:
            result.save_reports(outdir)
        return result

    # ------------------------------------------------------------- misc
    @property
    def ph(self) -> float:
        return self.spec.ph

    def describe(self) -> str:
        """What this object will do, in one readable block."""
        lines = [
            f"structure : {self.pdb}",
            f"force field: {self.ff}",
            f"pH        : {self.spec.ph:.2f} ({self.spec.pka_source} pKa)",
            f"hydrogens : {self.spec.hydrogens}",
        ]
        if self.spec.residues:
            lines.append("pinned    :")
            lines += [f"  {r.chain}:{r.resid}{r.icode} -> {r.state}"
                      for r in self.spec.residues]
        if self.spec.termini:
            lines.append("termini   :")
            lines += [f"  chain {t.chain} {t.end}-term -> {t.state}"
                      for t in self.spec.termini]
        return "\n".join(lines)

    def __repr__(self) -> str:  # pragma: no cover - convenience only
        return (f"Protonator({os.path.basename(self.pdb)!r}, ph={self.spec.ph}, "
                f"pka={self.spec.pka_source!r}, hydrogens={self.spec.hydrogens!r}, "
                f"fixed={len(self.spec.residues)})")


def protonate(
    pdb: str,
    outdir: str = "prepared",
    fix: Optional[Iterable[FixItem]] = None,
    nter: Optional[Iterable[Tuple[str, str]]] = None,
    cter: Optional[Iterable[Tuple[str, str]]] = None,
    ff: str = "amber-99sb-ildn.ff",
    ph: float = 7.0,
    pka: str = "propka",
    hydrogens: str = "all",
    **kwargs,
) -> Result:
    """One-shot helper for when a whole class feels like too much.

        result = protonate("6CFO.pdb", fix=["A:167:p", ("A", 63, "HIE")],
                           nter=[("A", "ACE")], ph=7.4)
    """
    prot = Protonator(pdb, ff=ff, ph=ph, pka=pka, hydrogens=hydrogens, **kwargs)
    if fix:
        prot.fix_many(fix)
    for chain, state in (nter or []):
        prot.nter(chain, state)
    for chain, state in (cter or []):
        prot.cter(chain, state)
    return prot.run(outdir)
