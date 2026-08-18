"""protprep - prepare protein structures for pdb2gmx with fixed protonation
states for selected residues.

The scheme: pdb2pqr/PROPKA computes local pKa values at the requested pH, while
for residues listed by the user those values are overridden with +-1000, which
pins the desired form. States pdb2pqr cannot produce (histidine tautomers,
deprotonated TYR/ARG, ASH/GLH at chain termini) are then finished off
geometrically, ACE/NME/NHE caps are built, and atom names are matched against
the nomenclature of the specific GROMACS force field.

Command line:

    protonate -f 6CFO.pdb -o prepared --ph 7.4 --fix A:167:p

From Python:

    from protprep import Protonator

    result = Protonator("6CFO.pdb", ph=7.4).fix("A", 167, "p").run("prepared")
    print(result.summary())
"""

from .api import Protonator, protonate
from .gmx import Pdb2gmxResult, find_gmx
from .pipeline import ForceFieldError, ResidueReport, Result, prepare
from .spec import ResidueSpec, Spec, SpecError, TerminusSpec, load_spec

__version__ = "1.2.0"

__all__ = [
    "Protonator",
    "protonate",
    "prepare",
    "Result",
    "ResidueReport",
    "ForceFieldError",
    "Spec",
    "SpecError",
    "ResidueSpec",
    "TerminusSpec",
    "load_spec",
    "Pdb2gmxResult",
    "find_gmx",
    "__version__",
]
