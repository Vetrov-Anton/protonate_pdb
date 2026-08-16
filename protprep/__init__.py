"""protprep - prepare protein structures for pdb2gmx with fixed protonation
states for selected residues.

The scheme: pdb2pqr/PROPKA computes local pKa values at the requested pH, while
for residues listed by the user those values are overridden with +-1000, which
pins the desired form. States pdb2pqr cannot produce (histidine tautomers,
deprotonated TYR/ARG, ASH/GLH at chain termini) are then finished off
geometrically, ACE/NME/NHE caps are built, and atom names are matched against
the nomenclature of the specific GROMACS force field.
"""

__version__ = "1.1.0"
