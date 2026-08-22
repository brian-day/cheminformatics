"""Molecule standardization: cleanup, salt stripping, charge neutralization, and tautomer canonicalization.

Structures drawn from different sources rarely agree, bit-for-bit, on how the "same" molecule should be
represented: one source might store a hydrochloride salt, another the free base; one might draw a
carboxylic acid as -COOH, another as the carboxylate -COO-; a compound with a tautomerizable group (e.g.
keto/enol) can be drawn either way. None of these differences change what the molecule *is*, but they do
change what RDKit computes for it (LogP, TPSA, formal charge, ...) unless everything is normalized to a
single representation first. `standardize_molecule` runs a molecule through that normalization pipeline.
"""

from rdkit.Chem import Mol
from rdkit.Chem.MolStandardize import rdMolStandardize

_TAUTOMER_ENUMERATOR = rdMolStandardize.TautomerEnumerator()
_UNCHARGER = rdMolStandardize.Uncharger()


def standardize_molecule(mol: Mol) -> Mol:
    """Normalize a molecule to a canonical representation.

    Applies, in order:

    1. **Cleanup** — recomputes valences and aromaticity and removes redundant explicit
       hydrogens, fixing up minor inconsistencies in how the input was drawn.
    2. **Parent fragment selection** — keeps only the largest covalently-bonded fragment,
       discarding counterions and other salt/solvate fragments (e.g. "drug.HCl" becomes
       just "drug"). Without this, a salt form and its free base would get different
       molecular weights and descriptor values despite being the "same" active compound.
    3. **Charge neutralization** — neutralizes atoms that can be neutralized without
       changing the molecule's connectivity (e.g. a carboxylate -COO- becomes -COOH).
       Formal charge state otherwise depends on how the input happened to be drawn or
       what pH it was modeled at, which isn't a meaningful structural difference here.
    4. **Tautomer canonicalization** — picks a single canonical tautomer, so that two
       inputs differing only in tautomeric form (e.g. keto vs. enol) standardize to the
       same structure rather than being treated as different molecules.

    Run this before computing descriptors or fingerprints on molecules from external
    sources: salt form, charge state, and tautomeric form all affect those calculations
    even though none of them reflect a genuine structural difference between compounds.
    """
    mol = rdMolStandardize.Cleanup(mol)
    mol = rdMolStandardize.FragmentParent(mol)
    mol = _UNCHARGER.uncharge(mol)
    mol = _TAUTOMER_ENUMERATOR.Canonicalize(mol)
    return mol
