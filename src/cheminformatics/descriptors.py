"""Molecular descriptor calculation and Lipinski-style drug-likeness filters."""

from dataclasses import dataclass

from rdkit.Chem import Crippen, Descriptors, Lipinski, Mol


@dataclass
class DescriptorSet:
    molecular_weight: float
    logp: float
    tpsa: float
    h_bond_donors: int
    h_bond_acceptors: int
    rotatable_bonds: int
    ring_count: int
    aromatic_ring_count: int
    heavy_atom_count: int

    def to_dict(self) -> dict:
        return self.__dict__


def compute_descriptors(mol: Mol) -> DescriptorSet:
    """Compute a standard set of physicochemical descriptors for a molecule."""
    return DescriptorSet(
        molecular_weight=round(Descriptors.MolWt(mol), 2),
        logp=round(Crippen.MolLogP(mol), 2),
        tpsa=round(Descriptors.TPSA(mol), 2),
        h_bond_donors=Lipinski.NumHDonors(mol),
        h_bond_acceptors=Lipinski.NumHAcceptors(mol),
        rotatable_bonds=Lipinski.NumRotatableBonds(mol),
        ring_count=Lipinski.RingCount(mol),
        aromatic_ring_count=Lipinski.NumAromaticRings(mol),
        heavy_atom_count=mol.GetNumHeavyAtoms(),
    )


def lipinski_violations(descriptors: DescriptorSet) -> list[str]:
    """Return a list of violated Lipinski rule-of-five criteria (empty if fully compliant)."""
    violations = []
    if descriptors.molecular_weight > 500:
        violations.append("molecular_weight > 500")
    if descriptors.logp > 5:
        violations.append("logp > 5")
    if descriptors.h_bond_donors > 5:
        violations.append("h_bond_donors > 5")
    if descriptors.h_bond_acceptors > 10:
        violations.append("h_bond_acceptors > 10")
    return violations
