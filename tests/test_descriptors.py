from rdkit import Chem

from cheminformatics.descriptors import compute_descriptors, lipinski_violations

ASPIRIN_SMILES = "CC(=O)OC1=CC=CC=C1C(=O)O"


def test_compute_descriptors_aspirin():
    mol = Chem.MolFromSmiles(ASPIRIN_SMILES)
    descriptors = compute_descriptors(mol)

    assert 179 < descriptors.molecular_weight < 181
    assert descriptors.h_bond_donors == 1
    assert descriptors.ring_count == 1


def test_lipinski_violations_aspirin_is_compliant():
    mol = Chem.MolFromSmiles(ASPIRIN_SMILES)
    descriptors = compute_descriptors(mol)

    assert lipinski_violations(descriptors) == []
