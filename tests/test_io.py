import pytest
from rdkit import Chem

from cheminformatics.io import load_molecule


def test_load_molecule_standardizes_charged_smiles():
    mol = load_molecule("CC(=O)[O-]")  # acetate anion

    assert Chem.GetFormalCharge(mol) == 0
    assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(Chem.MolFromSmiles("CC(=O)O"))


def test_load_molecule_invalid_smiles_raises():
    with pytest.raises(ValueError):
        load_molecule("not a smiles")
