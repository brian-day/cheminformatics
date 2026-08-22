import pytest
from rdkit import Chem

from cheminformatics.conformers import generate_3d_coordinates


def test_generate_3d_coordinates_produces_a_single_3d_conformer():
    mol = Chem.MolFromSmiles("CCO")
    conf_mol = generate_3d_coordinates(mol)

    assert conf_mol.GetNumConformers() == 1
    assert conf_mol.GetConformer().Is3D()


def test_generate_3d_coordinates_adds_explicit_hydrogens():
    mol = Chem.MolFromSmiles("CCO")  # 3 heavy atoms, 6 hydrogens once explicit
    conf_mol = generate_3d_coordinates(mol)

    assert conf_mol.GetNumAtoms() == 9


def test_generate_3d_coordinates_is_deterministic_given_a_seed():
    mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin

    conf_a = generate_3d_coordinates(mol, seed=1)
    conf_b = generate_3d_coordinates(mol, seed=1)

    positions_a = conf_a.GetConformer().GetPositions()
    positions_b = conf_b.GetConformer().GetPositions()
    assert (positions_a == positions_b).all()


def test_generate_3d_coordinates_raises_on_unembeddable_molecule():
    with pytest.raises(ValueError):
        generate_3d_coordinates(Chem.MolFromSmiles(""), num_conformers=1)
