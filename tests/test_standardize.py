from rdkit import Chem

from cheminformatics.standardize import standardize_molecule


def test_standardize_molecule_strips_salt():
    mol = Chem.MolFromSmiles("CN(C)C(=N)NC(=N)N.Cl")  # metformin hydrochloride
    std = standardize_molecule(mol)

    assert "." not in Chem.MolToSmiles(std)


def test_standardize_molecule_neutralizes_charge():
    mol = Chem.MolFromSmiles("CC(=O)[O-]")  # acetate anion
    std = standardize_molecule(mol)

    assert Chem.GetFormalCharge(std) == 0
    assert Chem.MolToSmiles(std) == Chem.MolToSmiles(Chem.MolFromSmiles("CC(=O)O"))


def test_standardize_molecule_canonicalizes_tautomers():
    keto = standardize_molecule(Chem.MolFromSmiles("CC(=O)C"))
    enol = standardize_molecule(Chem.MolFromSmiles("CC(O)=C"))

    assert Chem.MolToSmiles(keto) == Chem.MolToSmiles(enol)
