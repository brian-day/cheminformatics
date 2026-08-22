"""Compute physicochemical properties for a set of molecules."""

from rdkit import Chem

from cheminformatics.descriptors import compute_descriptors, lipinski_violations
from cheminformatics.io import load_molecule

SMILES = [
    "CCO",  # ethanol
    "c1ccccc1",  # benzene
    "CC(=O)Oc1ccccc1C(=O)O",  # aspirin
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # caffeine
]


def main() -> None:
    for smiles in SMILES:
        mol = load_molecule(smiles)
        descriptors = compute_descriptors(mol)

        print(f"SMILES: {Chem.MolToSmiles(mol)}")
        for key, value in descriptors.to_dict().items():
            print(f"  {key}: {value}")

        violations = lipinski_violations(descriptors)
        print(f"  lipinski_violations: {', '.join(violations) if violations else 'none'}")
        print()


if __name__ == "__main__":
    main()
