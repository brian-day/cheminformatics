"""Command-line interface for the cheminformatics platform."""

import typer
from rdkit import Chem

from cheminformatics.descriptors import compute_descriptors, lipinski_violations
from cheminformatics.fingerprints import tanimoto_similarity
from cheminformatics.io import load_molecule

app = typer.Typer(help="Cheminformatics tools for structure-based drug design.")


@app.command()
def describe(molecule: str) -> None:
    """Print physicochemical descriptors and Lipinski rule-of-five violations for a molecule."""
    mol = load_molecule(molecule)
    descriptors = compute_descriptors(mol)

    typer.echo(f"SMILES: {Chem.MolToSmiles(mol)}")
    for key, value in descriptors.to_dict().items():
        typer.echo(f"  {key}: {value}")

    violations = lipinski_violations(descriptors)
    if violations:
        typer.echo(f"Lipinski violations: {', '.join(violations)}")
    else:
        typer.echo("Lipinski violations: none")


@app.command()
def similarity(molecule_a: str, molecule_b: str) -> None:
    """Compute Tanimoto similarity between two molecules (Morgan fingerprints)."""
    mol_a = load_molecule(molecule_a)
    mol_b = load_molecule(molecule_b)
    score = tanimoto_similarity(mol_a, mol_b)
    typer.echo(f"Tanimoto similarity: {score:.4f}")


if __name__ == "__main__":
    app()
