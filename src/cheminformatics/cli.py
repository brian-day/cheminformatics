"""Command-line interface for the cheminformatics platform."""

from pathlib import Path

import typer
from rdkit import Chem

from cheminformatics.conformers import generate_3d_coordinates
from cheminformatics.descriptors import compute_descriptors, lipinski_violations
from cheminformatics.fingerprints import tanimoto_similarity
from cheminformatics.io import load_molecule, write_molecule

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


@app.command()
def conformer(
    molecule: str,
    output: Path = typer.Option(..., "--output", "-o", help="Output file (.mol, .sdf, .pdb)"),
    num_conformers: int = typer.Option(10, help="Number of candidate conformers to embed and minimize"),
    seed: int = typer.Option(0xF00D, help="Random seed for conformer embedding"),
) -> None:
    """Generate a 3D conformer for a molecule (ETKDG + MMFF94/UFF minimization) and write it to a file."""
    mol = load_molecule(molecule)
    conf_mol = generate_3d_coordinates(mol, num_conformers=num_conformers, seed=seed)
    write_molecule(conf_mol, output)
    typer.echo(f"Wrote 3D conformer to {output}")


@app.command()
def prep_protein(
    input_pdb: Path,
    output: Path = typer.Option(..., "--output", "-o", help="Output PDB file"),
    ph: float = typer.Option(7.4, help="pH used to assign hydrogen protonation states"),
    keep_water: bool = typer.Option(False, help="Keep crystallographic waters instead of stripping them"),
) -> None:
    """Clean a raw PDB structure for docking (missing atoms, hydrogens, heterogen removal).

    Requires the `protein` extra: `uv sync --extra protein`.
    """
    from cheminformatics.protein import prepare_protein

    prepare_protein(input_pdb, output, ph=ph, keep_water=keep_water)
    typer.echo(f"Wrote prepared protein to {output}")


if __name__ == "__main__":
    app()
