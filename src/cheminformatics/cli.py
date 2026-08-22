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


@app.command()
def dock(
    receptor_pdb: Path,
    ligand: str,
    output: Path = typer.Option(..., "--output", "-o", help="Output SDF file for ranked docked poses"),
    reference_ligand: str = typer.Option(
        None, help="A known-binding ligand (SMILES or file) to center the search box on"
    ),
    center_x: float = typer.Option(None, help="Box center x (Angstroms); alternative to --reference-ligand"),
    center_y: float = typer.Option(None, help="Box center y (Angstroms)"),
    center_z: float = typer.Option(None, help="Box center z (Angstroms)"),
    box_size: float = typer.Option(20.0, help="Cubic box side length (Angstroms), if not using a reference ligand"),
    num_modes: int = typer.Option(9, help="Number of poses to generate"),
    exhaustiveness: int = typer.Option(8, help="Vina search exhaustiveness"),
    vina_executable: str = typer.Option("vina", help="Path to (or name of) the vina executable"),
) -> None:
    """Dock a ligand against a prepared receptor with AutoDock Vina and write ranked poses to an SDF file.

    Requires the `docking` extra (`uv sync --extra docking`) and the `vina` binary on PATH.
    """
    import tempfile

    from cheminformatics.docking import (
        DockingBox,
        docking_box_from_ligand,
        prepare_ligand_pdbqt,
        prepare_receptor_pdbqt,
        run_vina,
    )

    if reference_ligand is not None:
        reference_mol = generate_3d_coordinates(load_molecule(reference_ligand))
        box = docking_box_from_ligand(reference_mol)
    elif center_x is not None and center_y is not None and center_z is not None:
        box = DockingBox(center=(center_x, center_y, center_z), size=(box_size, box_size, box_size))
    else:
        raise typer.BadParameter("Provide either --reference-ligand or --center-x/--center-y/--center-z")

    ligand_mol = load_molecule(ligand)

    with tempfile.TemporaryDirectory() as tmp_dir:
        receptor_pdbqt = prepare_receptor_pdbqt(receptor_pdb, Path(tmp_dir) / "receptor.pdbqt")
        ligand_pdbqt = prepare_ligand_pdbqt(ligand_mol, Path(tmp_dir) / "ligand.pdbqt")
        results = run_vina(
            receptor_pdbqt,
            ligand_pdbqt,
            box,
            num_modes=num_modes,
            exhaustiveness=exhaustiveness,
            vina_executable=vina_executable,
        )

    writer = Chem.SDWriter(str(output))
    for result in results:
        result.mol.SetProp("vina_affinity", f"{result.affinity:.2f}")
        writer.write(result.mol)
    writer.close()

    typer.echo(f"Wrote {len(results)} docked poses to {output}")
    for i, result in enumerate(results, start=1):
        typer.echo(f"  pose {i}: affinity {result.affinity:.2f} kcal/mol")


if __name__ == "__main__":
    app()
