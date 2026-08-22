"""Loading and writing molecules from common cheminformatics file formats."""

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Mol


def load_molecule(source: str) -> Mol:
    """Load a single molecule from a SMILES string or a file path (.smi, .mol, .sdf, .pdb).

    For .sdf files containing multiple records, the first molecule is returned;
    use `load_molecules` to get all of them.
    """
    path = Path(source)
    if not path.exists():
        mol = Chem.MolFromSmiles(source)
        if mol is None:
            raise ValueError(f"Could not parse '{source}' as a SMILES string or find it as a file")
        return mol

    suffix = path.suffix.lower()
    if suffix in (".mol",):
        mol = Chem.MolFromMolFile(str(path))
    elif suffix == ".sdf":
        mols = load_molecules(path)
        if not mols:
            raise ValueError(f"No parsable molecules found in '{path}'")
        return mols[0]
    elif suffix == ".pdb":
        mol = Chem.MolFromPDBFile(str(path))
    elif suffix in (".smi", ".smiles"):
        with path.open() as f:
            first_line = f.readline().strip().split()[0]
        mol = Chem.MolFromSmiles(first_line)
    else:
        raise ValueError(f"Unsupported file extension: '{suffix}'")

    if mol is None:
        raise ValueError(f"RDKit failed to parse molecule from '{path}'")
    return mol


def load_molecules(path: str | Path) -> list[Mol]:
    """Load all molecules from a multi-record file (currently .sdf)."""
    path = Path(path)
    if path.suffix.lower() != ".sdf":
        return [load_molecule(str(path))]

    supplier = Chem.SDMolSupplier(str(path))
    return [mol for mol in supplier if mol is not None]


def write_molecule(mol: Mol, destination: str | Path) -> None:
    """Write a molecule to a file, inferring the format from the extension (.sdf, .mol, .pdb)."""
    path = Path(destination)
    suffix = path.suffix.lower()

    if suffix in (".mol", ".sdf"):
        Chem.MolToMolFile(mol, str(path))
    elif suffix == ".pdb":
        Chem.MolToPDBFile(mol, str(path))
    else:
        raise ValueError(f"Unsupported output extension: '{suffix}'")
