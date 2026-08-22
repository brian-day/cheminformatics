# Cheminformatics

A cheminformatics platform for structure-based drug design, built on [RDKit](https://www.rdkit.org/).

Currently implements the core molecule-handling layer (file I/O, physicochemical
descriptors, fingerprint-based similarity search, 3D conformer generation) plus the
first structure-based feature: protein preparation for docking. Docking prep/execution
and molecule generation are next.

## Repository structure

```
src/cheminformatics/
    io.py            Reading/writing molecules (SMILES, .mol, .sdf, .pdb, .smi) as RDKit Mol objects
    standardize.py   Molecule standardization (salt stripping, charge neutralization, tautomer canonicalization)
    descriptors.py   Physicochemical descriptors (MW, LogP, TPSA, ...) and Lipinski rule-of-five checks
    fingerprints.py  Morgan fingerprints, Tanimoto similarity, ranked similarity search
    conformers.py    3D conformer generation (ETKDG embedding + MMFF94/UFF minimization)
    protein.py       Protein structure prep for docking (PDBFixer: missing atoms, hydrogens, heterogen removal)
    cli.py           Typer CLI (`chem ...`) wiring the above into commands
    web.py           FastAPI app (`chem-web`) wiring the above into a browser UI
    templates/       Jinja2 templates for the web app
examples/            Standalone example scripts demonstrating the library
tests/               pytest tests, one file per module above
```

Each module is a thin, independent layer over RDKit (or, for `protein.py`, over PDBFixer/OpenMM)
— `io` produces `Mol` objects, `descriptors` and `fingerprints` consume them, and `cli` exposes
each as a command. As more structure-based features are added (docking, molecule generation),
they'll follow the same pattern: one module per concern under `src/cheminformatics/`, with a
matching test file and CLI command.

## Setup

```bash
uv sync
```

## CLI

```bash
uv run chem describe "CC(=O)OC1=CC=CC=C1C(=O)O"           # descriptors + Lipinski check
uv run chem similarity "c1ccccc1" "Cc1ccccc1"              # Tanimoto similarity
uv run chem conformer "CCO" -o ethanol.sdf                 # 3D conformer (ETKDG + MMFF94/UFF)
uv run chem prep-protein receptor.pdb -o receptor_prepped.pdb  # docking-ready protein (requires the `protein` extra)
```

Molecule arguments accept a SMILES string or a path to a `.mol`, `.sdf`, `.pdb`, or `.smi` file.

## Web app

A small FastAPI app lets you paste SMILES or upload a molecule file and view each
molecule's 2D and 3D structure (toggle button per molecule, rendered with
[3Dmol.js](https://3dmol.org/)) alongside its calculated properties in the browser.

```bash
uv sync --extra web
uv run chem-web
```

Then open http://127.0.0.1:8080 in a browser.

## Protein prep

Cleans a raw PDB structure for docking — fills in missing residues/atoms, adds hydrogens
at a target pH, and strips waters/heterogens — using [PDBFixer](https://github.com/openmm/pdbfixer).

```bash
uv sync --extra protein
uv run chem prep-protein receptor.pdb -o receptor_prepped.pdb
```

## Library

```python
from cheminformatics.io import load_molecule
from cheminformatics.descriptors import compute_descriptors, lipinski_violations
from cheminformatics.fingerprints import tanimoto_similarity, similarity_search
from cheminformatics.conformers import generate_3d_coordinates

mol = load_molecule("CC(=O)OC1=CC=CC=C1C(=O)O")  # standardized automatically (see standardize.py)
descriptors = compute_descriptors(mol)
lipinski_violations(descriptors)

conformer_mol = generate_3d_coordinates(mol)  # lowest-energy 3D conformer
```

## Tests

```bash
uv run pytest
```

## Docs

Docs are built with [MkDocs Material](https://squidfunk.github.io/mkdocs-material/) from
the `docs/` directory. There's no hosted version — run them locally:

```bash
uv sync --group docs
uv run mkdocs serve --dev-addr 127.0.0.1:8090
```

## Roadmap

- Docking prep and execution: receptor/ligand → PDBQT, search box definition, subprocess to
  AutoDock Vina, and parsing scores/poses back into RDKit `Mol` objects
- Combinatorial candidate generation: scaffold + substituent enumeration, filtered through the
  existing standardization/descriptor/Lipinski pipeline
- Binding pocket detection and pose scoring/reporting
