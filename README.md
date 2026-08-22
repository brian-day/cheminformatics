# Cheminformatics

A cheminformatics platform for structure-based drug design, built on [RDKit](https://www.rdkit.org/).

Currently implements the core molecule-handling layer: file I/O, physicochemical
descriptors, and fingerprint-based similarity search. Structure-based features
(protein prep, pocket detection, docking) are planned as the next layer.

## Setup

```bash
uv sync
```

## CLI

```bash
uv run chem describe "CC(=O)OC1=CC=CC=C1C(=O)O"   # descriptors + Lipinski check
uv run chem similarity "c1ccccc1" "Cc1ccccc1"      # Tanimoto similarity
```

Molecule arguments accept a SMILES string or a path to a `.mol`, `.sdf`, `.pdb`, or `.smi` file.

## Library

```python
from cheminformatics.io import load_molecule
from cheminformatics.descriptors import compute_descriptors, lipinski_violations
from cheminformatics.fingerprints import tanimoto_similarity, similarity_search

mol = load_molecule("CC(=O)OC1=CC=CC=C1C(=O)O")
descriptors = compute_descriptors(mol)
lipinski_violations(descriptors)
```

## Tests

```bash
uv run pytest
```

## Roadmap

- Protein structure preparation (cleanup, protonation) and binding pocket detection
- Docking / virtual screening integration (e.g. AutoDock Vina / smina)
- Pose scoring and reporting
