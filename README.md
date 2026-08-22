# Cheminformatics

A cheminformatics platform for structure-based drug design, built on [RDKit](https://www.rdkit.org/).

Currently implements the core molecule-handling layer: file I/O, physicochemical
descriptors, and fingerprint-based similarity search. Structure-based features
(protein prep, pocket detection, docking) are planned as the next layer.

## Repository structure

```
src/cheminformatics/
    io.py            Reading/writing molecules (SMILES, .mol, .sdf, .pdb, .smi) as RDKit Mol objects
    descriptors.py   Physicochemical descriptors (MW, LogP, TPSA, ...) and Lipinski rule-of-five checks
    fingerprints.py  Morgan fingerprints, Tanimoto similarity, ranked similarity search
    cli.py           Typer CLI (`chem ...`) wiring the above into commands
    web.py           FastAPI app (`chem-web`) wiring the above into a browser UI
    templates/       Jinja2 templates for the web app
examples/            Standalone example scripts demonstrating the library
tests/               pytest tests, one file per module above
```

Each module is a thin, independent layer over RDKit — `io` produces `Mol` objects,
`descriptors` and `fingerprints` consume them, and `cli` exposes both as commands.
As structure-based features are added (protein prep, pocket detection, docking),
they'll follow the same pattern: one module per concern under `src/cheminformatics/`,
with a matching test file and CLI command.

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

## Web app

A small FastAPI app lets you paste SMILES or upload a molecule file and view each
molecule's 2D structure alongside its calculated properties in the browser.

```bash
uv sync --extra web
uv run chem-web
```

Then open http://127.0.0.1:8080 in a browser.

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

## Docs

Docs are built with [MkDocs Material](https://squidfunk.github.io/mkdocs-material/) from
the `docs/` directory. There's no hosted version — run them locally:

```bash
uv sync --group docs
uv run mkdocs serve --dev-addr 127.0.0.1:8090
```

## Roadmap

- Protein structure preparation (cleanup, protonation) and binding pocket detection
- Docking / virtual screening integration (e.g. AutoDock Vina / smina)
- Pose scoring and reporting
