# Cheminformatics

A cheminformatics platform for structure-based drug design, built on [RDKit](https://www.rdkit.org/).

Currently implements the core molecule-handling layer: file I/O, physicochemical
descriptors, and fingerprint-based similarity search. Structure-based features
(protein prep, pocket detection, docking) are planned as the next layer.

## Setup

```bash
uv sync
```

## Library

```python
from cheminformatics.io import load_molecule
from cheminformatics.descriptors import compute_descriptors, lipinski_violations
from cheminformatics.fingerprints import tanimoto_similarity, similarity_search

mol = load_molecule("CC(=O)OC1=CC=CC=C1C(=O)O")
descriptors = compute_descriptors(mol)
lipinski_violations(descriptors)
```

See the [CLI](cli.md) and [Web app](web-app.md) pages for the other ways to use this
library, or the API reference for full details on each module.
