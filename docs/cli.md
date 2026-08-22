# CLI

```bash
uv run chem describe "CC(=O)OC1=CC=CC=C1C(=O)O"   # descriptors + Lipinski check
uv run chem similarity "c1ccccc1" "Cc1ccccc1"      # Tanimoto similarity
```

Molecule arguments accept a SMILES string or a path to a `.mol`, `.sdf`, `.pdb`, or `.smi` file.
