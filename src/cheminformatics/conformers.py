"""3D conformer generation: embedding a plausible starting geometry and relaxing it to a nearby energy
minimum with a force field.

A 2D structure (a molecular graph) doesn't tell you where atoms actually sit in space. Getting there is a
two-stage process:

1. **Embedding (ETKDG)** — RDKit's ETKDG algorithm turns the graph into an initial 3D guess. It's a
   distance-geometry method (it satisfies known bond-length/angle distance bounds) augmented with torsion
   angle preferences mined from real crystal structures, so the starting guess is *plausible* but not
   necessarily energy-minimized.
2. **Optimization (force field)** — a force field (MMFF94, or UFF as a fallback when MMFF94 lacks
   parameters for some atom type) computes a potential energy from real physical terms — bond stretching,
   angle bending, torsional strain, van der Waals and electrostatic interactions — and locally minimizes
   that energy, relaxing the embedded guess toward a nearby low-energy geometry.

Because embedding is stochastic and torsion preferences are only statistical, a single embed-and-minimize
attempt can land in a mediocre local minimum. Generating several candidate conformers and keeping the
lowest-energy one after minimization gives a better chance of finding a geometry close to the true
minimum-energy conformation.
"""

from rdkit.Chem import AllChem, Mol

_DEFAULT_NUM_CONFORMERS = 10
_DEFAULT_SEED = 0xF00D


def generate_3d_coordinates(mol: Mol, num_conformers: int = _DEFAULT_NUM_CONFORMERS, seed: int = _DEFAULT_SEED) -> Mol:
    """Generate a 3D conformer for a molecule, returning the lowest-energy result.

    Adds explicit hydrogens (required for realistic 3D geometry and force field calculations), embeds
    `num_conformers` candidate 3D structures with ETKDG, minimizes each with the MMFF94 force field
    (falling back to UFF for atom types MMFF94 doesn't cover), and returns a copy of the molecule holding
    only the single lowest-energy conformer.

    Raises `ValueError` if no conformer could be embedded (can happen for heavily constrained or unusual
    structures, e.g. large macrocycles).
    """
    mol = AllChem.AddHs(mol)  # type: ignore[attr-defined]

    params = AllChem.ETKDGv3()  # type: ignore[attr-defined]
    params.randomSeed = seed
    conformer_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, params=params)  # type: ignore[attr-defined]
    if not conformer_ids:
        raise ValueError("Could not embed any 3D conformer for this molecule")

    mmff_properties = AllChem.MMFFGetMoleculeProperties(mol)  # type: ignore[attr-defined]
    energies: dict[int, float] = {}
    for conformer_id in conformer_ids:
        force_field = (
            AllChem.MMFFGetMoleculeForceField(mol, mmff_properties, confId=conformer_id)  # type: ignore[attr-defined]
            if mmff_properties is not None
            else None
        )
        if force_field is None:
            force_field = AllChem.UFFGetMoleculeForceField(mol, confId=conformer_id)  # type: ignore[attr-defined]
        force_field.Minimize()
        energies[conformer_id] = force_field.CalcEnergy()

    best_conformer_id = min(energies, key=lambda conformer_id: energies[conformer_id])
    return Mol(mol, confId=best_conformer_id)
