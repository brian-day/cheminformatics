"""Docking: PDBQT preparation, a search box helper, and running/parsing AutoDock Vina.

This module deliberately doesn't reimplement docking search or scoring — that's a
substantial, actively-researched optimization problem, and AutoDock Vina already does it
well. Its job is the surrounding plumbing: converting a receptor and ligand into the
PDBQT format Vina expects, picking a reasonable search box, shelling out to the `vina`
binary, and parsing its output (binding affinities and docked poses) back into RDKit
`Mol` objects you can keep working with.

`vina` itself is not a Python dependency of this project — it's a separate binary you
install yourself and must have on PATH (see the README for install options). This mirrors
how the rest of the docking pipeline is scoped: heavy search/optimization work happens in
an external, purpose-built program, and this module prepares its inputs and reads its
outputs.

Ligand PDBQT preparation uses Meeko rather than a generic converter (e.g. Open Babel)
because of a real interoperability gap: AutoDock's PDBQT format uses "united atom" typing
for nonpolar hydrogens — it merges them into their carbon rather than listing them as
atoms — so a bare PDBQT pose can't be unambiguously converted back into a correct
molecule (a generic reader ends up with open valences instead of the original C-H bonds).
Meeko works around this by embedding a SMILES + atom-index mapping as a REMARK line when
it writes the ligand PDBQT, which is what lets a docked pose be reconstructed as a
correct, fully-protonated RDKit `Mol` afterward. The receptor side doesn't have this
problem (it's never converted back into a `Mol`), so it uses Open Babel, which is simpler
for a plain rigid conversion.
"""

import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from meeko import MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy, RDKitMolCreate
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import Mol

from cheminformatics.conformers import generate_3d_coordinates


@dataclass
class DockingBox:
    """An axis-aligned search box, in Angstroms, that Vina should search within."""

    center: tuple[float, float, float]
    size: tuple[float, float, float] = (20.0, 20.0, 20.0)


@dataclass
class DockingResult:
    """A single docked ligand pose and its predicted binding affinity (kcal/mol).

    More negative affinity means a more energetically favorable (predicted) binding pose;
    this is Vina's empirical scoring function, not a physically rigorous free energy.
    """

    mol: Mol
    affinity: float


def docking_box_from_ligand(ligand_mol: Mol, padding: float = 8.0) -> DockingBox:
    """Compute a search box centered on a reference ligand's 3D coordinates.

    A common way to define where Vina should search is to center the box on a ligand
    already known to bind the target (e.g. a co-crystallized inhibitor) and pad its
    bounding box by a margin. This gives the search room to explore nearby poses without
    having to search the entire protein surface, which would be both slower and more
    likely to find spurious high-scoring poses outside the real binding site.
    """
    positions = ligand_mol.GetConformer().GetPositions()
    mins = positions.min(axis=0)
    maxs = positions.max(axis=0)
    center = tuple(float(v) for v in (mins + maxs) / 2)
    size = tuple(float(v) for v in (maxs - mins) + 2 * padding)
    return DockingBox(center=center, size=size)


def prepare_receptor_pdbqt(receptor_pdb: str | Path, output_pdbqt: str | Path) -> Path:
    """Convert a prepared receptor PDB (see `protein.prepare_protein`) to a rigid PDBQT file."""
    output_pdbqt = Path(output_pdbqt)
    obmol = next(pybel.readfile("pdb", str(receptor_pdb)))
    obmol.calccharges("gasteiger")
    obmol.write("pdbqt", str(output_pdbqt), overwrite=True, opt={"r": None})
    return output_pdbqt


def prepare_ligand_pdbqt(mol: Mol, output_pdbqt: str | Path) -> Path:
    """Convert a ligand `Mol` to a flexible PDBQT file, generating a 3D conformer first if needed."""
    if mol.GetNumConformers() == 0 or not mol.GetConformer().Is3D():
        mol = generate_3d_coordinates(mol)
    mol = Chem.AddHs(mol, addCoords=True)

    molecule_setups = MoleculePreparation().prepare(mol)
    pdbqt_string, success, error_message = PDBQTWriterLegacy.write_string(molecule_setups[0])
    if not success:
        raise ValueError(f"Could not prepare ligand PDBQT: {error_message}")

    output_pdbqt = Path(output_pdbqt)
    output_pdbqt.write_text(pdbqt_string)
    return output_pdbqt


def run_vina(
    receptor_pdbqt: str | Path,
    ligand_pdbqt: str | Path,
    box: DockingBox,
    num_modes: int = 9,
    exhaustiveness: int = 8,
    vina_executable: str = "vina",
) -> list[DockingResult]:
    """Dock a ligand against a receptor with AutoDock Vina, returning poses ranked best-first.

    Runs `vina_executable` as a subprocess with the given search box, then reads back the
    output PDBQT it writes (which Vina annotates per-pose with a `REMARK VINA RESULT:`
    affinity line) and reconstructs each pose as a `DockingResult`.

    Raises `RuntimeError` if the executable can't be found on PATH, or if Vina itself fails.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        output_pdbqt = Path(tmp_dir) / "docked.pdbqt"
        command = [
            vina_executable,
            "--receptor",
            str(receptor_pdbqt),
            "--ligand",
            str(ligand_pdbqt),
            "--center_x",
            str(box.center[0]),
            "--center_y",
            str(box.center[1]),
            "--center_z",
            str(box.center[2]),
            "--size_x",
            str(box.size[0]),
            "--size_y",
            str(box.size[1]),
            "--size_z",
            str(box.size[2]),
            "--num_modes",
            str(num_modes),
            "--exhaustiveness",
            str(exhaustiveness),
            "--out",
            str(output_pdbqt),
        ]
        try:
            subprocess.run(command, capture_output=True, text=True, check=True)
        except FileNotFoundError as exc:
            raise RuntimeError(
                f"Could not find the '{vina_executable}' executable. Install AutoDock Vina and make "
                "sure it's on PATH (see the README)."
            ) from exc
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(f"vina failed:\n{exc.stderr}") from exc

        pdbqt_mol = PDBQTMolecule.from_file(output_pdbqt, skip_typing=True)
        affinities = []
        for _ in pdbqt_mol:
            affinities.append(pdbqt_mol.score)

    multi_pose_mol = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)[0]
    return [
        DockingResult(mol=Chem.Mol(multi_pose_mol, confId=conf_id), affinity=affinity)
        for conf_id, affinity in enumerate(affinities)
    ]
