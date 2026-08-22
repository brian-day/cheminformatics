import os
import stat
from pathlib import Path

import pytest
from rdkit import Chem

from cheminformatics.conformers import generate_3d_coordinates
from cheminformatics.docking import (
    DockingBox,
    docking_box_from_ligand,
    prepare_ligand_pdbqt,
    prepare_receptor_pdbqt,
    run_vina,
)

FIXTURE_DIR = Path(__file__).parent / "fixtures"
ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"

_FAKE_VINA = """\
#!/usr/bin/env python3
import sys

args = dict(zip(sys.argv[1::2], sys.argv[2::2]))
with open(args["--ligand"]) as f:
    lines = f.read().splitlines()

remarks = "\\n".join(l for l in lines if l.startswith(("REMARK SMILES", "REMARK H PARENT")))
body = "\\n".join(l for l in lines if not l.startswith(("REMARK SMILES", "REMARK H PARENT")))

with open(args["--out"], "w") as f:
    for i, affinity in enumerate([-6.3, -5.8], start=1):
        f.write(f"MODEL {i:8d}\\n")
        f.write(remarks + "\\n")
        f.write(f"REMARK VINA RESULT: {affinity:10.1f} {0.0:10.3f} {0.0:10.3f}\\n")
        f.write(body + "\\n")
        f.write("ENDMDL\\n")
"""


@pytest.fixture
def fake_vina(tmp_path, monkeypatch):
    script = tmp_path / "vina"
    script.write_text(_FAKE_VINA)
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", f"{tmp_path}{os.pathsep}{os.environ['PATH']}")
    return script


def test_prepare_receptor_pdbqt(tmp_path):
    output = prepare_receptor_pdbqt(FIXTURE_DIR / "1crn.pdb", tmp_path / "receptor.pdbqt")

    assert output.exists()
    assert "ROOT" not in output.read_text()  # rigid receptor: no torsion tree


def test_prepare_ligand_pdbqt_roundtrips_via_smiles_remark(tmp_path):
    mol = Chem.MolFromSmiles(ASPIRIN_SMILES)
    output = prepare_ligand_pdbqt(mol, tmp_path / "ligand.pdbqt")

    contents = output.read_text()
    assert "REMARK SMILES CC(=O)Oc1ccccc1C(=O)O" in contents
    assert "ROOT" in contents  # flexible ligand: has a torsion tree


def test_docking_box_from_ligand_centers_on_bounding_box():
    mol = generate_3d_coordinates(Chem.MolFromSmiles("CCO"))
    box = docking_box_from_ligand(mol, padding=5.0)

    positions = mol.GetConformer().GetPositions()
    assert box.center == pytest.approx(tuple(positions.mean(axis=0)), abs=2.0)
    assert all(s > 10.0 for s in box.size)  # at least 2x padding


def test_run_vina_returns_ranked_poses(tmp_path, fake_vina):
    receptor_pdbqt = prepare_receptor_pdbqt(FIXTURE_DIR / "1crn.pdb", tmp_path / "receptor.pdbqt")
    ligand_pdbqt = prepare_ligand_pdbqt(Chem.MolFromSmiles(ASPIRIN_SMILES), tmp_path / "ligand.pdbqt")

    results = run_vina(receptor_pdbqt, ligand_pdbqt, DockingBox(center=(0.0, 0.0, 0.0)))

    assert [r.affinity for r in results] == [-6.3, -5.8]
    assert Chem.MolToSmiles(Chem.RemoveHs(results[0].mol)) == ASPIRIN_SMILES


def test_run_vina_raises_if_executable_missing(tmp_path):
    receptor_pdbqt = prepare_receptor_pdbqt(FIXTURE_DIR / "1crn.pdb", tmp_path / "receptor.pdbqt")
    ligand_pdbqt = prepare_ligand_pdbqt(Chem.MolFromSmiles("CCO"), tmp_path / "ligand.pdbqt")

    with pytest.raises(RuntimeError, match="Could not find"):
        run_vina(
            receptor_pdbqt,
            ligand_pdbqt,
            DockingBox(center=(0.0, 0.0, 0.0)),
            vina_executable="definitely_not_a_real_vina_binary",
        )
