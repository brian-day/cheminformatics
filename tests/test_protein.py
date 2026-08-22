from pathlib import Path

from cheminformatics.protein import prepare_protein

FIXTURE_DIR = Path(__file__).parent / "fixtures"


def _count_atoms(pdb_path: Path, record: str = "ATOM") -> int:
    with pdb_path.open() as f:
        return sum(1 for line in f if line.startswith(record))


def test_prepare_protein_adds_hydrogens(tmp_path):
    output = prepare_protein(FIXTURE_DIR / "1crn.pdb", tmp_path / "1crn_prepped.pdb")

    heavy_atom_count = _count_atoms(FIXTURE_DIR / "1crn.pdb")
    total_atom_count = _count_atoms(output)

    assert total_atom_count > heavy_atom_count


def test_prepare_protein_removes_heterogens_by_default(tmp_path):
    pdb_with_water = tmp_path / "with_water.pdb"
    pdb_with_water.write_text(
        (FIXTURE_DIR / "1crn.pdb").read_text().replace("END", "")
        + "HETATM  999  O   HOH A 200      10.000  10.000  10.000  1.00  0.00           O\n"
        + "END\n"
    )

    output = prepare_protein(pdb_with_water, tmp_path / "prepped.pdb", keep_water=False)

    assert _count_atoms(output, record="HETATM") == 0
