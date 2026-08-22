"""FastAPI web app: upload or paste molecules, view them alongside their calculated properties."""

import json
import tempfile
from pathlib import Path

import uvicorn
from fastapi import FastAPI, Form, Request, UploadFile
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates
from rdkit import Chem
from rdkit.Chem import Mol
from rdkit.Chem.Draw import rdMolDraw2D

from cheminformatics.conformers import generate_3d_coordinates
from cheminformatics.descriptors import compute_descriptors, lipinski_violations
from cheminformatics.io import load_molecule, load_molecules

app = FastAPI(title="Cheminformatics")
_templates = Jinja2Templates(directory=Path(__file__).parent / "templates")


def render_2d_svg(mol: Mol, size: int = 300) -> str:
    """Render a molecule as a 2D structure depiction (SVG markup)."""
    drawer = rdMolDraw2D.MolDraw2DSVG(size, size)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def render_3d_molblock(mol: Mol) -> str | None:
    """Generate a 3D conformer and return it as a MOL block string, or None if generation fails."""
    try:
        conformer_mol = generate_3d_coordinates(mol)
    except ValueError:
        return None
    return Chem.MolToMolBlock(conformer_mol)


def _analyze_mol(mol: Mol, label: str) -> dict:
    descriptors = compute_descriptors(mol)
    return {
        "input": label,
        "error": None,
        "svg": render_2d_svg(mol),
        "molblock_3d_json": json.dumps(render_3d_molblock(mol)),
        "descriptors": descriptors.to_dict(),
        "violations": lipinski_violations(descriptors),
    }


def _analyze_smiles_text(smiles_text: str) -> list[dict]:
    results = []
    for line in smiles_text.splitlines():
        smiles = line.strip()
        if not smiles:
            continue
        try:
            mol = load_molecule(smiles)
            results.append(_analyze_mol(mol, smiles))
        except ValueError as exc:
            results.append({"input": smiles, "error": str(exc)})
    return results


def _analyze_upload(file: UploadFile) -> list[dict]:
    suffix = Path(file.filename or "").suffix
    with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as tmp:
        tmp.write(file.file.read())
        tmp_path = Path(tmp.name)

    try:
        mols = load_molecules(tmp_path)
    except ValueError as exc:
        return [{"input": file.filename, "error": str(exc)}]
    finally:
        tmp_path.unlink()

    return [_analyze_mol(mol, file.filename or "uploaded file") for mol in mols]


@app.get("/", response_class=HTMLResponse)
def index(request: Request) -> HTMLResponse:
    return _templates.TemplateResponse(request, "index.html", {"results": []})


@app.post("/analyze", response_class=HTMLResponse)
async def analyze(
    request: Request,
    smiles_text: str = Form(""),
    file: UploadFile | None = None,
) -> HTMLResponse:
    results = _analyze_smiles_text(smiles_text)
    if file is not None and file.filename:
        results += _analyze_upload(file)

    return _templates.TemplateResponse(request, "index.html", {"results": results, "smiles_text": smiles_text})


def main() -> None:
    """Run the web app with uvicorn (entry point for the `chem-web` script)."""
    uvicorn.run(app, host="127.0.0.1", port=8080)


if __name__ == "__main__":
    main()
