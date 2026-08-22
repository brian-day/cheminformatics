"""Protein structure preparation for docking: fixing up a raw crystal structure PDB.

Crystal structures are rarely usable for docking as downloaded. X-ray diffraction typically
can't resolve hydrogen atoms at all (they scatter X-rays too weakly), so PDB files ship with
none; low electron-density regions (flexible loops, disordered side chains) are often missing
one or more atoms, or even whole residues; and the file usually still contains waters, bound
ions, or crystallization buffer molecules that aren't part of the biological receptor you want
to dock against. None of this is optional cleanup — a docking engine needs a complete, sensibly
protonated structure to compute reasonable interaction energies, and skipping any of these steps
introduces artifacts (missing atoms create clashes or holes; the wrong protonation state changes
which functional groups are charged, which changes electrostatic scoring).

This module wraps `PDBFixer` (built on OpenMM) to perform the standard prep sequence: replace
non-standard residues, fill in missing residues/atoms, add hydrogens at a target pH, and strip
out heterogens (waters, ions, ligands) that aren't part of the receptor itself.
"""

from pathlib import Path

from openmm.app import PDBFile
from pdbfixer import PDBFixer


def prepare_protein(
    input_pdb: str | Path,
    output_pdb: str | Path,
    ph: float = 7.4,
    keep_water: bool = False,
) -> Path:
    """Clean a raw PDB structure for docking and write the result to `output_pdb`.

    Applies, in order:

    1. **Non-standard residue replacement** — swaps modified residues (e.g. a chemically
       altered amino acid used to aid crystallization) for their standard counterparts.
    2. **Missing residue/atom completion** — fills in residues or individual atoms absent
       from the structure due to weak electron density, so the chain is structurally complete.
    3. **Hydrogen addition at `ph`** — adds hydrogens with protonation states appropriate for
       the given pH (this determines, e.g., whether histidine or a carboxylic acid side chain
       is charged — physiological pH 7.4 is the usual default for drug-target docking).
    4. **Heterogen removal** — strips waters, ions, and other non-receptor molecules, keeping
       water only if `keep_water` is set (occasionally a structural water mediates a binding
       interaction you want to preserve for docking).

    Returns the path the cleaned structure was written to.
    """
    fixer = PDBFixer(filename=str(input_pdb))

    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()

    fixer.removeHeterogens(keepWater=keep_water)

    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)

    output_pdb = Path(output_pdb)
    with output_pdb.open("w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    return output_pdb
