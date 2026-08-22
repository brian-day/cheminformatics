"""Molecular fingerprinting and similarity search."""

from rdkit import DataStructs
from rdkit.Chem import Mol
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

_MORGAN_GENERATOR = GetMorganGenerator(radius=2, fpSize=2048)


def morgan_fingerprint(mol: Mol) -> DataStructs.ExplicitBitVect:
    """Compute a Morgan (ECFP-like) fingerprint, radius 2, 2048 bits.

    A Morgan fingerprint encodes a molecule's structure by hashing the
    circular atom neighborhoods around each atom (up to the given radius)
    into a fixed-length bit vector. Two molecules that share substructures
    tend to set many of the same bits, which is what makes the fingerprint
    useful for similarity comparisons.
    """
    return _MORGAN_GENERATOR.GetFingerprint(mol)


def tanimoto_similarity(mol_a: Mol, mol_b: Mol) -> float:
    """Tanimoto similarity between two molecules' Morgan fingerprints.

    Tanimoto similarity is the ratio of bits set in both fingerprints
    (intersection) to bits set in either fingerprint (union), giving a
    score from 0 (no shared bits) to 1 (identical fingerprints).
    """
    fp_a = morgan_fingerprint(mol_a)
    fp_b = morgan_fingerprint(mol_b)
    return DataStructs.TanimotoSimilarity(fp_a, fp_b)


def similarity_search(query: Mol, candidates: list[Mol], threshold: float = 0.7) -> list[tuple[int, float]]:
    """Rank candidate molecules by Tanimoto similarity to a query, filtering by threshold.

    Returns a list of (index_into_candidates, similarity) sorted by descending similarity.
    """
    query_fp = morgan_fingerprint(query)
    candidate_fps = [morgan_fingerprint(mol) for mol in candidates]
    similarities = DataStructs.BulkTanimotoSimilarity(query_fp, candidate_fps)

    hits = [(i, sim) for i, sim in enumerate(similarities) if sim >= threshold]
    hits.sort(key=lambda pair: pair[1], reverse=True)
    return hits
