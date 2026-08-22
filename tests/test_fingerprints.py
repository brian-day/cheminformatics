from rdkit import Chem

from cheminformatics.fingerprints import similarity_search, tanimoto_similarity

BENZENE = "c1ccccc1"
TOLUENE = "Cc1ccccc1"
METHANE = "C"


def test_tanimoto_similarity_identical_molecules():
    mol = Chem.MolFromSmiles(BENZENE)
    assert tanimoto_similarity(mol, mol) == 1.0


def test_tanimoto_similarity_related_molecules_above_unrelated():
    benzene = Chem.MolFromSmiles(BENZENE)
    toluene = Chem.MolFromSmiles(TOLUENE)
    methane = Chem.MolFromSmiles(METHANE)

    related_score = tanimoto_similarity(benzene, toluene)
    unrelated_score = tanimoto_similarity(benzene, methane)
    assert related_score > unrelated_score


def test_similarity_search_ranks_and_filters():
    query = Chem.MolFromSmiles(BENZENE)
    candidates = [Chem.MolFromSmiles(s) for s in (TOLUENE, METHANE)]

    hits = similarity_search(query, candidates, threshold=0.0)
    assert [index for index, _ in hits] == [0, 1]
