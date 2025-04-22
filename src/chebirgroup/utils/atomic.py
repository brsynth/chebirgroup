from rdkit import Chem

def is_biochemical(mol):
    for atom in mol.GetAtoms():
        # https://www.ncbi.nlm.nih.gov/books/NBK26883/
        if atom.GetAtomicNum() not in [0, 1, 6, 7, 8, 11, 12, 15, 16, 17, 19, 20]:
            return False
    return True

def is_candidate(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return smiles.count("*") > 0 and mol.GetNumHeavyAtoms() > 5 and is_biochemical(mol=mol)
