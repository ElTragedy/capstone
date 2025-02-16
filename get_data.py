import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

def smiles_morgan_fingerprint(smiles, radius = 2, n_bits = 2048):
  mol = Chem.MolFromSmiles(smiles)
  if mol:
    return np.array(AllChem.GetMorganFingerprintAsBitVect(mol,radius,nBits=n_bits))
  else:
    return None
  
def smiles_to_graph(smiles):
    mol = Chem.MolFromSmiles(smiles)
    num_nodes = mol.GetNumAtoms()

    node_features = []
    for atom in mol.GetAtoms():
        node_features.append([atom.GetAtomicNum()])

    adj_matrix = Chem.GetAdjacencyMatrix(mol)

    return adj_matrix, node_features

qm9 = pd.read_csv("https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/qm9.csv")

qm9["Morgan_fingerprint"] = qm9["smiles"].apply(lambda x: smiles_morgan_fingerprint(x))
qm9[['adj_matrix', 'node_features']] = qm9.apply(lambda row: pd.Series(smiles_to_graph(row['smiles'])), axis=1)


qm9.to_csv('qm9.csv', index=False)