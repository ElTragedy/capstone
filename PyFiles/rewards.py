from rdkit import Chem
from rdkit.Chem import QED
import tensorflow as tf
import warnings
warnings.filterwarnings('ignore')

def validity_reward(mol):
    if mol is None:
        print("Molecule is None")
        return -0.5

    try:
        Chem.SanitizeMol(mol)
    except Exception:
        print("Sanitization failed")
        return -0.5

    valence_score = 0
    for atom in mol.GetAtoms():
        explicit_valence = atom.GetExplicitValence()
        max_valence = atom.GetTotalValence()
        valence_score += max_valence - explicit_valence  

    return max(0, 5 - abs(valence_score))  

def stray_hydros_reward(mol):
    if mol is None:
        return -0.5
    elif any(atom.GetAtomicNum() == 1 and atom.GetDegree() == 0 for atom in mol.GetAtoms()):
        return -0.5
    else:
        return 0.5

def uniqueness_reward(gen_smiles, curr_mol):
  curr_smiles = Chem.MolToSmiles(curr_mol)
  if validity_reward(curr_mol) == 0:
      print("Validity failed")
      return -0.5
  occurrence_count = gen_smiles.count(curr_smiles)
  total_molecules = len(gen_smiles) + 1  

  uniqueness_score = 1 - (occurrence_count / total_molecules)

  uniqueness_score = max(0, min(uniqueness_score, 1))

  return uniqueness_score
  
def novelty_reward(curr_mol, train_smiles, generated_smiles):
    unique_train = set(train_smiles)
    if validity_reward(curr_mol) == 0:
        print("Validity failed")
        return -0.3
    elif uniqueness_reward(generated_smiles, curr_mol) == 0:
        print("Uniqueness failed")
        return -0.3
    curr_smiles = Chem.MolToSmiles(curr_mol)
    if curr_smiles in unique_train:
        print("Novelty failed")
        return -0.7
    else:
        return 1
    
def drug_like_reward(mol):
    if validity_reward(mol):
        drug_score = QED.qed(mol)
        return drug_score
    else: 
        return -0.3
