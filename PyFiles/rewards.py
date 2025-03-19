import tensorflow as tf
import numpy as np
from rdkit import Chem
from tensorflow.keras import optimizers

def validity_reward(mol):
    if mol is None:
      print("Molecule is None")
      return -5
    # Sanitize
    try:
      Chem.SanitizeMol(mol)
    except Exception:
      print("Sanitization failed")
      return -5
    # Check Kekulization
    try:
      Chem.Kekulize(mol, clearAromaticFlags = True)
    except Exception:
      print("Kekulization failed")
      return -5
    # Check Valency
    for atom in mol.GetAtoms():
      explicit_valence = atom.GetExplicitValence()
      if explicit_valence > atom.GetTotalValence():
        print("Valency failed")
        return -5
    return 5

def uniqueness_reward(gen_smiles, curr_mol):
  curr_smiles = Chem.MolToSmiles(curr_mol)
  if validity_reward(curr_mol) == 0:
      print("Validity failed")
      return -5
  occurrence_count = gen_smiles.count(curr_smiles)
  total_molecules = len(gen_smiles) + 1  

  uniqueness_score = 1 - (occurrence_count / total_molecules)

  uniqueness_score = max(0, min(uniqueness_score, 1))

  return uniqueness_score
  
def novelty_reward(curr_mol, train_smiles, generated_smiles):
    unique_train = set(train_smiles)
    if validity_reward(curr_mol) == 0:
        print("Validity failed")
        return -5
    elif uniqueness_reward(generated_smiles, curr_mol) == 0:
        print("Uniqueness failed")
        return -3
    curr_smiles = Chem.MolToSmiles(curr_mol)
    if curr_smiles in unique_train:
        print("Novelty failed")
        return -10
    else:
        return 10