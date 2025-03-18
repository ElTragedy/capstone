import tensorflow as tf
import numpy as np
from rdkit import Chem
from tensorflow.keras import optimizers

def validity_reward(mol):
    if mol is None:
      return 0
    # Sanitize
    try:
      Chem.SanitizeMol(mol)
    except Exception:
      return 0
    # Check Kekulization
    try:
      Chem.Kekulize(mol, clearAromaticFlags = True)
    except Exception:
      return 0
    # Check Valency
    for atom in mol.GetAtoms():
      explicit_valence = atom.GetExplicitValence()
      if explicit_valence > atom.GetTotalValence():
        return 0
    return 2

def uniqueness_reward(gen_smiles, curr_mol):
  curr_smiles = Chem.MolToSmiles(curr_mol)
  if validity_reward(curr_mol) == 0:
      return 0
  if curr_smiles in set(gen_smiles):
      return 0
  else:
      return 1
  
def novelty_reward(curr_mol, train_smiles, generated_smiles):
    unique_train = set(train_smiles)
    if validity_reward(curr_mol) == 0:
        return 0
    elif uniqueness_reward(generated_smiles, curr_mol) == 0:
        return 0
    curr_smiles = Chem.MolToSmiles(curr_mol)
    if curr_smiles in unique_train:
        return 0
    else:
        return 3