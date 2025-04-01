import tensorflow as tf
import numpy as np
from rdkit import Chem

def adjacency_matrix_to_mol(adj, nodes):
  mol = Chem.RWMol()
  atom_types = np.unique(nodes)
  matrix = adj
  atom_map = {i: mol.AddAtom(Chem.Atom(int(atom_types[i]))) for i in range(len(atom_types))}

  if isinstance(matrix, tf.Tensor):
      matrix = matrix.numpy()
  if isinstance(atom_types, tf.Tensor):
      atom_types = atom_types.numpy()

  matrix = np.array(matrix)
  atoms = np.array(nodes)
  atom_types = np.array(atom_types)

  if len(matrix) != len(atoms):
    raise ValueError("NUMBER OF ATOM TYPES DOES NOT MATCH MATRIX DIMENSIONS")

  # Step 1: Add Atoms
  for i, atom_num in enumerate(atom_types):
    atom = Chem.Atom(int(atom_num))
    mol_idx = mol.AddAtom(atom)
    atom_map[i] = mol_idx

  # Step 2: Add Bonds
  for i in range(len(matrix)):
    for j in range(i + 1, len(matrix)):
        if j not in atom_map:
            continue

        value = int(np.argmax(matrix[i, j]))
        value = min(max(value, 0), 4)

        if value != 0:
            bond_type = {
                1: Chem.BondType.SINGLE,
                2: Chem.BondType.DOUBLE,
                3: Chem.BondType.TRIPLE,
                4: Chem.BondType.AROMATIC
            }.get(value, None)

            if bond_type is None:
                raise ValueError(f"INVALID BOND TYPE DETECTED: {value}")

            print(f"Adding bond: {i}-{j} Type: {bond_type}")
            mol.AddBond(atom_map[i], atom_map[j], bond_type)

  return mol