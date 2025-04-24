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

  mol = Chem.RemoveHs(mol) # Remove extraneous hydrogens

  return mol

"""
import tensorflow as tf
import numpy as np
from rdkit import Chem

def adjacency_matrix_to_mol(adj, nodes):
    # Convert TensorFlow tensors to NumPy arrays if necessary
    if isinstance(adj, tf.Tensor):
        adj = adj.numpy()
    if isinstance(nodes, tf.Tensor):
        nodes = nodes.numpy()

    adj = np.array(adj)
    nodes = np.array(nodes)

    # Validate dimensions
    if adj.shape[0] != adj.shape[1] or adj.shape[0] != len(nodes):
        raise ValueError("Adjacency matrix dimensions must match the number of nodes")

    # Initialize a read-write molecule
    mol = Chem.RWMol()

    # Step 1: Add atoms
    atom_map = {}
    for i, atom_num in enumerate(nodes):
        atom = Chem.Atom(int(atom_num))
        mol_idx = mol.AddAtom(atom)
        atom_map[i] = mol_idx

    # Step 2: Add bonds (assuming binary adjacency matrix: 0 = no bond, 1 = single bond)
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            if adj[i, j] > 0.5:  # Threshold for bond presence (adjust based on your generator)
                bond_type = Chem.BondType.SINGLE  # Assume single bonds for now
                # If your generator supports bond orders, modify this logic accordingly
                print(f"Adding bond: {i}-{j} Type: {bond_type}")
                mol.AddBond(atom_map[i], atom_map[j], bond_type)

    # Step 3: Check for isolated hydrogens
    isolated_hydrogens = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1 and len(atom.GetBonds()) == 0:
            isolated_hydrogens = True
            print(f"Warning: Isolated hydrogen detected at atom index {atom.GetIdx()}")
            break

    if isolated_hydrogens:
        raise ValueError("Generated molecule contains isolated hydrogens, which is chemically invalid")

    # Step 4: Sanitize the molecule
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        raise ValueError(f"Molecule sanitization failed: {e}")

    # Do NOT call Chem.RemoveHs here; let the caller handle it if needed
    return mol"""