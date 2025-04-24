import tensorflow as tf
import numpy as np
from rdkit import Chem

"""
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

def adjacency_matrix_to_mol(adj, nodes, bond_order_matrix=None, chiral_tags=None, bond_directions=None):
    adj = np.array(adj)
    nodes = np.array(nodes)
    if bond_order_matrix is not None:
        bond_order_matrix = np.array(bond_order_matrix)

    if adj.shape[0] != adj.shape[1] or adj.shape[0] != len(nodes):
        raise ValueError("Adjacency matrix dimensions must match the number of nodes")

    mol = Chem.RWMol()
    atom_map = {}
    for i, atom_num in enumerate(nodes):
        atom = Chem.Atom(int(atom_num))
        if chiral_tags and i in chiral_tags:
            atom.SetChiralTag(chiral_tags[i])
        mol_idx = mol.AddAtom(atom)
        atom_map[i] = mol_idx

    print("Adjacency matrix:")
    print(adj)
    print("Nodes:", nodes)
    if bond_order_matrix is not None:
        print("Bond order matrix:")
        print(bond_order_matrix)

    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            if adj[i, j] > 0.5:
                bond_order = int(bond_order_matrix[i, j]) if bond_order_matrix is not None and bond_order_matrix[i, j] > 0 else 1
                bond_type = Chem.BondType.SINGLE if bond_order == 1 else Chem.BondType.DOUBLE if bond_order == 2 else Chem.BondType.TRIPLE
                print(f"Adding bond: {i}-{j} Type: {bond_type}")
                bond = mol.AddBond(atom_map[i], atom_map[j], bond_type)
                if bond_directions and (i, j) in bond_directions:
                    mol.GetBondBetweenAtoms(atom_map[i], atom_map[j]).SetBondDir(bond_directions[(i, j)])
                elif bond_directions and (j, i) in bond_directions:
                    mol.GetBondBetweenAtoms(atom_map[i], atom_map[j]).SetBondDir(bond_directions[(j, i)])

    # Fix isolated hydrogens
    isolated_hydrogens = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1 and len(atom.GetBonds()) == 0:
            isolated_hydrogens.append(atom.GetIdx())

    if isolated_hydrogens:
        print(f"Warning: Found {len(isolated_hydrogens)} isolated hydrogens. Attempting to connect them.")
        for h_idx in isolated_hydrogens:
            # Find the nearest non-hydrogen atom with available valence
            best_candidate = None
            best_score = float('inf')
            for i in range(len(nodes)):
                if i == h_idx or nodes[i] == 1:
                    continue
                atom = mol.GetAtomWithIdx(i)
                # Estimate available valence (simplified)
                current_bonds = sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds())
                max_valence = 4 if nodes[i] == 6 else 3 if nodes[i] == 7 else 2 if nodes[i] == 8 else 1
                if current_bonds >= max_valence:
                    continue
                # Use adjacency matrix value as a "distance" metric
                score = -adj[h_idx, i]  # Higher adj value means closer
                if score < best_score:
                    best_score = score
                    best_candidate = i

            if best_candidate is not None:
                print(f"Connecting isolated hydrogen {h_idx} to atom {best_candidate}")
                mol.AddBond(h_idx, best_candidate, Chem.BondType.SINGLE)
            else:
                print(f"Warning: Could not find a suitable atom to connect isolated hydrogen {h_idx}. Removing it.")
                mol.RemoveAtom(h_idx)

    print("Molecule before sanitization:")
    for atom in mol.GetAtoms():
        print(f"Atom {atom.GetIdx()}: {atom.GetSymbol()}, Bonds: {len(atom.GetBonds())}")
    for bond in mol.GetBonds():
        print(f"Bond {bond.GetIdx()}: {bond.GetBeginAtomIdx()}-{bond.GetEndAtomIdx()} Type: {bond.GetBondType()}")

    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        raise ValueError(f"Molecule sanitization failed: {e}")

    return mol