"""import os
import sys
import tensorflow as tf
import pandas as pd
import re
import ast
import numpy as np
from rdkit import Chem
#from PyFiles import hierarchical_discriminator as hd
#from PyFiles import hierarchical_discriminator_nd as hd_nd
#from PyFiles import rewards as rds
#from PyFiles import preprocess as pp
#from PyFiles import test_graph_gen as gen

sys.path.insert(0, os.path.abspath('./PyFiles_docker_test'))

from matrix_to_mol import adjacency_matrix_to_mol as adj
import preprocess as pp


def normalize_graph_features(x):
    x = re.sub(r'(?<=[0-9])\s+(?=[0-9])', ',', x)
    x = re.sub(r'\]\s+\[', '],[', x)
    return x

def smiles_to_graph(smiles):
    mol = Chem.MolFromSmiles(smiles)
    num_nodes = mol.GetNumAtoms()

    node_features = []
    for atom in mol.GetAtoms():
        node_features.append([atom.GetAtomicNum()])

    adj_matrix = Chem.GetAdjacencyMatrix(mol)

    return adj_matrix, node_features, num_nodes

def get_dataset():
    data = pd.read_csv("data/qm9.csv")
    #data = pd.read_csv("https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/qm9.csv")
    #data[['adj_matrix', 'node_features', 'node_count']] = data.apply(lambda row: pd.Series(smiles_to_graph(row['smiles'])), axis=1)
    #print(data[['adj_matrix', 'node_features', 'node_count']].tail())

    data['adj_matrix'] = data["adj_matrix"].apply(normalize_graph_features)
    data['node_features'] = data["node_features"].apply(normalize_graph_features)

    data['adj_matrix'] = data['adj_matrix'].apply(lambda x: np.array(ast.literal_eval(x)))
    data['node_features'] = data['node_features'].apply(lambda x: ast.literal_eval(x))

    valid_atom_types = [1, 6, 7, 8, 9]
    num_features = len(valid_atom_types)

    max_size = max(matrix.shape[0] for matrix in data['adj_matrix'])

    one_hot_node_features = [
        pp.preprocess_node_features(nf, valid_atom_types, num_nodes)
        for nf, num_nodes in zip(data['node_features'], data['node_count'])
    ]
    one_hot_node_features = [
        pp.pad_one_hot_features(features, max_size, num_features)
        for features in one_hot_node_features
    ]

    node_features = tf.argmax(one_hot_node_features, axis = -1)
    node_features = np.array(node_features)
    node_counts = data['node_count']
    sliced_features = []
    train_smiles = []

    for node_feature, node_count in zip(node_features, node_counts):
        sliced_features.append(node_feature[:node_count])

    feature_map = {
        0 : valid_atom_types[0],
        1 : valid_atom_types[1],
        2 : valid_atom_types[2],
        3 : valid_atom_types[3],
        4 : valid_atom_types[4]
    }
    mapped_features = [0] * len(sliced_features)

    for i in range(len(sliced_features)):
        mapped_features[i] = [feature_map[value] for value in sliced_features[i]]

    adj_matrices = data['adj_matrix'].values

    for i in range(len(adj_matrices)):
        train_smiles.append(Chem.MolToSmiles(adj(adj_matrices[i], mapped_features[i])))

    data['adj_matrix'] = data['adj_matrix'].apply(lambda x: pp.pad_adj_matrix(x, max_size))

    adj_matrices = np.stack(data['adj_matrix'].values)
    one_hot_node_features = np.stack(one_hot_node_features)

    dataset = pp.create_dataset(adj_matrices, one_hot_node_features, batch_size = 32)

    return dataset, train_smiles, max_size, num_features, adj_matrices, one_hot_node_features"""



import os
import sys
import tensorflow as tf
import pandas as pd
import re
import ast
import numpy as np
from rdkit import Chem
import preprocess as pp
import traceback  

sys.path.insert(0, os.path.abspath('./PyFiles_docker_test'))

from matrix_to_mol import adjacency_matrix_to_mol as adj

def normalize_graph_features(x):
    x = re.sub(r'(?<=[0-9])\s+(?=[0-9])', ',', x)
    x = re.sub(r'\]\s+\[', '],[', x)
    return x

def smiles_to_graph(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)  # Ensure hydrogens are present
    num_nodes = mol.GetNumAtoms()

    node_features = []
    for atom in mol.GetAtoms():
        node_features.append([atom.GetAtomicNum()])

    # Get the adjacency matrix
    adj_matrix = Chem.GetAdjacencyMatrix(mol)

    # Create a bond order matrix
    bond_order_matrix = np.zeros((num_nodes, num_nodes), dtype=np.float32)
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        bond_type = bond.GetBondType()
        if bond_type == Chem.BondType.SINGLE:
            bond_order = 1
        elif bond_type == Chem.BondType.DOUBLE:
            bond_order = 2
        elif bond_type == Chem.BondType.TRIPLE:
            bond_order = 3
        else:
            bond_order = 1  # Default to single bond for other types (e.g., aromatic)
        bond_order_matrix[i, j] = bond_order
        bond_order_matrix[j, i] = bond_order

    return adj_matrix, node_features, num_nodes, bond_order_matrix

def get_dataset(sample_size=10):
    data = pd.read_csv("https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/qm9.csv")
    
    """data = data.sample(n=sample_size, random_state=23).reset_index(drop=True)
    print(f"Subsampled dataset to {sample_size} molecules")"""

    data[['adj_matrix', 'node_features', 'node_count', 'bond_order_matrix']] = data['smiles'].apply(
        lambda s: pd.Series(smiles_to_graph(s))
    )

    valid_atom_types = [1, 6, 7, 8, 9] 
    num_features = len(valid_atom_types)

    max_size = max(matrix.shape[0] for matrix in data['adj_matrix'])

    one_hot_node_features = [
        pp.preprocess_node_features(nf, valid_atom_types, num_nodes)
        for nf, num_nodes in zip(data['node_features'], data['node_count'])
    ]
    one_hot_node_features = [
        pp.pad_one_hot_features(features, max_size, num_features)
        for features in one_hot_node_features
    ]

    node_features = tf.argmax(one_hot_node_features, axis=-1)
    node_features = np.array(node_features)
    node_counts = data['node_count']
    sliced_features = []
    train_smiles = []

    for node_feature, node_count in zip(node_features, node_counts):
        sliced_features.append(node_feature[:node_count])

    feature_map = {
        0: valid_atom_types[0],
        1: valid_atom_types[1],
        2: valid_atom_types[2],
        3: valid_atom_types[3],
        4: valid_atom_types[4]
    }
    mapped_features = [0] * len(sliced_features)

    for i in range(len(sliced_features)):
        mapped_features[i] = [feature_map[value] for value in sliced_features[i]]

    adj_matrices = data['adj_matrix'].values
    bond_order_matrices = data['bond_order_matrix'].values
    original_smiles = data['smiles'].values

    for i in range(len(adj_matrices)):
        #print(f"\nProcessing molecule {i}, Original SMILES: {original_smiles[i]}")
        try:
            node_count = node_counts[i]
            adj_matrix = adj_matrices[i][:node_count, :node_count]
            bond_order_matrix = bond_order_matrices[i][:node_count, :node_count]
            node_feature = mapped_features[i][:node_count]
            mol = adj(adj_matrix, node_feature, bond_order_matrix)
            smiles = Chem.MolToSmiles(mol, canonical=True, allHsExplicit=False)
            #print(f"Converted SMILES: {smiles}")
            original_mol = Chem.MolFromSmiles(original_smiles[i])
            original_canonical = Chem.MolToSmiles(original_mol, canonical=True)
            converted_canonical = Chem.MolToSmiles(mol, canonical=True)
            """if original_canonical == converted_canonical:
                print("Validation: SMILES strings are equivalent")
            else:
                print("Validation: SMILES strings are NOT equivalent")"""
            train_smiles.append(smiles)
        except Exception as e:
            print(f"Failed to convert adj matrix to SMILES at index {i}:")
            print(traceback.format_exc())
            train_smiles.append(None)

    data['adj_matrix'] = data['adj_matrix'].apply(lambda x: pp.pad_adj_matrix(x, max_size))

    adj_matrices = np.stack(data['adj_matrix'].values)
    one_hot_node_features = np.stack(one_hot_node_features)

    dataset = pp.create_dataset(adj_matrices, one_hot_node_features, batch_size=32)

    return dataset, train_smiles, max_size, num_features, adj_matrices, one_hot_node_features, node_counts, mapped_features

"""dataset, train_smiles, max_size, num_features, adj_matrices, one_hot_node_features, node_counts, mapped_features = get_dataset(sample_size=10)

print("First few train_smiles:", train_smiles[:5])
print("Max size:", max_size)
print("Num features:", num_features)
print("Adjacency matrices shape:", adj_matrices.shape)
print("One-hot node features shape:", one_hot_node_features.shape)

# Validate SMILES strings
for smiles in train_smiles[:5]:
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        print(f"SMILES: {smiles}, Valid: {mol is not None}")"""