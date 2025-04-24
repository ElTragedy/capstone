import os
import sys
import tensorflow as tf
import pandas as pd
import re
import ast
import numpy as np
from rdkit import Chem
import json

sys.path.insert(0, os.path.abspath('./PyFiles_docker_test'))

from matrix_to_mol import adjacency_matrix_to_mol as adj
import preprocess as pp


def normalize_graph_features(x):
    x = re.sub(r'(?<=[0-9])\s+(?=[0-9])', ',', x)
    x = re.sub(r'\]\s+\[', '],[', x)
    return x

def get_dataset():
    data = pd.read_csv("data/qm9.csv")

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

    """for i in range(len(sliced_features)):
        mapped_features[i] = [feature_map[value] for value in sliced_features[i]]

    adj_matrices = data['adj_matrix'].values

    for i in range(len(adj_matrices)):
        train_smiles.append(Chem.MolToSmiles(adj(adj_matrices[i], mapped_features[i])))"""

    for i in range(len(sliced_features)):
        mapped_features[i] = [feature_map[value] for value in sliced_features[i]]

    adj_matrices = data['adj_matrix'].values

    for i in range(len(adj_matrices)):
        try:
            node_count = len(mapped_features[i])
            trimmed_adj = adj_matrices[i][:node_count, :node_count]
            mol = adj(trimmed_adj, mapped_features[i])
            if mol is not None:
                smiles = Chem.MolToSmiles(mol)
                train_smiles.append(smiles)
        except Exception as e:
            print(f"[Warning] Failed to convert molecule at index {i}: {e}")

        data['adj_matrix'] = data['adj_matrix'].apply(lambda x: pp.pad_adj_matrix(x, max_size))

    adj_matrices = np.stack(data['adj_matrix'].values)
    one_hot_node_features = np.stack(one_hot_node_features)

    dataset = pp.create_dataset(adj_matrices, one_hot_node_features, batch_size = 32)

    return dataset, train_smiles, max_size, num_features

if __name__ == "__main__":
    dataset, train_smiles, max_size, num_features = get_dataset()

    smiles_df = pd.DataFrame({'smiles': train_smiles})
    smiles_df.to_csv('data/train_smiles.csv', index=False)
    dataset = pd.DataFrame(dataset)
    dataset.to_csv('data/dataset.csv', index=False)


    with open('data/graph_config.json', 'w') as f:
        json.dump({'max_size': max_size, 'num_features': num_features}, f)
    print("[Saved] graph_config.json")

