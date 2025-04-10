import tensorflow as tf
import pandas as pd
import re
import ast
import numpy as np
from rdkit import Chem
from PyFiles import hierarchical_discriminator as hd
from PyFiles import hierarchical_discriminator_nd as hd_nd
from PyFiles import rewards as rds
from PyFiles.matrix_to_mol import adjacency_matrix_to_mol as adj
from PyFiles import preprocess as pp
from PyFiles import test_graph_gen as gen

def normalize_graph_features(x):
    x = re.sub(r'(?<=[0-9])\s+(?=[0-9])', ',', x)
    x = re.sub(r'\]\s+\[', '],[', x)
    return x

def get_dataset():
    data = pd.read_csv("gan-auto-train/data/qm9.csv")

    data['adj_matrix'] = data["adj_matrix"].apply(normalize_graph_features)
    data['node_features'] = data["node_features"].apply(normalize_graph_features)

    data['adj_matrix'] = data['adj_matrix'].apply(lambda x: np.array(ast.literal_eval(x)))
    data['node_features'] = data['node_features'].apply(lambda x: ast.literal_eval(x))

    valid_atom_types = [1, 6, 7, 8, 9]
    num_features = len(valid_atom_types)

    max_size = max(matrix.shape[0] for matrix in qm9['adj_matrix'])

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

    one_hot_node_features = np.stack(one_hot_node_features)

    dataset = pp.create_dataset(adj_matrices, one_hot_node_features, batch_size = 32)

    return dataset