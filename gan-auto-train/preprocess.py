import numpy as np
import tensorflow as tf

def preprocess_node_features(node_features_list, valid_atom_types, num_nodes):
    num_features = len(valid_atom_types)
    one_hot_features = np.zeros((num_nodes, num_features))

    atom_type_to_index = {atom: idx for idx, atom in enumerate(valid_atom_types)}

    for i, atom in enumerate(node_features_list):
        if atom[0] in atom_type_to_index:
            one_hot_features[i, atom_type_to_index[atom[0]]] = 1

    return one_hot_features

def pad_adj_matrix(matrix, max_size):
    padded_matrix = np.zeros((max_size, max_size))
    n, m = matrix.shape
    padded_matrix[:n, :m] = matrix
    return padded_matrix

def create_dataset(adj_matrices, one_hot_node_features, batch_size=32):
    dataset = tf.data.Dataset.from_tensor_slices((one_hot_node_features, adj_matrices))
    dataset = dataset.shuffle(buffer_size=len(one_hot_node_features)).batch(batch_size)
    return dataset

def pad_one_hot_features(features, max_size, num_features):
    padded_features = np.zeros((max_size, num_features))
    n, _ = features.shape
    padded_features[:n, :] = features
    return padded_features