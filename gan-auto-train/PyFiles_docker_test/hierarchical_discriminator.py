import tensorflow as tf
from keras import layers, models

#base graph discriminator without pooling
class GraphDiscriminator(tf.keras.Model):
    def __init__(self, num_nodes, node_features):
        super(GraphDiscriminator, self).__init__()

        self.num_nodes = num_nodes
        self.node_features = node_features

        #dense layers to process node features
        self.conv1 = layers.Dense(64, activation = 'leaky_relu')
        self.conv2 = layers.Dense(128, activation = 'leaky_relu')

        #global pooling to combine features across nodes
        self.pool = layers.GlobalAveragePooling1D()

        #fully connected output layer (sigmoid for binary classification)
        self.fc = layers.Dense(1, activation = 'sigmoid')

        #optmizer for training the discriminator
        self.optimizer = tf.keras.optimizers.SGD(learning_rate=0.000099, momentum=0.5)
        
    #forward pass for discriminator; adg - adjacency matrix, node_features: node feature matrix
    def call(self, adj, node_features):
        x = self.conv1(node_features)
        x = self.conv2(x)
        x = self.pool(x)
        
        output = self.fc(x)
        return output
    """def call(self, z):
        batch_size = tf.shape(z)[0]

        output = self.mlp(z)

        adj_flat = output[:, :self.num_nodes * self.num_nodes]
        adj_flat = tf.nn.sigmoid(adj_flat)
        node_feature_flat = output[:, self.num_nodes * self.num_nodes:]

        adj = tf.reshape(adj_flat, (batch_size, self.num_nodes, self.num_nodes))
        adj = (adj + tf.transpose(adj, perm=[0, 2, 1])) / 2
        adj = tf.nn.softmax(adj, axis=-1)

        valid_atom_types = tf.constant([1, 6, 7, 8, 9], dtype=tf.int64)

        # Compute node features for the discriminator (one-hot encoded)
        node_features = tf.reshape(node_feature_flat, (batch_size, self.num_nodes, self.node_features))
        node_features = tf.nn.softmax(node_features, axis=-1)  # Probabilities over atom types
        # Convert to indices for atomic numbers
        node_indices = tf.argmax(node_features, axis=-1)  # Shape: (batch_size, num_nodes)
        atomic_numbers = tf.gather(valid_atom_types, node_indices)  # Shape: (batch_size, num_nodes)
        # Convert node_features to one-hot encoded format for the discriminator
        node_features = tf.one_hot(node_indices, depth=self.node_features, dtype=tf.float32)  # Shape: (batch_size, num_nodes, node_features)

        # Post-process adjacency matrix to ensure hydrogens are connected
        hydrogen_mask = tf.cast(tf.equal(atomic_numbers, 1), tf.float32)
        for b in range(batch_size):
            for i in range(self.num_nodes):
                if hydrogen_mask[b, i] == 1:
                    connections = tf.reduce_sum(adj[b, i, :])
                    if connections < 0.5:
                        adj_row = adj[b, i, :]
                        adj_row = tf.where(tf.range(self.num_nodes) == i, tf.zeros_like(adj_row), adj_row)
                        max_idx = tf.argmax(adj_row, axis=-1)
                        adj = tf.tensor_scatter_nd_update(
                            adj,
                            [[b, i, max_idx], [b, max_idx, i]],
                            [0.9, 0.9]
                        )

        adj = tf.cast(adj, tf.float32)
        atomic_numbers = tf.cast(atomic_numbers, tf.float32)

        return adj, node_features, atomic_numbers
"""