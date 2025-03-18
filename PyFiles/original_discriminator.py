import tensorflow as tf
from tensorflow.keras import layers, models

class GraphDiscriminator(tf.keras.Model):
    def __init__(self, num_nodes, node_features):
        super().__init__()
        self.input_dim = num_nodes * num_nodes
        self.condition_dim = num_nodes * node_features

        self.mlp = models.Sequential([
            layers.Dense(256, input_shape = (self.input_dim + self.condition_dim,), activation='leaky_relu'),
            layers.Dense(128, activation='leaky_relu'),
            layers.BatchNormalization(),
            layers.Dense(1, activation='sigmoid')
        ])

        self.optimizer = tf.keras.optimizers.Adam(learning_rate=0.0002, beta_1=0.5)

    def call(self, adj, node_features):
        batch_size = tf.shape(adj)[0]

        adj_flat = tf.reshape(adj, (batch_size, self.input_dim))
        node_features_flat = tf.reshape(node_features, (batch_size, self.condition_dim))

        combined = tf.concat([adj_flat, node_features_flat], axis = 1)

        validity = self.mlp(combined)

        return validity