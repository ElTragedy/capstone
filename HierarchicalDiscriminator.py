import tensorflow as tf
from tensorflow.keras import layers

class HierarchicalDiscriminator(tf.keras.Model):
    def __init__(self, num_nodes, node_features, pool_ratio=0.5):
        super().__init__()
        self.num_nodes = num_nodes
        self.node_features = node_features
        self.pool_ratio = pool_ratio
        
        self.gcn1 = layers.Dense(128, activation="relu")  # Simulating graph conv layer
        self.gcn2 = layers.Dense(64, activation="relu")
        
        self.pooling = layers.Dense(int(num_nodes * pool_ratio), activation="softmax")  # DiffPool
        
        self.mlp = tf.keras.Sequential([
            layers.Dense(64, activation="relu"),
            layers.Dense(1, activation="sigmoid")  # Output validity score
        ])

        self.optimizer = tf.keras.optimizers.Adam(learning_rate=0.0002, beta_1=0.5)
    
    def call(self, adj, node_features):
        x = self.gcn1(node_features)
        x = self.gcn2(x)
        
        pool_scores = self.pooling(x)  # Compute node grouping

        adj = tf.cast(adj, dtype=tf.float32)
        pool_scores = tf.cast(pool_scores, dtype=tf.float32)


        pooled_x = tf.linalg.matmul(tf.transpose(pool_scores, perm = [0, 2, 1]), x)  # Pool node features
        pooled_adj = tf.linalg.matmul(tf.transpose(pool_scores, perm = [0, 2, 1]), tf.linalg.matmul(adj, pool_scores))  # Pool adjacency
        
        combined = tf.concat([tf.reshape(pooled_x, [x.shape[0], -1]), tf.reshape(pooled_adj, [x.shape[0], -1])], axis=-1)
        validity = self.mlp(combined)
        
        return validity