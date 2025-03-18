import tensorflow as tf
from tensorflow.keras import Model, layers
import numpy as np
from spektral.layers import GCNConv, DiffPool
from spektral.data import BatchLoader, Graph
from spektral.datasets import qm9
from spektral.transforms import GCNFilter

# -----------------------------
# Gumbel-Softmax Utility
# -----------------------------
def sample_gumbel(shape, eps=1e-20):
    U = tf.random.uniform(shape)
    return -tf.math.log(-tf.math.log(U + eps) + eps)

def gumbel_softmax_sample(logits, temperature):
    y = logits + sample_gumbel(tf.shape(logits))
    return tf.nn.softmax(y / temperature)

def gumbel_softmax(logits, temperature=0.5):
    return gumbel_softmax_sample(logits, temperature)


# -----------------------------
# Graph Generator
# -----------------------------
class GraphGenerator(Model):
    def __init__(self, num_nodes, num_node_types, num_edge_types, latent_dim=128):
        super().__init__()
        self.num_nodes = num_nodes
        self.num_node_types = num_node_types
        self.num_edge_types = num_edge_types
        self.latent_dim = latent_dim

        self.dense_proj = layers.Dense(256, activation='relu')
        self.node_logits = layers.Dense(num_nodes * num_node_types)
        self.adj_logits = layers.Dense(num_nodes * num_nodes * num_edge_types)

    def call(self, z):
        h = self.dense_proj(z)
        node_logits = self.node_logits(h)
        adj_logits = self.adj_logits(h)

        node_logits = tf.reshape(node_logits, (-1, self.num_nodes, self.num_node_types))
        adj_logits = tf.reshape(adj_logits, (-1, self.num_nodes, self.num_nodes, self.num_edge_types))

        node_sample = gumbel_softmax(node_logits)
        adj_sample = gumbel_softmax(adj_logits)

        return node_sample, adj_sample


# -----------------------------
# DiffPool Discriminator
# -----------------------------
class DiffPoolDiscriminator(Model):
    def __init__(self, num_node_features, hidden_channels=64, num_clusters=5):
        super().__init__()
        self.pre_gcn = GCNConv(hidden_channels, activation='relu')
        self.pool_gcn = GCNConv(num_clusters, activation='softmax')
        self.embed_gcn = GCNConv(hidden_channels, activation='relu')
        self.diffpool = DiffPool(hidden_channels, num_clusters, activation='relu')
        self.final_gcn = GCNConv(hidden_channels, activation='relu')
        self.global_pool = layers.GlobalAveragePooling1D()
        self.out_layer = layers.Dense(1)  # binary classification

    def call(self, inputs):
        x, a, i = inputs
        z = self.pre_gcn([x, a])
        s = self.pool_gcn([x, a])
        z, a, i = self.diffpool([z, a, i], s)
        z = self.final_gcn([z, a])
        pooled = self.global_pool(z)
        return self.out_layer(pooled)


# -----------------------------
# Training Step
# -----------------------------
def train_step(generator, discriminator, opt_g, opt_d, batch_size, latent_dim):
    z = tf.random.normal((batch_size, latent_dim))
    fake_nodes, fake_adjs = generator(z)
    fake_adj = tf.reduce_max(fake_adjs, axis=-1)
    fake_i = tf.repeat(tf.range(batch_size), repeats=[fake_nodes.shape[1]] * batch_size)

    real_data = next(iter(load_qm9_dataset(batch_size=batch_size)))
    real_x, real_a, real_i = real_data

    with tf.GradientTape() as tape_d:
        real_logits = discriminator((real_x, real_a, real_i))
        fake_logits = discriminator((fake_nodes, fake_adj, fake_i))
        d_loss = tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(
            labels=tf.ones_like(real_logits), logits=real_logits)) +                  tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(
            labels=tf.zeros_like(fake_logits), logits=fake_logits))
    grads_d = tape_d.gradient(d_loss, discriminator.trainable_variables)
    opt_d.apply_gradients(zip(grads_d, discriminator.trainable_variables))

    with tf.GradientTape() as tape_g:
        fake_nodes, fake_adjs = generator(z)
        fake_adj = tf.reduce_max(fake_adjs, axis=-1)
        fake_logits = discriminator((fake_nodes, fake_adj, fake_i))
        g_loss = tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(
            labels=tf.ones_like(fake_logits), logits=fake_logits))
    grads_g = tape_g.gradient(g_loss, generator.trainable_variables)
    opt_g.apply_gradients(zip(grads_g, generator.trainable_variables))

    return d_loss.numpy(), g_loss.numpy()


# -----------------------------
# Load QM9 Dataset
# -----------------------------
def load_qm9_dataset(sample_size=1000, batch_size=16):
    dataset = qm9.QM9(amount=sample_size, transforms=GCNFilter())
    loader = BatchLoader(dataset, batch_size=batch_size, epochs=1)
    return loader


# -----------------------------
# Run Training Loop
# -----------------------------
def run_training(epochs=10, batch_size=16, latent_dim=128):
    num_nodes = 9
    num_node_types = 5
    num_edge_types = 3
    num_clusters = 5

    generator = GraphGenerator(num_nodes, num_node_types, num_edge_types, latent_dim)
    discriminator = DiffPoolDiscriminator(num_node_types, hidden_channels=64, num_clusters=num_clusters)

    opt_g = tf.keras.optimizers.Adam(1e-4)
    opt_d = tf.keras.optimizers.Adam(1e-4)

    for epoch in range(epochs):
        d_loss, g_loss = train_step(generator, discriminator, opt_g, opt_d, batch_size, latent_dim)
        print(f"Epoch {epoch + 1} - D Loss: {d_loss:.4f} | G Loss: {g_loss:.4f}")


if __name__ == "__main__":
    run_training()
