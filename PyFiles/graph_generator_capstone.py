import tensorflow as tf
from tensorflow.keras import layers, models
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from matrix_to_mol import adjacency_matrix_to_mol
from rewards import validity_reward, uniqueness_reward, novelty_reward, drug_like_reward
import warnings

warnings.filterwarnings('ignore')
tf.config.run_functions_eagerly(True)

class GraphGenerator(tf.keras.Model):
    def __init__(self, num_nodes, node_features, latent_dim):
        super().__init__()

        self.num_nodes = num_nodes
        self.node_features = node_features
        self.latent_dim = latent_dim

        self.mlp = models.Sequential([
            layers.Dense(128, input_shape=(latent_dim,), activation="relu"),
            layers.Dense(256, activation="relu"),
            layers.Dense((num_nodes * num_nodes) + (num_nodes * node_features), activation="softmax")
        ])

        self.optimizer = tf.keras.optimizers.Adam(learning_rate=0.0002, beta_1=0.5)

    def call(self, z):
        batch_size = tf.shape(z)[0]

        output = self.mlp(z)

        adj_flat = output[:, :self.num_nodes * self.num_nodes]
        adj_flat = tf.nn.sigmoid(adj_flat)
        node_feature_flat = output[:, self.num_nodes * self.num_nodes:]

        adj = tf.reshape(adj_flat, (batch_size, self.num_nodes, self.num_nodes))
        adj = (adj + tf.transpose(adj, perm=[0, 2, 1])) / 2
        adj = tf.nn.softmax(adj, axis=-1)
        adj = tf.cast(adj, tf.float32)

        valid_atom_types = tf.constant([1, 6, 7, 8, 9], dtype=tf.int64)

        node_features = tf.reshape(node_feature_flat, (batch_size, self.num_nodes, self.node_features))
        node_features = tf.nn.softmax(node_features, axis=-1)
        node_features = tf.clip_by_value(node_features, 0, 4)
        node_features = tf.cast(node_features, tf.int32)
        node_features = tf.gather(valid_atom_types, node_features)
        node_features = tf.cast(node_features, tf.float32)

        return adj, node_features

    def loss_function(self, real_output, fake_output):
        loss_func = tf.keras.losses.BinaryCrossentropy(from_logits=False)
        return loss_func(real_output, fake_output)

    def fit(self, dataset, discriminator, train_smiles, epochs=10):
        d_loss_list = []
        g_loss_list = []
        r_loss_list = []
        scl_loss_list = []
        gen_smiles_list = []
        for epoch in range(epochs):
            for batch_idx, (real_node_features, real_adj) in enumerate(dataset):
                z = tf.random.normal([tf.shape(real_node_features)[0], self.latent_dim])
                gen_adj, gen_node_features = self.call(z)

                with tf.GradientTape() as tape:
                    fake_output = discriminator.call(gen_adj, gen_node_features)
                    real_output = discriminator.call(real_adj, real_node_features)

                    # Label smoothing test
                    real_labels = tf.ones_like(real_output) * 0.9
                    fake_labels = tf.zeros_like(fake_output) + 0.1

                    # Label switching test
                    """if np.random.rand() < 0.05:
                        real_labels, fake_labels = fake_labels, real_labels"""

                    d_loss = self.loss_function(real_labels, real_output) + \
                             self.loss_function(fake_labels, fake_output)

                    d_loss_list.append(d_loss.numpy())

                    gradients = tape.gradient(d_loss, discriminator.trainable_variables)
                    discriminator.optimizer.apply_gradients(zip(gradients, discriminator.trainable_variables))

                z = tf.random.normal([tf.shape(real_node_features)[0], self.latent_dim])

                with tf.GradientTape() as tape:
                    gen_adj, gen_node_features = self(z)
                    gen_combined = [gen_adj, gen_node_features]
                    gen_smiles_list.append(Chem.MolToSmiles(adjacency_matrix_to_mol(gen_combined)))
                    curr_mol = adjacency_matrix_to_mol(gen_combined)
                    validity = tf.convert_to_tensor(validity_reward(curr_mol), dtype = tf.float32)
                    uniqueness = tf.convert_to_tensor(uniqueness_reward(gen_smiles_list, curr_mol), dtype = tf.float32) 
                    novelty = tf.convert_to_tensor(novelty_reward(curr_mol, train_smiles, gen_smiles_list), dtype = tf.float32) 
                    drug_like = tf.convert_to_tensor(drug_like_reward(curr_mol), dtype = tf.float32)

                    total_reward = tf.convert_to_tensor(validity + uniqueness + novelty + drug_like, dtype = tf.float32)
                    r_loss = total_reward  
                    real_r_loss = r_loss.numpy()
                    r_loss = tf.reshape(r_loss, [1])  

                    fake_output = discriminator(gen_adj, gen_node_features)

                    g_loss = self.loss_function(tf.ones_like(fake_output), fake_output)
                    scaled_loss = g_loss * tf.exp(-r_loss)
                    real_scaled_loss = g_loss.numpy() * np.exp(real_r_loss)
                    g_loss_list.append(g_loss.numpy())
                    scl_loss_list.append(real_scaled_loss)

                    gradients = tape.gradient(scaled_loss, self.trainable_variables)
                    self.optimizer.apply_gradients(zip(gradients, self.trainable_variables))

                    if epoch % 10 == 0:
                        print(f"Epoch {epoch+1}/{epochs}",
                            f"D Loss: {d_loss.numpy():.4f}",
                            f"G Loss: {g_loss.numpy():.4f}", 
                            f"R Loss: {real_r_loss:.4f}", 
                            f"Scaled Loss: {real_scaled_loss:.4f}")

        return d_loss_list, g_loss_list, r_loss_list, scl_loss_list