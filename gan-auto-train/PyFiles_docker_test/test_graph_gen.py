import tensorflow as tf
from tensorflow.keras import layers, models
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from matrix_to_mol import adjacency_matrix_to_mol
import checkpoint as cp
import rewards
import warnings

warnings.filterwarnings('ignore')
tf.config.run_functions_eagerly(True)

class test_GraphGenerator(tf.keras.Model):
    def __init__(self, num_nodes, node_features, latent_dim):
        super().__init__()

        self.num_nodes = num_nodes
        self.node_features = node_features
        self.latent_dim = latent_dim

        self.mlp = models.Sequential([
            layers.Dense(128, activation="leaky_relu"),
            layers.Dropout(0.3),
            layers.BatchNormalization(),
            layers.Dense(256, activation="leaky_relu"),
            layers.Dropout(0.3),
            layers.BatchNormalization(),
            layers.Dense(512, activation="leaky_relu"),
            layers.Dropout(0.3),
            layers.BatchNormalization(),
            layers.Dense((num_nodes * num_nodes) + (num_nodes * node_features), activation="sigmoid")])

        self.optimizer = tf.keras.optimizers.Adam(learning_rate=0.0002, beta_1=0.5)

    """def call(self, z):
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

        return adj, node_features"""
    
    def call(self, z):
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

    def loss_function(self, real_output, fake_output):
        loss_func = tf.keras.losses.BinaryCrossentropy(from_logits=False)
        return loss_func(real_output, fake_output)

    def log_prob(self, z):
        logits = self.mlp(z)
        return tf.nn.log_softmax(logits)

    def fit(self, dataset, discriminator, train_smiles, SAVE_EVERY, OUTPUT_DIR, epochs=10):
        d_loss_list = []
        g_loss_list = []
        r_loss_list = []
        scl_loss_list = []
        gen_smiles_list = []
        for epoch in range(epochs):
            for batch_idx, (real_node_features, real_adj) in enumerate(dataset):
                z = tf.random.normal([tf.shape(real_node_features)[0], self.latent_dim])
                gen_adj, gen_node_features, gen_atomic_numbers = self.call(z)

                with tf.GradientTape() as tape:
                    fake_output = discriminator.call(gen_adj, gen_node_features)
                    real_output = discriminator.call(real_adj, real_node_features)

                    real_labels = tf.ones_like(real_output) * 0.9
                    fake_labels = tf.zeros_like(fake_output) + 0.1

                    d_loss = self.loss_function(real_labels, real_output) + \
                            self.loss_function(fake_labels, fake_output)

                    d_loss_list.append(d_loss.numpy())

                    gradients = tape.gradient(d_loss, discriminator.trainable_variables)
                    discriminator.optimizer.apply_gradients(zip(gradients, discriminator.trainable_variables))

                for _ in range(2):
                    z = tf.random.normal([tf.shape(real_node_features)[0], self.latent_dim])

                    with tf.GradientTape() as tape:
                        gen_adj, gen_node_features, gen_atomic_numbers = self(z)
                        batch_rewards = []
                        batch_size = gen_adj.shape[0]
                        gen_adj_np = gen_adj.numpy()
                        gen_atomic_numbers_np = gen_atomic_numbers.numpy()

                        for i in range(batch_size):
                            adj_i = gen_adj_np[i]
                            nodes_i = gen_atomic_numbers_np[i]

                            try:
                                curr_mol = adjacency_matrix_to_mol(adj_i, nodes_i)
                                if curr_mol is None or curr_mol.GetNumAtoms() == 0:
                                    print("Invalid molecule generated.")
                                    continue
                                curr_mol = Chem.RemoveHs(curr_mol)
                                gen_smiles = Chem.MolToSmiles(curr_mol, canonical=True, allHsExplicit=False)
                                gen_smiles_list.append(gen_smiles)

                                total_reward = rewards.calculate_total_reward(curr_mol, gen_smiles_list, train_smiles)
                                batch_rewards.append(total_reward)
                            except Exception as e:
                                print(f"Failed to process molecule {i} in batch: {e}")
                                continue

                        if not batch_rewards:
                            continue

                        batch_reward = np.mean(batch_rewards)
                        alpha = 0.99
                        if hasattr(self, "reward_baseline"):
                            self.reward_baseline = alpha * self.reward_baseline + (1 - alpha) * batch_reward
                        else:
                            self.reward_baseline = batch_reward

                        advantage = batch_reward - self.reward_baseline
                        r_loss = tf.convert_to_tensor(advantage, dtype=tf.float32)
                        real_r_loss = r_loss.numpy()
                        r_loss = tf.reshape(r_loss, [1])

                        fake_output = discriminator(gen_adj, gen_node_features)

                        g_loss = self.loss_function(tf.ones_like(fake_output), fake_output)
                        log_probs = self.log_prob(z)
                        lambda_g = 3.0
                        lambda_r = 7.0
                        scaled_loss = lambda_g * g_loss - lambda_r * tf.reduce_mean(log_probs * r_loss)
                        g_loss_list.append(g_loss.numpy())
                        scl_loss_list.append(scaled_loss.numpy())

                        gradients = tape.gradient(scaled_loss, self.trainable_variables)
                        self.optimizer.apply_gradients(zip(gradients, self.trainable_variables))

                if (epoch + 1) % SAVE_EVERY == 0:
                    cp.save_checkpoint(self, discriminator, epoch + 1, OUTPUT_DIR)
                    print(f"Epoch {epoch+1}/{epochs}",
                        f"D Loss: {d_loss.numpy():.4f}",
                        f"G Loss: {g_loss.numpy():.4f}",
                        f"R Loss: {real_r_loss:.4f}",
                        f"Scaled Loss: {scaled_loss.numpy():.4f}")

        return d_loss_list, g_loss_list, r_loss_list, scl_loss_list