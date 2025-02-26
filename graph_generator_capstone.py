"""graph_generator_capstone

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1EJMhGqLKkjUYUDSNVYljXQ7cWS720it-
"""

def mount_drive(dir):
  import os
  colab = 1
  if colab == 1:
    from google.colab import drive
    drive.mount('/content/drive', force_remount = True)
    current_folder = dir
    dest_folder = '/content/drive/My Drive/' + current_folder
    os.chdir(dest_folder)
    print('\n Current path: ' + os.getcwd())

mount_drive('')

!pip install rdkit selfies

import numpy as np
import pandas as pd
from rdkit import Chem
import tensorflow as tf
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt
import tensorflow.keras as keras
from tensorflow.keras import layers, models
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

qm9 = pd.read_csv("qm9.csv")
qm9.head()

class GraphGenerator_test(tf.keras.Model):
  def __init__(self, num_nodes, node_features, latent_dim):
     super().__init__()

     self.num_nodes = num_nodes
     self.node_features = node_features
     self.latent_dim = latent_dim

     self.mlp = models.Sequential([
         layers.Dense(128, input_shape = (latent_dim,), activation = "leaky_relu"),
         layers.Dense(256, activation = "leaky_relu"),
         layers.Dense((num_nodes * num_nodes) + (num_nodes * node_features), activation = "softmax")
     ])

     self.optimizer = tf.keras.optimizers.Adam(learning_rate = 0.001)

  def call(self, z):
    batch_size = tf.shape(z)[0]

    output = self.mlp(z)

    adj_flat = output[:, :self.num_nodes * self.num_nodes]
    adj_flat = tf.nn.sigmoid(adj_flat)
    node_feature_flat = output[:, self.num_nodes * self.num_nodes:]

    adj = tf.reshape(adj_flat, (batch_size, self.num_nodes, self.num_nodes))
    adj = (adj + tf.transpose(adj, perm = [0, 2, 1])) / 2
    adj = tf.nn.softmax(adj, axis=-1)
    adj = tf.cast(adj, tf.float32)

    valid_atom_types = tf.constant([1, 6, 7, 8, 9], dtype = tf.int64)

    node_features = tf.reshape(node_feature_flat, (batch_size, self.num_nodes, self.node_features))
    node_features = tf.nn.softmax(node_features, axis = -1)
    node_features = tf.clip_by_value(node_features, 0, 4)
    node_features = tf.cast(node_features, tf.int32)
    node_features = tf.gather(valid_atom_types, node_features)
    node_features = tf.cast(node_features, tf.float32)

    return adj, node_features

  def loss_function(self, real_output, fake_output):
    loss_func = tf.keras.losses.BinaryCrossentropy(from_logits = False)
    return loss_func(real_output, fake_output)

  def fit(self, data, discriminator, epochs = 10):
    for epoch in range(epochs):
      total_loss = 0
      num_batches = 0
      for z in data:
        with tf.GradientTape() as tape:
          gen_adj, gen_z = self(z)

          fake_output = discriminator(gen_adj, gen_z)
          loss = self.loss_function(tf.ones_like(fake_output), fake_output)

        gradients = tape.gradient(loss, self.trainable_variables)
        self.optimizer.apply_gradients(zip(gradients, self.trainable_variables))

        total_loss += loss
        num_batches += 1
      average_loss = total_loss / num_batches
      print(f"Epoch {epoch+1}/{epochs}, Loss: {average_loss.numpy():.4f}")

class GraphDiscriminator_test(tf.keras.Model):
    def __init__(self, num_nodes, node_features):
        super().__init__()
        self.input_dim = num_nodes * num_nodes
        self.condition_dim = num_nodes * node_features

        self.mlp = models.Sequential([
            layers.Dense(256, input_shape = (self.input_dim + self.condition_dim,), activation = 'leaky_relu'),
            layers.Dense(128, activation = 'leaky_relu'),
            layers.BatchNormalization(),
            layers.Dense(1, activation = 'sigmoid')
        ])

    def call(self, adj, node_features):
        batch_size = tf.shape(adj)[0]

        adj_flat = tf.reshape(adj, (batch_size, self.input_dim))
        node_features_flat = tf.reshape(node_features, (batch_size, self.condition_dim))

        combined = tf.concat([adj_flat, node_features_flat], axis = 1)

        validity = self.mlp(combined)

        return validity

def adjacency_matrix_to_mol(inp_matrix):
  mol = Chem.RWMol()
  atom_types = np.unique(inp_matrix[1])
  matrix = inp_matrix[0]
  atom_map = {i: mol.AddAtom(Chem.Atom(int(atom_types[i]))) for i in range(len(atom_types))}

  if isinstance(matrix, tf.Tensor):
        matrix = matrix.numpy()
  if isinstance(atom_types, tf.Tensor):
      atom_types = atom_types.numpy()

  matrix = np.array(matrix)
  atom_types = np.array(atom_types)

  if len(matrix) != len(atom_types):
    raise ValueError("NUMBER OF ATOM TYPES DOES NOT MATCH MATRIX DIMENSIONS")

  # Step 1: Add Atoms
  for i, atom_num in enumerate(atom_types):
    atom = Chem.Atom(int(atom_num))
    mol_idx = mol.AddAtom(atom)
    atom_map[i] = mol_idx

  # Step 2: Add Bonds
  for i in range(len(matrix)):
    for j in range(i + 1, len(matrix)):
        if j not in atom_map:
            continue

        value = int(np.argmax(matrix[i, j]))
        value = min(max(value, 0), 4)

        if value != 0:
            bond_type = {
                1: Chem.BondType.SINGLE,
                2: Chem.BondType.DOUBLE,
                3: Chem.BondType.TRIPLE,
                4: Chem.BondType.AROMATIC
            }.get(value, None)

            if bond_type is None:
                raise ValueError(f"INVALID BOND TYPE DETECTED: {value}")

            print(f"Adding bond: {i}-{j} Type: {bond_type}")
            mol.AddBond(atom_map[i], atom_map[j], bond_type)

  return mol

condition_features = ['mu', 'homo', 'gap']
condition_data = qm9.loc[:4, condition_features].values.astype(np.float32)

scaler = StandardScaler()
condition_data = scaler.fit_transform(condition_data)

print("Condition Data Shape:", condition_data.shape)
print(condition_data)

gen = GraphGenerator_test(10,5,32)
z = tf.random.normal((5, 32))
gen_out = gen(z)

gen.summary()

gen_out

mol = adjacency_matrix_to_mol(gen_out)

mol
