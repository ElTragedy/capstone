import os
import sys
import tensorflow as tf
from rdkit import Chem

sys.path.insert(0, os.path.abspath('./PyFiles_docker_test'))
#print("PYTHONPATH =", sys.path)
#print("Checking for hierarchical_discriminator.py:", os.path.exists('./PyFiles/hierarchical_discriminator.py'))

import hierarchical_discriminator as hd
import hierarchical_discriminator_nd as hd_nd
import rewards as rds
import test_graph_gen as gen
import checkpoint as cp
import get_dataset2

DATA_DIR = "./data"
OUTPUT_DIR = "./output"
os.makedirs(OUTPUT_DIR, exist_ok = True)

EPOCHS = 1000
BATCH_SIZE = 32
SAVE_EVERY = 25
LATENT_DIM = 256

dataset, train_smiles, NUM_NODES, NUM_FEATURES, adj_matrices, one_hot_node_features, node_counts, mapped_features = get_dataset2.get_dataset()

generator = gen.test_GraphGenerator(NUM_NODES, NUM_FEATURES, LATENT_DIM)
no_diff_gen = gen.test_GraphGenerator(NUM_NODES, NUM_FEATURES, LATENT_DIM)
diff_discriminator = hd.GraphDiscriminator(NUM_NODES, NUM_FEATURES)
nodiff_discriminator = hd_nd.GraphDiscriminatorNoDiff(NUM_NODES, NUM_FEATURES)

gen_optimizer = tf.keras.optimizers.Adam(learning_rate = 0.0002, beta_1 = 0.5)
diff_disc_optimizer = tf.keras.optimizers.SGD(learning_rate = 0.0002, momentum = 0.5)
nodiff_disc_optimizer = tf.keras.optimizers.SGD(learning_rate = 0.0002, momentum = 0.5)

d_loss_diff, g_loss_diff, r_loss_diff, scl_loss_diff = generator.fit(dataset, diff_discriminator, train_smiles, SAVE_EVERY, OUTPUT_DIR, EPOCHS)
d_loss_nodiff, g_loss_nodiff, r_loss_nodiff, scl_loss_nodiff = no_diff_gen.fit(dataset, nodiff_discriminator, train_smiles, SAVE_EVERY, OUTPUT_DIR, EPOCHS)