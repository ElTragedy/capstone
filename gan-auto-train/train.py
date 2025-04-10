import os
import tensorflow as tf
from rdkit import Chem
from PyFiles import hierarchical_discriminator as hd
from PyFiles import hierarchical_discriminator_nd as hd_nd
from PyFiles import rewards as rds
from PyFiles import test_graph_gen as gen
from PyFiles import checkpoint as cp
import get_dataset

DATA_DIR = "./data"
OUTPUT_DIR = "./output"
os.makedirs(OUTPUT_DIR, exist_ok = True)

EPOCHS = 1000
BATCH_SIZE = 32
SAVE_EVERY = 25
LATENT_DIM = 256

dataset, train_smiles = get_dataset.get_dataset()

generator = gen.test_GraphGenerator()
no_diff_gen = gen.test_GraphGenerator()
diff_discriminator = hd.GraphDiscriminatorDiff()
nodiff_discriminator = hd_nd.GraphDiscriminatorNoDiff()

gen_optimizer = tf.keras.optimizers.Adam(learning_rate = 0.0002, beta_1 = 0.5)
diff_disc_optimizer = tf.keras.optimizers.SGD(learning_rate = 0.0002, momentum = 0.5)
nodiff_disc_optimizer = tf.keras.optimizers.SGD(learning_rate = 0.0002, momentum = 0.5)

d_loss_diff, g_loss_diff, r_loss_diff, scl_loss_diff = generator.fit(dataset, diff_discriminator, train_smiles, SAVE_EVERY, OUTPUT_DIR, EPOCHS)
d_loss_nodiff, g_loss_nodiff, r_loss_nodiff, scl_loss_nodiff = no_diff_gen.fit(dataset, nodiff_discriminator, train_smiles, SAVE_EVERY, OUTPUT_DIR, EPOCHS)