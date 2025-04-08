import os
import tensorflow as tf
import pandas as pd
import hierarchical_discriminator as hd
import hierarchical_discriminator_nd as hd_nd
import rewards as rds
from matrix_to_mol import adjacency_matrix_to_mol as adj
import preprocess as pp
import test_graph_gen as gen

DATA_DIR = "./data"
OUTPUT_DIR = "./output"
os.makedirs(OUTPUT_DIR, exist_ok = True)

EPOCHS = 1000
BATCH_SIZE = 32
SAVE_EVERY = 25

data = pd.read_csv("gan-auto-train/data/qm9.csv")

