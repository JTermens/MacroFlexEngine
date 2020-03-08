from tensorsHelper import pre_process_pdb_into_vectors, normalize_vector_of_int

import numpy as np
import os
from tensorflow import keras
from keras import layers
import matplotlib.pyplot as plt
import tensorflow as tf

protein_example_dir = './Files/PDBs/Chains/1a4e'
protein_example_output_dir = './Files/PDBs/Chains/1a4m'
amino_acids_dict = {'-': -1, 'A': 0, 'R': 2, 'N': 2, 'D': 3, 'C': 4, 'E': 5, 'Q': 6, 'G': 7, 'H': 8, 'I': 9, 'L': 10,
                    'K': 11, 'M': 12, 'F': 13, 'P': 14, 'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19}

pssm = []
pssm_with_seq = []
for a_file in os.listdir(protein_example_dir):
    if a_file.endswith('.pssm'):
        with open(f'{protein_example_dir}/{a_file}', 'r') as f:
            pssm = []
            for line in f:
                aux = line.split()
                aux_float = [float(i) for i in aux[1:21]]
                aux_float = [float(i) / max(aux_float) for i in aux_float]
                aux_float.append(float(amino_acids_dict.get(aux[0])) / max(aux_float))
                pssm.append(aux_float)
            pssm_with_seq.append(pssm)

pssm_np = np.array(pssm, dtype=float)  # shape (488,21) (files, residues, pssm for each residue)

pdb_vectors = pre_process_pdb_into_vectors(protein_example_dir, True)
pdbs_np = np.array(normalize_vector_of_int(pdb_vectors),
                   dtype=float)  # shape (15728, 12) (lines per file, columns for ATOM line)

model = keras.Sequential([keras.layers.Dense(units=1, input_shape=[1])])
model.compile(optimizer='sgd', loss='mean_squared_error')

model.fit(np.reshape(pssm_np, pdbs_np), pdbs_np, epochs=100)

# print("predicted Answer : ", model.predict(pssm_output))
