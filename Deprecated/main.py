import numpy as np
from Deprecated.tensorsHelper import pre_process_pdb_into_vectors, normalize_vector_of_int

pdb_vectors = pre_process_pdb_into_vectors()

pdb_vectors_normalized = []
for pdb_vector in pdb_vectors:
    pdb_vec_normalized = normalize_vector_of_int(pdb_vector)
    pdb_vectors_normalized.append(pdb_vec_normalized)

np_pdb_normalized_vector = np.array(pdb_vectors_normalized)
print(np_pdb_normalized_vector.ndim)