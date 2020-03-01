from tensorsHelper import pre_process_pdb_into_vectors, normalize_vector_of_int

pdb_vectors = pre_process_pdb_into_vectors()
print(pdb_vectors)

pdb_vectors_normalized = []
for pdb_vector in pdb_vectors:
    pdb_vec_normalized = normalize_vector_of_int(pdb_vector)
    print(pdb_vec_normalized)
    pdb_vectors_normalized.append(pdb_vec_normalized)

print(pdb_vectors_normalized)

