from __future__ import print_function

import numpy as np

import rmsd

from Deprecated.tensorsHelper import pre_process_pdb_into_vectors


def rotation_matrix(sigma):
    """
    https://en.wikipedia.org/wiki/Rotation_matrix
    """

    radians = sigma * np.pi / 180.0

    #r11 = 1
    #r12 = 0
    #r13 = 0

    #r21 = 0
    #r22 = np.cos(radians)
    #r23 = -np.sin(radians)

    #r31 = 0
    #r32 = np.sin(radians)
    #r33 = np.cos(radians)

    r11 = np.cos(radians)
    r12 = -np.sin(radians)
    r13 = 0

    r21 = np.sin(radians)
    r22 = np.cos(radians)
    r23 = 0

    r31 = 0
    r32 = 0
    r33 = 1

    R = np.array([[r11, r12,r13], [r21, r22,r23], [r31, r32,r33]])

    return R


pdb_vectors = pre_process_pdb_into_vectors()

vector1 = np.asarray(pdb_vectors[0])
vector2 = np.asarray(pdb_vectors[1])

A = vector1[:,6:9]
B = vector2[:,6:9]


# Rotate
B = np.dot(B, rotation_matrix(90))

print("Normal RMSD", rmsd.rmsd(A, B))


# Manipulate
A -= rmsd.centroid(A)
B -= rmsd.centroid(B)

print("Translated RMSD", rmsd.rmsd(A, B))


U = rmsd.kabsch(A, B)
A = np.dot(A, U)

print("Rotated RMSD", rmsd.rmsd(A, B))
