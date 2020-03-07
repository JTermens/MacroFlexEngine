from __future__ import print_function

import numpy 
import Bio.PDB

from RMSD_Kabsch import *

class ResSelect(Bio.PDB.Select):

    def accept_chain(self, chain):
        if chain.get_id() == sp_chain:
            return False
        else:
            return True


parser = Bio.PDB.PDBParser()
structure1 = parser.get_structure('4g83A-B', "4g83A-B.pdb")
structure2 = parser.get_structure('4g83A-EF', "4g83A-EF.pdb")


sp_chain = "A"

sp1 = []
sd1 = []
sp2 = []
sd2 = []

for model in structure1:
    for chain in model:
        if chain.get_id() == sp_chain:
            for residue in chain:
                for atom in residue:
                     sp1.append([atom.get_vector()[0],atom.get_vector()[1],atom.get_vector()[2]])
        else:
           for residue in chain:
                for atom in residue:
                     sd1.append([atom.get_vector()[0],atom.get_vector()[1],atom.get_vector()[2]])

                      
for model in structure2:
    for chain in model:
        if chain.get_id() == sp_chain:
            for residue in chain:
                for atom in residue:
                     sp2.append([atom.get_vector()[0],atom.get_vector()[1],atom.get_vector()[2]])
        else:
            for residue in chain:
                for atom in residue:
                     sd2.append([atom.get_vector()[0],atom.get_vector()[1],atom.get_vector()[2]])

                      
sp1 = numpy.asarray(sp1)
sp2 = numpy.asarray(sp2)
sd1 = numpy.asarray(sd1)
sd2 = numpy.asarray(sd2)


#Centroid calculation 
cen1 = centroid(sp1)
cen2 = centroid(sp2)

#Rotation matrix construction
rotation_matrix = kabsch(sp1,sp2)


#Apply translation and rotation on all atoms
for model in structure1:
    for chain in model:
            for residue in chain:
                for atom in residue:
                     atom.transform(rotation_matrix,cen1)
                      
for model in structure2:
    for chain in model:
            for residue in chain:
                for atom in residue:
                    atom.transform(rotation_matrix,cen1)
                     

pdb_io = Bio.PDB.PDBIO()
pdb_io.set_structure(structure1)
pdb_io.save("s1.pdb")

pdb_io.set_structure(structure2)
pdb_io.save("s2.pdb",ResSelect())

