import Bio.PDB


parser = Bio.PDB.PDBParser()
structure = parser.get_structure('4g83', "data/pdb/4g83.pdb")

#iterate through pdb 

for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                print (atom)

#separate pdb by chains 

for chain in structure.get_chains():
	io = Bio.PDB.PDBIO()
	io.set_structure(chain)
	io.save(structure.get_id()+chain.get_id()+'.pdb')
