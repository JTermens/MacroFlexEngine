
import Bio.PDB
import Bio.pairwise2
from Bio.PDB.Polypeptide import PPBuilder


def align(chain1,chain2):

	ppb = PPBuilder()
	for pp in ppb.build_peptides(chain1):
		sequence1 = pp.get_sequence()	
	for pp in ppb.build_peptides(chain2):
		sequence2 = pp.get_sequence()
	alignment = Bio.pairwise2.align.globalxx(sequence1, sequence2) # just proteins
	return alignment

def superimposition(homo_chains,structure1,structure2):

	sp1 = []
	sp2 = []

	for model in structure1:
		for chain in model:
			if chain.get_id() == homo_chains[0][0].id:
					for residue in chain:
						for atom in residue:
							sp1.append(atom)

	for model in structure2:
		for chain in model:
			if chain.get_id() == homo_chains[0][1].id:
				for residue in chain:
					for atom in residue:
						sp2.append(atom)

	super_imposer = Bio.PDB.Superimposer()
	super_imposer.set_atoms(sp1,sp2)
	super_imposer.apply(structure2.get_atoms())

	# Save the aligned version of 1UBQ.pdb
	io = Bio.PDB.PDBIO()
	io.set_structure(structure2) 
	io.save("models/s2.pdb")
	io.set_structure(structure1) 
	io.save("models/s1.pdb")

	return(super_imposer.rms)


if __name__ == "__main__":

	parser = Bio.PDB.PDBParser()
	structure1 = parser.get_structure('6om3A-B', "input/6om3A-B.pdb")
	structure2 = parser.get_structure('6om3B-C', "input/6om3B-C.pdb")

	threshold = 0.9
	homo_chains = []

	for chain1 in (structure1.get_chains()):
		for chain2 in (structure2.get_chains()):
			Alignment = align(chain1,chain2)
			score = Alignment[0][2] / len(Alignment[0][0])
			if score > threshold:
				homo_chains.append((chain1,chain2))
	print(Alignment)

#Found no similar chain:
	if len(homo_chains) == 0:
		print("No chains are similar enough to be superimposed")

#Only found one superimposable chain:
	elif len(homo_chains) == 1:
		print(superimposition(homo_chains,structure1,structure2))

#Found more than one similar chain we are superimposing two homodimers
	else:
		pass
		



	
