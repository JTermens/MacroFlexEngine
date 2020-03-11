# The MIT License
# 
# Copyright (c) 2010-2016 Anders S. Christensen
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

from Bio import PDB

sp_chain = "A"

parser = PDB.PDBParser()
structure1 = parser.get_structure('4g83A-B', "pdb/rot_4g83A-B.pdb")
structure2 = parser.get_structure('4g83A-EF', "pdb/4g83A-EF.pdb")

sp1 = []
sp2 = []

for model in structure1:
    for chain in model:
        if chain.get_id() == sp_chain:
            for residue in chain:
                for atom in residue:
                    sp1.append(atom)

for model in structure2:
    for chain in model:
        if chain.get_id() == sp_chain:
            for residue in chain:
                for atom in residue:
                    sp2.append(atom)

# Now we initiate the superimposer:
super_imposer = PDB.Superimposer()
super_imposer.set_atoms(sp1, sp2)
super_imposer.apply(structure2.get_atoms())

# Print RMSD:
print(super_imposer.rms)

# Save the aligned version of 1UBQ.pdb
io = PDB.PDBIO()
io.set_structure(structure2)
io.save("s2.pdb")
io.set_structure(structure1)
io.save("s1.pdb")
