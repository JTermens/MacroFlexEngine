import numpy as np
from Bio.PDB import NeighborSearch, PDBParser, Selection, Superimposer
from Bio.PDB.Polypeptide import PPBuilder

import center_of_mass as cm


class ProtComplex(object):
    """Class to model a protein & DNA/RNA macrocomplex (cuaternary structure)"""
    def __init__(self, verbose=False):
        """Constructor of the class"""
        self.model = None
        self.verbose = verbose

    def is_clashed(self, chain, max_clashes=10):
        """This function returns True if the number of detected clashed is grater than max_clashes 
        and False otherwise"""

        num_clashes = 0

        if verbose:
            print("Analysing possible chain clashes")

        model_atoms = Selection.unfold_entities(list(self.model.get_atoms()), 'A')
        neighbor_search = NeighborSearch(model_atoms)

        center = np.array(cm.center_of_mass_chain(chain, 'ATOM'), dtype=float)
        max_radius = max([self.atom_distance(center, atom) for atom in chain.get_atoms()])

        close_atoms_chain = neighbor_search.search(center, max_radius+2) #add two amstrongs to the max radius

        if close_atoms_chain:
            if verbose:
                print(len(close_atoms_chain)," clash candidates found")

            neighbor_search = NeighborSearch(close_atoms_chain)

            for residue in chain:
                for atom in residue:
                    atom_coord= tuple(atom.get_coord())
                    close_atoms_atom = neighbor_search.search(atom_coord, 0.8) #CH distance 1.09
                    if close_atoms_atom:
                        num_clashes +=1
                        if num_clashes >= max_clashes:
                            if verbose:
                                print("WARNING: ",num_clashes," clashes found, CHANGEING CHAIN")
                            return True
        return False
