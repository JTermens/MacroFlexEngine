import copy
import numpy as np

from Bio.PDB import NeighborSearch, PDBParser, Selection

import Final.lib.center_of_mass as cm
import Final.lib.interactions as intrc

from Final.lib.dictionary import chain_dict
from Final.lib.utils import output_print


class MacroFlexEngine(object):
    """Engine to create a model"""

    def __init__(self, verbose=False):
        """Constructor MacroFlexEngine class"""
        self.model = None
        self.to_evaluate = None
        self.interactions = intrc.ChainInteractions()
        self.codes_used = None
        self.verbose = verbose

    def construct_engine(self, input_folder, max_chains=200):
        """
        Retruns a Bio.PDB.Structure composed by complex pairs.
        Arguments:
            - input_folder - string, a folder containing a list of complex pairs
            - max_chains (optional) - int, maximum number of chains
        """

        #------------------------------------------------------------------------------------------
        # PROPUESTA IMPLEMENTACION MODELER
        # 
        # if fasta files in input folder:
        #   for each fasta file:
        #       identify chains without structure -> run modeler
        #       identify chains with structure
        #       populate interactions with the newly modeled chains and the ones with structure
        #------------------------------------------------------------------------------------------

        self.interactions.populate_interactions(input_folder)
        self.interactions.populate_homologous_chains(score_limit=9.8)

        remaining_complexes = self.interactions.get_complexes_list()

        parser = PDBParser()
        chain_number = 2

        output_print("Initializing the builder engine...",verbose)

        structure = None
        while chain_number == 2:
            first_complex = self.interactions.interactions[remaining_complexes.pop(0)]
            first_chains = first_complex.get_chain_list()
            self.to_evaluate = [first_chains[0], first_chains[1]]
            structure = parser.get_structure("final_structure", first_complex.filename)
            self.model = structure[0]
            self.codes_used = None

            chain_number = self.__process_chains(max_chains)
        output_print("Model ended", verbose)
        return structure


    def __process_chains(self, max_chains):
        """
        Add chains to a model and return the number of final chains.
        Arguments:
            - max_chains - int (optional), maximum number of chains of the model
        """
        while self.to_evaluate:
            evaluating_chain = self.to_evaluate.pop(0)

            for candidate_chain in evaluating_chain.homologous_chains:
                candidate_complex = candidate_chain.parent
                complementary_candidate_chain = candidate_complex.complementary_chain(candidate_chain)

                parser = PDBParser()
                candidate_model = parser.get_structure(candidate_complex.id, candidate_complex.filename)

                fixed_chain = self.model[evaluating_chain.label]
                target_chain = candidate_model[0][candidate_chain.original_label]
                moved_chain = candidate_model[0][complementary_candidate_chain.original_label]
                # -----------------------------------------------------------------------------
                # TODO: Here we need to connect Miguel logic for superimposition, proposed connection commented below
                # self.__Miguel_superimposition(fixed_chain, target_chain, moved_chain)
                self.do_superimpose(fixed_chain, target_chain, moved_chain)
                # -----------------------------------------------------------------------------

                # -----------------------------------------------------------------------------
                # RESOLVED! TODO: Here we need to connect Joan logic for clashes, proposed connection commented below
                # if not self.__determine_clashes(moved_chain, limit=10):
                if not self.__determine_clashes(moved_chain):
                # -----------------------------------------------------------------------------
                    moved_chain.id = self.__move_next_residue({target_chain.id})
                    self.model.add(moved_chain)

                    next_to_evaluate = copy.copy(complementary_candidate_chain)
                    next_to_evaluate.label = moved_chain.id
                    self.to_evaluate.append(next_to_evaluate)

                    output_print(f"Chain num {len(self.model)} added.", verbose)

                    if len(self.model) >= max_chains:
                        return len(self.model)
        return len(self.model)

    def __do_superimpose(self,fixed_chain,target_chain,moved_chain):
        # Here Migue should put his adapted code to handle 3 chains as input

        fixed_atoms = []
        target_atoms = []  
        
        for atom in fixed_chain.get_atoms():
            fixed_atoms.append(atom)

        for atom in target_chain.get_atoms():
            target_atoms.append(atom)

        if len(fixed_atoms) > len(target_atoms):
            fixed_atoms = fixed_atoms[0:len(target_atoms)]

        if len(fixed_atoms) < len (target_atoms):
            target_atoms = target_atoms[0:len(fixed_atoms)]

        super_imposer = Superimposer()
        super_imposer.set_atoms(fixed_atoms, target_atoms)
        
        moved_chain_atoms = []

        for residues in moved_chain:
            moved_chain_atoms.append(residues)

        super_imposer.apply(moved_chain_atoms)

    # -------------------Joan logic from chain_clashed.py
    def __determine_clashes(self, chain, max_clashes=10):
        """This function returns True if the number of detected clashed is grater than max_clashes
        and False otherwise"""

        num_clashes = 0 # starts a counter
        output_print("Analysing possible chain clashes...")

        model_atoms = Selection.unfold_entities(list(self.model.get_atoms()), 'A')
        neighbor_search = NeighborSearch(model_atoms)

        center = np.array(cm.center_of_mass_chain(chain, 'ATOM'), dtype=float)
        max_radius = max([self.__atom_distance(center, atom) for atom in chain.get_atoms()])

        close_atoms_chain = neighbor_search.search(center, max_radius+2) #add two amstrongs to the max radius

        if close_atoms_chain:
            output_print(f"{len(close_atoms_chain)} clashed candidates found", self.verbose)

            neighbor_search = NeighborSearch(close_atoms_chain)

            for residue in chain:
                for atom in residue:
                    atom_coord= tuple(atom.get_coord())
                    close_atoms_atom = neighbor_search.search(atom_coord, 0.8) #CH distance 1.09
                    if close_atoms_atom:
                        num_clashes +=1
                        if num_clashes >= max_clashes:
                            output_print(f"WARNING: {num_clashes} clashes found, CHANGING CHAIN...", self.verbose)
                            return True
        return False

    def __atom_distance(self, coord, atom):
        """Returns the distance between a coordinate and an atom
        Arguments:
            - coord - np.Array, a float numpy array
            - atom - Bio.PDB.Atom, the output dir where the file is created
        """
        return np.linalg.norm(atom.coord - coord)

    def __move_next_residue(self, excluded=None):
        """Returns the next avaliable chain label in the model
        Arguments:
            - excluded - set, additional set of labels to be omitted
        """
        if(self.codes_used == None):
            self.codes_used = [x._id for x in self.model.get_chains()]

        for residue in chain_dict:
            if residue not in self.codes_used and residue not in excluded:
                return residue
        return None