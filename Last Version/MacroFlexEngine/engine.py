import copy
import numpy as np
import sys

from Bio.PDB import NeighborSearch, PDBParser, Selection, Superimposer

import MacroFlexEngine.lib.center_of_mass as cm
import MacroFlexEngine.lib.interactions as intrc

from MacroFlexEngine.lib.dictionary import chain_dict
from MacroFlexEngine.lib.utils import output_print, get_ca_atoms


class MacroFlexEngine(object):
    """Engine to create a model"""

    def __init__(self, verbose=False):
        """Constructor MacroFlexEngine class"""
        self.model = None
        self.to_evaluate = None
        self.interactions = intrc.Interactions()
        self.codes_used = []
        self.verbose = verbose

    def construct_engine(self, input_folder, max_chains=200):
        """
        Retruns a Bio.PDB.Structure composed by complex pairs.
        Arguments:
            - input_folder - string, a folder containing a list of complex pairs
            - max_chains (optional) - int, maximum number of chains
        """

        output_print("Initializing the builder engine...",self.verbose)

        self.interactions.populate_interactions(input_folder,self.verbose)

        output_print(f"Complexes found in PDBs:\n{self.interactions}",self.verbose)
        output_print(f"Computing chain homologs...",self.verbose)

        self.interactions.populate_homologous_chains(identity=1)

        #sys.exit()

        remaining_complexes = self.interactions.get_complexes_list()

        parser = PDBParser()
        chain_number = 2

        structure = None
        while chain_number == 2:
            first_complex = self.interactions.interactions[remaining_complexes.pop(0)]
            first_chains = first_complex.get_chain_list()
            self.to_evaluate = [first_chains[0], first_chains[1]]
            structure = parser.get_structure("final_structure", first_complex.filename)
            self.model = structure[0]
            self.codes_used = None

            chain_number = self.__process_chains(max_chains)
        output_print("Model ended", self.verbose)
        return structure


    def __process_chains(self, max_chains):
        """
        Add chains to a model and return the number of final chains.
        Arguments:
            - max_chains - int (optional), maximum number of chains of the model
        """
        while self.to_evaluate:
            evaluating_chain = self.to_evaluate.pop(0)

            output_print(f"Evaluating chain: {evaluating_chain.original_label} from "+\
                evaluating_chain.get_filename(), self.verbose)

            for candidate_chain in evaluating_chain.homologous_chains:
                candidate_complex = candidate_chain.parent
                complementary_candidate_chain = candidate_complex.complementary_chain(candidate_chain)

                output_print(f"\t-Candidate chain: {candidate_chain.label} form "+\
                    candidate_chain.get_filename(), self.verbose)
                output_print(f"\t-Complem.  chain: {complementary_candidate_chain.label} from "+\
                    complementary_candidate_chain.get_filename(), self.verbose)

                parser = PDBParser()
                candidate_model = parser.get_structure(candidate_complex.id, candidate_complex.filename)

                fixed_chain = self.model[evaluating_chain.label]
                target_chain = candidate_model[0][candidate_chain.original_label]
                moved_chain = candidate_model[0][complementary_candidate_chain.original_label]
                superimp_done = self.__do_superimpose(fixed_chain, target_chain, moved_chain,1.0)
                if not superimp_done:
                    continue
                if not self.__determine_clashes(moved_chain):
                    moved_chain.id = self.__next_chain_id(set(moved_chain.parent.child_dict))
                    self.model.add(moved_chain)

                    next_to_evaluate = copy.copy(complementary_candidate_chain)
                    next_to_evaluate.label = moved_chain.id
                    self.to_evaluate.append(next_to_evaluate)

                    output_print(f"Chain num {len(self.model)} added.", self.verbose)

                    if len(self.model) >= max_chains:
                        return len(self.model)
        return len(self.model)

    def __do_superimpose(self,fixed_chain,target_chain,moved_chain,max_rmsd_100):

        fixed_atoms = get_ca_atoms(fixed_chain)
        target_atoms = get_ca_atoms(target_chain)

        super_imposer = Superimposer()
        super_imposer.set_atoms(fixed_atoms, target_atoms)

        rmsd = super_imposer.rms
        N = min(len(fixed_atoms),len(target_atoms))
        rmsd_100 = rmsd/(1+np.log(np.sqrt(N/100)))

        if rmsd_100 > max_rmsd_100:
            output_print(f"WARNING: rmsd_100 = {rmsd_100} > threshold ({max_rmsd_100}), CHANGING CHAIN...",\
             self.verbose)
            return False
        
        moved_chain_atoms = []

        for residues in moved_chain:
            moved_chain_atoms.append(residues)

        super_imposer.apply(moved_chain_atoms)
        return True

    def __determine_clashes(self, chain, max_clashes=10):
        """This function returns True if the number of detected clashed is grater than max_clashes
        and False otherwise"""

        num_clashes = 0 # starts a counter
        output_print("Analysing possible chain clashes...",self.verbose)

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

    def __next_chain_id(self, excluded=None):
        """
        Returns the next avaliable chain label in the model
        Arguments:
            - excluded - set, additional set of labels to be omitted
        """
        
        self.codes_used = [x._id for x in self.model.get_chains()]

        for residue in chain_dict:
            if residue not in self.codes_used and residue not in excluded:
                return residue
        return None

    def get_model_profile(output_file,file_format):

        # read model file
        mdl = complete_pdb(output_file+file_format)

        # Assess with DOPE:
        s = selection(mdl)  # all atom selection
        s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=output_file+".profile",\
          normalize_profile=True, smoothing_window=15)