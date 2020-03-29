import copy
import os
import MacroFlexEngine.lib.utils as utils
from MacroFlexEngine.lib.chain import Chain
from MacroFlexEngine.lib.complex import Complex


class Interactions(object):
    """Class to save the interactions between chains"""

    def __init__(self):
        """Creator of Interactions class"""
        self.interactions = {}

    def get_complexes_list(self):
        """Returns a list of the interactions complexes"""
        return list(self.interactions.keys())

    def __get_chains(self, input_file):
        """
        Returns a set with all the chains found in a PDB file
        Arguments:
         - input_file - string, the PDB file to look for chains
        """
        set_chains = set([])
        with open(input_file, 'r') as fh:
            for line in fh:
                if line.startswith("ATOM"):
                    set_chains.add(line[21])

        return sorted(list(set_chains))

    def populate_interactions(self, input_folder):
        """
        Populate interactions dictionary with empty homologous
        Arguments:
         - input_folder - string, folder to look at for PDB files
        """
        input_files = utils.get_files(input_path=input_folder, allowed_formats={"pdb"})
        if not input_files:
            raise FileNotFoundError(f"The folder {input_folder} not found or no pdb files found on it")

        i = 1
        for filename in input_files:
            chains = self.__get_chains(filename)
            for chain in chains:
                complex_id = f"{i}:{i+1}"
                try:
                    self.interactions[complex_id].chain_dict[chain] = Chain(chain, self.interactions[complex_id])
                except KeyError:
                    self.interactions[complex_id] = Complex(complex_id, filename)
                    self.interactions[complex_id].chain_dict[chain] = Chain(chain, self.interactions[complex_id])
            i += 2

    def populate_homologous_chains(self, identity):
        """Populate the homologous chains from the interactions dict"""
        for complex_item in self.interactions.values():
            for chain_item in complex_item.chain_dict.values():
                if not chain_item.homologous_chains:
                    homologous_set = chain_item.get_homologous_chains(self.interactions, identity)
                    homologous_set.add(chain_item)
                    for chain in homologous_set:
                        own_homologous_set = copy.copy(homologous_set)
                        own_homologous_set.remove(chain)
                        chain.homologous_chains = own_homologous_set
