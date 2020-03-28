import copy
import os
import Final.lib.utils as utils
from Final.lib.chain import Chain
from Final.lib.complex import Complex


class Interactions(object):
    """Class to save the interactions between chains"""

    def __init__(self):
        """Creator of Interactions class"""
        self.interactions = {}

    def get_complexes_list(self):
        """Returns a list of the interactions complexes"""
        return list(self.interactions.keys())

    def __get_chains(self, input_file):
        """Returns a set with all the chains found in a PDB file
        Arguments:
         - input_file - string, the PDB file to look for chains
        """
        pre_process = os.popen("cat" + input_file + " | grep ATOM | awk '{print $5}' | sort | uniq").read()
        set_chains = pre_process.split('\n')

        return set_chains[:len(set_chains) - 1]

    def populate_interactions(self, input_folder):
        """Populate interactions dictionary with empty homologous
        Arguments:
         - input_folder - string, folder to look at for PDB files
        """
        input_files = utils.get_files(input_path=input_folder, allowed_formats={"pdb"})
        if not input_files:
            raise FileNotFoundError("Folder not found or no files found for the given formats")

        i = 1
        for filename in input_files:
            chains = self.__get_chains(filename)
            for chain in chains:
                complex_id = "%d:%d" % (i, i + 1)
                try:
                    self.interactions[complex_id].chain_dict[chain] = Chain(chain, self.interactions[complex_id])
                except KeyError:
                    self.interactions[complex_id] = Complex(complex_id, filename)
                    self.interactions[complex_id].chain_dict[chain] = Chain(chain, self.interactions[complex_id])
            i += 2

    def populate_homologous_chains(self, score_limit):
        """Populate the homologous chains from the interactions dict"""
        for complex_item in self.interactions.values():
            for chain_item in complex_item.chain_dict.values():
                if not chain_item.homologous_chains:
                    homologous_set = chain_item.get_homologous_chains(self.interactions, score_limit)
                    homologous_set.add(chain_item)
                    for chain in homologous_set:
                        homologous_specific_set = copy.copy(homologous_set)
                        homologous_specific_set.remove(chain)
                        chain.homologous_chains = homologous_specific_set
