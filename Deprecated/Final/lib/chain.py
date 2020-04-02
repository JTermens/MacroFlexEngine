import Levenshtein
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder

class Chain(object):
    """Class to save the homologous chains of the chain"""

    def __init__(self, label, parent):
        """
        Creator of Chain class
        Arguments:
         - label - string, the label of the chain
         - parent - Complex, the complex instance of the parent
        """
        self.label = label
        self.original_label = label
        self.parent = parent
        self.homologous_chains = set([])

    def __len__(self):
        """Returns the number of homologous chains"""
        return len(self.homologous_chains)

    def __get_nucleotides(self, chain):
        """
        Returns the nucleotide sequence of a chain
        Arguments:
         - chain - Bio.PDB.Chain, the chain containing nucleotides
        """
        return "".join([x.resname.strip() for x in chain.child_list])

    def __get_homologous_tuples(self, interactions, identity=1.0):
        """
        Returns a tuple with the type of chain and a list of its homologous chains,
        with an identity treshold
        Arguments:
         - interactions - dict, the dict conteining all the interactions
         - identity - float, the limit filter of identity of sequences
        """
        ppb = PPBuilder()

        tuple_list = [(self.parent.filename, self.parent.id, self.original_label)]
        chain_type = "prot"

        evaluating_structure = PDBParser().get_structure(self.original_label, self.parent.filename)
        evaluating_polypeptide = ppb.build_peptides(evaluating_structure[0][self.original_label])
        evaluating_sequence = None
        if evaluating_polypeptide:
            evaluating_sequence = evaluating_polypeptide[0].get_sequence()
        else:
            chain_type = "dna/rna"
            evaluating_sequence = self.__get_nucleotides(evaluating_structure[0][self.original_label])

        for complex_item in interactions.keys():
            complex_structure = PDBParser().get_structure(complex_item, interactions[complex_item].filename)
            for label in interactions[complex_item].chain_dict.keys():
                if complex_item != self.parent.id or label != self.label:
                    chain_polypeptide = ppb.build_peptides(complex_structure[0][label])
                    chain_sequence = None
                    if chain_polypeptide:
                        chain_sequence = chain_polypeptide[0].get_sequence()
                    else:
                        chain_sequence = self.__get_nucleotides(complex_structure[0][label])
                    if Levenshtein.ratio(str(evaluating_sequence), str(chain_sequence)) >= identity:
                        tuple_list.append((interactions[complex_item].filename, complex_item, label))
        return (chain_type, tuple_list)

    def get_homologous_chains(self, interactions, score_limit=9.8):
        """
        Returns a tuple with the type of chain and a list of its homologous chains,
        with an identity treshold
        Arguments:
         - interactions - dict, the dict conteining all the interactions
         - score_limit - float, the limit for score similarity
        """
        homologous_set = set([])

        homologous_tuples = self.__get_homologous_tuples(interactions, 0.95)
        if len(homologous_tuples[1]) != 1 and homologous_tuples[0] == "prot":
            # With STAMP
            sp = stamp.STAMPParser()
            sp.create_stamp_input(homologous_tuples[1], "/tmp/")
            stamp_id = self.parent.id + "_" + self.original_label
            scores_dict = sp.get_stamp_scores(input_folder="/tmp/", chain_id=stamp_id, limit=score_limit)
            if len(scores_dict) != 0:
                homologous_list = [x[0] for x in scores_dict[stamp_id]]
                for homologous_chain in homologous_list:
                    complex_id = homologous_chain.split("_")[0]
                    chain_label = homologous_chain.split("_")[1]
                    homologous_set.add(interactions[complex_id].chain_dict[chain_label])
        else:
            # #Without STAMP
            homologous_list = self.__get_homologous_tuples(interactions, 1)
            homologous_list[1].pop(0)
            for homologous_chain in homologous_list[1]:
                complex_id = homologous_chain[1]
                chain_label = homologous_chain[2]
                homologous_set.add(interactions[complex_id].chain_dict[chain_label])
        return homologous_set