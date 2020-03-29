#import Levenshtein
from Bio.PDB import PDBParser
import Bio.pairwise2 as pwise
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

    def get_homologous_chains(self, interactions, identity=1):

        homologous_set = set([])

        ppb = PPBuilder()

        eval_structure = PDBParser().get_structure(self.original_label, self.parent.filename)
        eval_polypeptide = ppb.build_peptides(eval_structure[0][self.original_label])
        eval_sequence = None
        if eval_polypeptide:
            eval_chain_type = "prot"
            eval_sequence = eval_polypeptide[0].get_sequence()
        else:
            eval_chain_type = "dna/rna"
            eval_sequence = self.__get_nucleotides(eval_structure[0][self.original_label])

        for complex_id in interactions.keys():
            complex_structure = PDBParser().get_structure(complex_id, interactions[complex_id].filename)
            for chain_label in interactions[complex_id].chain_dict.keys():
                if complex_id != self.parent.id or chain_label != self.label:
                    target_polypeptide = ppb.build_peptides(complex_structure[0][chain_label])
                    target_sequence = None
                    if target_polypeptide:
                        target_chain_type = "prot"
                        target_sequence = target_polypeptide[0].get_sequence()
                    else:
                        target_chain_type = "dna/rna"
                        target_sequence = self.__get_nucleotides(complex_structure[0][chain_label])

                    if target_chain_type == eval_chain_type:
                        aln = pwise.align.globalxx(str(eval_sequence), str(target_sequence))
                        score = aln[0][2] / len(aln[0][0])
                        if score >= identity:
                            homologous_set.add(interactions[complex_id].chain_dict[chain_label])
        return homologous_set
                            
