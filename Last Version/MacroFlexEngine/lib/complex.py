class Complex(object):
    """Class to save the complex file"""
    def __init__(self, id, filename):
        """Creator of Complex class
        Arguments:
         - id - string, the id of the complex
         - filename - string, the file name of the complex
        """
        self.id = id
        self.chain_dict = {}
        self.filename = filename

    def get_chain_list(self):
        """Return a list of the child chains"""
        return list(self.chain_dict.values())

    def complementary_chain(self, chain):
        """Return the complementary chain"""
        # We supose that every complex has two chain childs
        for chain_item in self.chain_dict.values():
            if chain_item != chain:
                return chain_item
        return None

    def __str__(self):
        complex_id = self.id
        filename = self.filename.split("/")[-1]
        complex_chains = [str(chain) for chain in self.chain_dict.values()]
        return 'complex_id: '+complex_id+', filename: '+filename+', complex_chains: '+", ".join(complex_chains)