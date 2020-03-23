import os
from Bio.PDB import PDBList


def downloadPDBFiles(pdbs, dir):
    """Downloads PDB files from a list of PDB codes.

    Parameters
    ----------
    pdbs : list of str
       list containing all PDB codes (str) for download
    dir : str
       output directory
    """
    pdb_list = PDBList()
    for i in pdbs:
        pdb_id = i[:4]
        pdb_list.retrieve_pdb_file(pdb_id, file_format="pdb", pdir=dir, obsolete=False)

    # Delete unused folder obsolete
    os.system("rm -d obsolete")

    files = os.listdir(dir)
    for index, file in enumerate(files):
        if file.endswith(".ent"):
            new_name = file.replace('pdb', '')
            os.rename(os.path.join(dir, file), os.path.join(dir, f'{new_name.split(".")[0]}.pdb'))
