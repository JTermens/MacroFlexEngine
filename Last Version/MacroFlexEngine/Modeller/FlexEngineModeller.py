import os
import re
from Bio import SeqIO
from helper import downloadPDBFiles
from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb

from MacroFlexEngine.lib.utils import make_profile, show_profile

log.verbose()
base_dir = 'aux/'
pdb_pir_file = 'aux/pdb_95.pir'
pdb_bin_file = 'aux/pdb_95.bin'
profile_file = 'aux/build_profile.prf'
profile_file_out = 'aux/build_profile.out'
profile_aln = 'aux/build_profile.ali'


class FlexEngineModeller:
    """
    Creates models from incomplete or non existent proteins structures, using the sequence file only
    """

    def __init__(self, file):
        """
        Parameters
        ----------
        file : str
            file in PIR format containing the sequence target
        Attributes
        ----------
        __sequenceFile : str
            file in PIR format containing the sequence target
        __model_alignment : str
            file for the alignment created to the template selected
        __templates : str
            list containing the best scored templates
        __unique_template : str
            template selected after dendrogram creation, assisting for crystallographic quality and sequence similarity
        __env : modeller.environ
            modeller environment for initial creation of DBs
        __sdb : modeller.sequence_db
            modeller object containing a database of sequences
        """
        self.__sequenceFile = self.__convertFASTAtoPIR(file)
        self.__sequenceCode = os.path.splitext(file)[0]
        self.__model_alignment = None
        self.__templates = []
        self.__unique_template = None

        self.__env = environ()
        self.__sdb = None
        self.__initializePDB()

    def __convertFASTAtoPIR(self, file):
        """Convert a fasta file into a PIR file format
        Arguments:
            - file - the fasta file you want to convert
        """
        output_pir_file = f"{os.path.splitext(file)[0]}.pir"
        SeqIO.convert(file, "fasta", output_pir_file, "pir")

        return output_pir_file

    def __initializePDB(self):
        """Initializes PDB databases needed for the following processes
        """
        self.__sdb = sequence_db(self.__env)
        self.__sdb.read(seq_database_file=pdb_pir_file, seq_database_format='PIR',
                        chains_list='ALL', minmax_db_seq_len=(30, 4000), clean_sequences=True)
        self.__sdb.write(seq_database_file=pdb_bin_file, seq_database_format='BINARY',
                         chains_list='ALL')
        self.__sdb.read(seq_database_file=pdb_bin_file, seq_database_format='BINARY',
                        chains_list='ALL')

    def __getProfile(self, matrix='blosum62', evalue=0.01):
        """Creates a profile file for the alignments made with the target sequence
        Parameters
        ----------
        matrix : str, optional
            The matrix used. Default to BLOSUM62
        evalue : float, optional
            Maximum alignment cut-off. Default to 0.01
        """
        aln = alignment(self.__env)
        aln.append(file=self.__sequenceFile, alignment_format='PIR', align_codes='ALL')

        prf = aln.to_profile()
        prf.build(self.__sdb, matrix_offset=-450, rr_file='${LIB}/' + matrix + '.sim.mat',
                  gap_penalties_1d=(-500, -50), n_prof_iterations=1,
                  check_profile=False, max_aln_evalue=evalue)
        prf.write(file=profile_file, profile_format='TEXT')
        aln.write(file=profile_aln, alignment_format='PIR')

    def __preProcessProfileFile(self):
        """Process profile file to a format easy to read and interpret for later use
       """
        os.system("cat " + profile_file + " | tail -n +8 | awk '{print $2, $12}' > " + profile_file_out)

    def __getBestTemplatesFromProfile(self):
        """Obtains a list of the better templates from the profile file
        """
        self.__preProcessProfileFile()

        myFile = open(profile_file_out, 'r')
        myFileLines = myFile.readlines()

        for line in myFileLines:
            content = line.split(' ')
            if (content[1].replace('\n', '') == '0.0'):
                self.__templates.append(content[0])

    def __performSequenceIdentityComparison(self):
        """Downloads the PDBs required to perform, then a structure comparison with them
        """
        downloadPDBFiles(self.__templates, base_dir)
        local_env = environ()
        aln = alignment(local_env)

        for pdb in self.__templates:
            pdb_id = pdb[:4]
            pdb_chain = pdb[4]
            m = model(local_env, file=f'{base_dir}{pdb[:4]}', model_segment=('FIRST:' + pdb_chain, 'LAST:' + pdb_chain))
            aln.append_model(m, atom_files=pdb_id, align_codes=pdb_id + pdb_chain)

        aln.malign()
        aln.malign3d()
        aln.compare_structures()
        aln.id_table(matrix_file='aux/family.mat')
        local_env.dendrogram(matrix_file='aux/family.mat', cluster_cut=-1.0)

        self.__unique_template = input(
            'Please, select the best template to perform the model from aux/family.mat file: ')
        print(self.__unique_template + ' selected! Building the model...')

    def __alignTargetSequenceWithTemplate(self):
        """ Aligns the target sequence with the template selected by the user
        """
        env = environ()
        aln = alignment(env)
        mdl = model(env, file=f'{base_dir}{self.__unique_template}', model_segment=('FIRST:A', 'LAST:A'))
        aln.append_model(mdl, align_codes=f'{self.__unique_template}A',
                         atom_files=f'{base_dir}{self.__unique_template}.pdb')
        aln.append(file=self.__sequenceFile, align_codes=self.__sequenceCode)
        aln.align2d()

        self.__model_alignment = f'{base_dir}{self.__sequenceCode}-{self.__unique_template}A.ali'
        aln.write(file=self.__model_alignment, alignment_format='PIR')
        aln.write(file=f'{base_dir}{self.__sequenceCode}-{self.__unique_template}A.pap', alignment_format='PAP')

    def __buildModel(self,file):
        """ Builds the final Model
        """
        env3 = environ()
        a = automodel(env3, alnfile=self.__model_alignment,
                      knowns=f'{self.__unique_template}A', sequence=self.__sequenceCode,
                      assess_methods=(assess.DOPE,
                                      # soap_protein_od.Scorer(),
                                      assess.GA341))
        a.starting_model = 1
        a.ending_model = 5
        a.make()

        #change output folders
        pattern = re.compile(r"^(.+)\/")
        fasta_folder = pattern.match(file)

        pattern = re.compile(r"\/(.+)")
        filename = pattern.match(file)
        
        os.system(f"cd {fasta_folder} && mkdir modeller && mkdir pdbs")
        os.system(f"mv *.pdb pdbs && mv {filename}.* modeller")

    def beginProcess(self):
        self.__getProfile()
        self.__getBestTemplatesFromProfile()
        self.__performSequenceIdentityComparison()
        self.__alignTargetSequenceWithTemplate()
        self.__buildModel()


def model_fasta(fasta_folder, pdb_folder):
    """ Main function to start processing fasta files from the given folder, and outputting them on the give output folder
    Arguments:
       - fasta_folder - the folder containing all the fasta files
       - pdb_folder - the output folder
    """
    fasta_files = utils.get_files(input_path=fasta_folder, allowed_formats={"fasta","FASTA","fa"})
    if not fasta_files:
        raise FileNotFoundError(f"Folder {fasta_folder} not found or no fasta files found in it")

    for fasta_file in fasta_files:
        builder = FlexEngineModeller(fasta_file)
        builder.beginProcess()

    confirm_selection = False
    while !confirm_selection:
        selected_pdb = input('Please, select the best model to build the complex: ')

        if selected_pdb+".pdb" not in fasta_folder+"/pdbs/":
            print(f"WARNING: {selected_pdb} not found, chose another one")
            continue

        make_profile(fasta_folder+"/pdbs/"+selected_pdb+".pdb", fasta_folder+"/pdbs")
        show_profile(fasta_folder+"/pdbs/"+selected_pdb+".profile")

        confirmation = input(f"Do you confirm the {selected_pdb}.pdb file (y/n): ")
        if confirmation == 'y':
            confirm_selection = True

    os.system(f"cp {fasta_folder}/pdbs/{selected_pdb}.pdb {pdb_folder}")
    print(f"Selection confirmed, pdb added to {pdb_folder}")
