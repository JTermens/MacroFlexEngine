# read all pdb files
# remove extra info from header
# normalize with an integer dictionary the columns that are not numbers (to integers)
# normalize all data to contain float32 values between 0 and 1 (SD = 1, Mean= 0)
import os
from Deprecated.dictionary import atoms_name_dict, residues_code_dict, symbols_dict

directory = './Files/PDBs'
files = os.listdir(directory)


def read_pdb_to_create_dict():
    atom_name_data = {}
    residues_data = {}
    element_symbol_data = {}

    for pdb_file in files:
        if pdb_file.endswith('.pdb'):
            with open(f'{directory}/{pdb_file}', 'r') as f:
                for line in f:
                    if line.startswith('ATOM'):
                        atom = line[0:4]
                        name = line[12:17]
                        residues = line[17:20]
                        element_symbol = line[76:78]
                        atom_name_data[name] = True
                        residues_data[residues] = True
                        element_symbol_data[element_symbol] = True
    atoms_name_keys = list(atom_name_data.keys())
    residues_name_keys = list(residues_data.keys())
    symbols_name_keys = list(element_symbol_data.keys())
    dict_atoms_name = {atoms_name_keys[i].strip(): i for i in range(0, len(atoms_name_keys))}
    dict_residues_code = {residues_name_keys[i].strip(): i for i in range(0, len(residues_name_keys))}
    dict_symbols = {symbols_name_keys[i].strip(): i for i in range(0, len(symbols_name_keys))}

    print(dict_symbols)
    print(dict_residues_code)
    print(dict_atoms_name)


def pre_process_pdb_into_vectors():
    pdb_files = []
    for pdb_file in files[:10]:
        if pdb_file.endswith('.pdb'):
            pdb_vectors = []
            current_pdb_file = open(f'{directory}/{pdb_file}', 'r')
            pdb_content = current_pdb_file.readlines()

            for line in pdb_content:
                if line.startswith('ATOM'):
                    pdb_vector = get_line_data_as_vector(line)
                    pdb_vectors.append(pdb_vector)

        if pdb_vectors:
            pdb_files.append(pdb_vectors)

    return pdb_files


def get_line_data_as_vector(line):
    atom = 1  # ATOM = 1, HETATM = 2, TER = 3 ...
    atom_serial_number = float(line[6:12])
    atom_name = atoms_name_dict.get(line[12:17].strip())
    residue_name = residues_code_dict.get(line[17:20].strip())
    chain = ord(line[21:22].strip())
    residue_seq_number = float(line[22:26])
    x_orthogonal = float(line[30:38])
    y_orthogonal = float(line[38:46])
    z_orthogonal = float(line[46:54])
    occupancy = float(line[54:60])
    temp_factor = float(line[60:66])
    element_symbol = symbols_dict.get(line[76:78].strip())
    pdb_vector = [atom, atom_serial_number, atom_name, residue_name, chain, residue_seq_number,
                  x_orthogonal, y_orthogonal, z_orthogonal, occupancy, temp_factor, element_symbol]

    return pdb_vector


def normalize_vector_of_int(data):
    normalized_vector = []
    for vec in data:
        normalized_vector.append([float(i) / sum(vec) for i in vec])

    return normalized_vector
