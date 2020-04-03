import os
import os.path
import re
from modeller import *
from modeller.scripts import complete_pdb
import pylab


def get_file_extension(filename):
    """Return the extension of a file.
        Arguments:
         - filename - string, the file to extract the extension

    """
    # If we want to get compound extensions, we just need to add it in the reg exp
    return re.compile(r'^.*?[.](?P<ext>fa\.gz|fasta\.gz|vcf\.gz|\w+)$').match(filename).group('ext')


def get_files(input_path, allowed_formats):
    """Return a list of path files.
        Arguments:
         - input_path - string, the path where will look for the files
         - allowed_formats - set, the allowed formats for the files
    """
    result = []
    folder = "."
    if os.path.isdir(input_path):
        folder = input_path
        result = [os.path.abspath(os.path.join(folder, f)) for f in os.listdir(folder) if
                  os.path.isfile(os.path.join(folder, f)) and get_file_extension(
                      os.path.join(folder, f)) in allowed_formats]
        result.sort()
    else:
        filename = os.path.join(folder, input_path)
        if os.path.isfile(filename) and get_file_extension(filename) in allowed_formats:
            result = [os.path.abspath(filename)]
    return result

def get_ca_atoms(chain):
    """
    Returns the CA atoms of a chain, or all the atoms in case it doesn't have CA atoms.
    Arguments:
        - chain - Bio.PDB.Chain, the chain to get the atoms
    """
    result = []
    if chain.child_list[0].has_id("CA"):
        for fixed_res in chain:
            if fixed_res.has_id("CA"):
                result.append(fixed_res['CA'])
    else:
        result = list(chain.get_atoms())
    return result


def output_print(message, verbose):
    if(verbose):
        print(message)

def make_profile(model,output_file):
    env = environ()
    env.libs.topology.read(file='$(LIB)/top_heav.lib')  # read topology
    env.libs.parameters.read(file='$(LIB)/par.lib')  # read parameters

    # read model file
    mdl = complete_pdb(env, output_file)

    # Assess with DOPE:
    s = selection(mdl)  # all atom selection
    s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=output_file+".profile",
                  normalize_profile=True, smoothing_window=15)

def show_profile(profile_file):
    """Read `profile_file` into a Python array"""
    # Read all non-comment and non-blank lines from the file:
    f = open(profile_file)
    profile_vals = []
    for line in f:
        if not line.startswith('#') and len(line) > 10:
            spl = line.split()
            profile_vals.append(float(spl[-1]))
    # Add a gap at position '0', so that we effectively count from 1:
    profile_vals.insert(0, None)

    # Plot the template and model profiles in the same plot for comparison:
    pylab.figure(1, figsize=(10,6))
    pylab.xlabel('Residue position')
    pylab.ylabel('DOPE per-residue score')
    pylab.plot(profile_vals, color='red', linewidth=2, label='Model')
    pylab.legend()
    pylab.savefig(profile_file+".png", dpi=65)

