import argparse

parser = argparse.ArgumentParser()

# Mandatory arguments
parser.add_argument("pdb_one", help="First PDB file to model the complex", type=str)
parser.add_argument("pdb_two", help="Second PDB file to model the complex", type=str)

# Optional arguments
parser.add_argument("-f", "--fasta", help="FASTA file for uncompleted models", type=str)
parser.add_argument("-v", "--verbosity", help="Increase output verbosity", action="store_true")

args = parser.parse_args()

if args.verbosity:
    print(f"For now, echoing things... PDB 1: {args.pdb_one} -  PDB 2: {args.pdb_two} - "
          f"Optional FASTA file sequence: {args.fasta}")
else:
    print(f"PDB 1: {args.pdb_one} -  PDB 2: {args.pdb_two} - FASTA: { args.fasta }")
