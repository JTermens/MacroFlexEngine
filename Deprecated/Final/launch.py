import argparse
import os
from Deprecated.Final.MacroFlexEngine import MacroFlexEngine

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This program provides a Complex Builder for a given interaction files")

    # Mandatory arguments
    parser.add_argument('-i', '--input',
                        dest="input",
                        action="store",
                        required=True,
                        help="The input folder containing the complexes files")
    parser.add_argument('-o', '--output',
                        dest="output",
                        action="store",
                        default="./final_complex",
                        help="The output file containing the final complex")
    parser.add_argument('-v', '--verbose',
                        dest="verbose",
                        action="store_true",
                        default=False,
                        help="More detailed process information")
    parser.add_argument("-f", "--fasta",
                        dest="fasta",
                        action="store_true",
                        required=False,
                        help="FASTA file for uncompleted models")

    options = parser.parse_args()
    if not os.path.isdir(os.path.abspath(options.input)):
        raise IOError("The input is not a folder")
    if os.path.isdir(os.path.abspath(options.output)):
        raise IOError("The output is not a file")

    engine = MacroFlexEngine(options.verbose)
    structure = None

    try:
        structure = engine.construct_engine(options.input)

    except KeyboardInterrupt:
        print("Creation of complex interrupted by the user...")

    if len(structure[0]) > 52:
        io = MMCIFIO()
        io.set_structure(structure)
        io.save(filepath=options.output + ".cif", select=exclude.ExcludeWaterSelect())
    else:
        io = PDBIO()
        io.set_structure(structure)
        io.save(file=options.output + ".pdb", select=exclude.ExcludeWaterSelect())
