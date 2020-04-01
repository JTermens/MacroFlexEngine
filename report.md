## 1. Introduction
MacroFlexEngine is a python application that allows the construction of protein and dna macrocomplexes from pdb information of the dimers that form this complex. Furthermore if a fasta file of the whole complex is inputed it will detect which parts of the complex are given and which are missing and complete the missing parts using modeller.


## 2. Running the application
#### 2.1 Installation
The application can be installed by running the following commands:
For linux distributions and Mac ios:
```bash
python setup.py sdist
python setup.py install
```
For windows:
```bash
python setup.py bdist_wininst
```
The application uses the BioPython library and it must be installed in your local python folder in order to use it.

#### 2.2 Execution

To execute the application the following command can be run
```bash
python MFEngine-launch -i absolute_path_to_pdb_folder
```
As well as the mandatory input flag there are several optional ones:
```
-o –output: The output file containing the final complex 
-v verbose: More detailed process information
-f -fasta: FASTA file for uncompleted models
```
By default the program will output a pdb file with the whole macrocomplex. If a folder is specified the output will be saved there. 


## 3. Theoric background
The main approach for the building of the complex used in this application is superimposition. The program finds the repeated chain in the two dimers by alignment,superimposes them and then moves one of the remaining chain. This process is repeated for all the dimers inputed as the complex grows. 
#### 3.1 Data treatment
First of all the input dimer pdbs are read and the chains and interactions are saved. With these variables saved they can now be compared to the fasta file if it is inputed to check if there is a missing chain.  
#### 3.2 Modeller
If there are chains missing the structure modeller is use to obtain it. The pdb database is used to obtain templates for the missing fasta sequence using ???blastp???. Then the alignment is done between the best templates. With these alignment the model can be created with modeller.    
#### 3.3 Homolog finder
To know which chains of the two structures must be superimposed first a pairwise sequence alignment of all of the chains is done. If the identity is higher than a ???certain score??? it is considered an homolog and will be superimposed. In these process a structural alignment is also done ?????.    
#### 3.4 Superimposition

#### 3.5 Clash Checking  


## 4.Examples

Example of execution and result analysis  

	

