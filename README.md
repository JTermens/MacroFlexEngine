# PPI prediction project

## Options

* Modeller, IMP, X3DNA
* Structural superimposition – random sampling
* Deep learning

## Papers to read

* [Exploring the protein-docking problem using convolutional neural networks](https://upcommons.upc.edu/handle/2117/115303): TFG, muy largo
* [Tensorflow Based Deep Learning Model and Snakemake Workflow for Peptide-Protein Binding Predictions](https://www.biorxiv.org/content/10.1101/410928v3)
* [Comparing two deep learning sequence-based models for protein-protein interaction prediction](https://arxiv.org/abs/1901.06268)
* [OnionNet: a Multiple-Layer Intermolecular-Contact-Based Convolutional Neural Network for Protein–Ligand Binding Affinity Prediction](https://pubs.acs.org/doi/full/10.1021/acsomega.9b01997#)
* [Predicting protein–protein interactions through sequence-based deep learning](https://academic.oup.com/bioinformatics/article/34/17/i802/5093239)
* [Protein-Protein Interaction Interface Residue Pair Prediction Based on Deep Learning Architecture](https://ieeexplore.ieee.org/document/7932134)

## Packages

### IMP 
* [reference guide documentation](https://integrativemodeling.org/2.12.0/doc/ref/), parece ser que IMP ya utiliza modeller de por si.
* [IMP tutorial](https://integrativemodeling.org/2.4.0/doc/tutorial/library.html) 
* Units: angstrom for all distances,kcal/molÅ for forces/derivatives,kcal/mol for energies,radians for angles,all angles are counterclockwise,all charges are in units of the elementary charge,all times are in femtoseconds

## PDB
* [PDB format structure](https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html)
* [BioPDB python module documentation](https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ)


Las compañeras de segundo me han dicho (literalmente): Nosotras al final desarrollamos un algoritmo bastante chulo (en mi opinión) pero lo único que al final importó fue que funcione con estructuras cuanto más complejas y que lo haga en un tiempo razonable tipo un par de horas (por Baldo) y que se usen todas las herramientas aprendidas en la asignatura de python (para Javi).

### Deep Learning approach
As far as I read, the paper that came out to be the bette solution above all was the one using as a starting point PSSM profiles obtained from PSI-BLAST. This can be a possibility, we need to manage data in a different manner to use it in a deep learning network (I suggest use Keras, as provides a framework easy to use that covers the majority of different kinds of NN [Neural Networks]). 

For our particular case in which we will be using PDB files as inputs, I was thinking the following:

#### Training the Model

##### OPTION 1
1) We can try to train our network with PDB files from the separated proteins that forms the complex (protein and ligand) and the PDB of the final complex
2) With that information we need to create tensors (vector information representation) to input our model training in a way that makes sense (could be matrices of 1s and 0s, it's the most used solution)
3) Separate sets of test data and training data *which can not be the same*, using training data to test the model (or viceversa) can cause overfitting 

##### OPTION 2 (https://github.com/hashemifar/DPPI)
1) Obtained the 2 separated sequences from each part of the complex (i.e. the protein and ligand) construct a tensor with a fixed length (the length of the longer sequence, we need to define which can be a good cut off) x 20 (the amount of aminoacids) 
2) To represent each of the peptide sequences, by using the PSSM files generated from PSI-BLAST with the option **-out_ascii_pssm** to obtain the specific sequence profile [like this one](./Files/ascii_target.pssm)
3) Clean-up a bit the file (removing headers, amino acid references of columns and footer information) before use it for training
4) As we can not use tensors with integers, the best approach is to transform them (vectorization) in a improved representative manner
