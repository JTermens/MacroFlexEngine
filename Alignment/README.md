# Rigid Body Alignment with Kabsch Algorithm 

* Read input dimers
* Detect which chain must be superimposed (not implemented)
* Apply BioPython superimposition function 
* Output pdb with the build complex (outputs two separate files - should be only one pdb)

## Notes
* Alignment can only be done when the two chains that are superimposed are exactly the same size. 
In most cases this will be the case since they should be the same chain but a condition should be added to check for length (and only align the first n atoms)
 
