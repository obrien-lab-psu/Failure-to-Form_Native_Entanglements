# Generating native entanglements of protein structures
This analysis deconstructs a protein structures topology into a series of loops closed by non-covalent lassos and threading segments. 
We then cluster the deconstructed topology to assign unique entanglements for each protein structure in either our set of high quality PDBs or AF(v4) structures. 
You can read more about this method in these papers [PAPER link PLACEHOLDER]. 

## General workflow 
```mermaid
graph TD
    A[Cleaning experimental PDBs] --> B[Deconstruct topology into raw entanglements]
    B -- EXP structure --> C[Map resids to Uniprot canon] --> D[Remove pure slipknots] --> E[Cluster into unique entanglements]
    B -- AF structure --> F[Remove low-quality entanglments] --> D
``` 

