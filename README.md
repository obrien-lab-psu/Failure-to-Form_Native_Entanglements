# Failure-to-Form_Native_Entanglements

The repository contains all code necessary to reproduce the results in [PAPER PLACEHOLDER].   
As many of the source data files are very large this repo only contains the code and knowledge to run the analysis.  
To recreate the data from the paper please download the tar-ball located here [TARBALL PLACEHOLDER] and place it in a directory that will be sourced by each script.  

## Environment setup
Please install a new Miniconda environment using the provided [environment.yml](environment.yml) file.
```python
conda env create -f environment.yml -n new_env_name
conda activate new_env_name
```
if you need help installing a Miniconda you can find instructions [here](https://docs.anaconda.com/miniconda/miniconda-install/) 

## Analysis sections 

[Generating native entanglements of protein structures](Native_Entanglements_in_PDBs/)  
  
[Processing LiP-MS data to elucidate changes in structure](Processing_LiP-MS_data/)  
  
[Making protein feature files and gene lists](Make_Protein_Feature_Files/)  
  
[Association analysis between native entanglement presence and misfolding](Association_Native_Entanglements_and_Misfolding/)  
  
[Modeling odds of misfolding involving native entanglement rather than not](Modeling_Odds_of_Misfolding/)  
  
[Chaperone clients and associations with gene essentiality and entanglements](Chaperone_Client_Associations/)  
  
[Entanglement topological complexity measures and discrimination analysis](Entanglement_Topological_Complexity_and_Discrimination/)  
  
[Sequence complexity measures and discrimination analysis](Sequence_Complexity_and_Discrimination/)  
  
[Selecting candidates for CG temperature quenching simulations](Candidate_Selection_for_CG_T-quench_Sims/)  
  
[Association analysis between native entanglement presence and protein essentiality](Association_Native_Entanglements_and_Essentiality/)
  
[Simulations of native entanglement misfolding](Simulations_of_Native_Entanglement_Misfolding/)  
