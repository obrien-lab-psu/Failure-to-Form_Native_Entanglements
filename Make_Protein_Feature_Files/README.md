# Prepare protein feature files and make gene lists 

## General workflow
```mermaid
graph TD
    A[Gen Feature Files] --> B[Make Gene lists]
    click A "https://github.com/obrien-lab-psu/Failure-to-Form_Native_Entanglements/tree/main/Make_Protein_Feature_Files#generate-feature-files"
    click B "https://github.com/obrien-lab-psu/Failure-to-Form_Native_Entanglements/tree/main/Make_Protein_Feature_Files#generate-gene-lists-for-all-other-analysis"
``` 

## Generate Feature Files
Here we collate information from various databases to make a static feature file for each protein in our dataset(s). 
Structural classification [SCOP](https://www.ebi.ac.uk/pdbe/scop/)  
Protein essentiality [DEG](http://origin.tubic.org/deg/public/index.php) - MG1655 II [[1]](#1)  
Disordered content [DisProt](https://disprot.org/)   
Native entanglements found by [user](https://github.com/obrien-lab-psu/Failure-to-Form_Native_Entanglements/tree/main/Native_Entanglements_in_PDBs)  
Secondary structure content [STRIDE](https://webclu.bio.wzw.tum.de/stride/)  
Local changes in conformation observed in LiPMS experiments  


### Usage of [Gen_proteome_features.py](src/data/Gen_proteome_features.py)
```
usage: Gen_proteome_features.py [-h] -f CONTROL_FILE -o OUTPATH -u UNIPROT_ID -p PDB_ID -c CHAIN -l LOG_FILE -m LIPMS_FILES -s SKIP_CONTACT_LIB

Process user specified arguments

options:
  -h, --help            show this help message and exit
  -f CONTROL_FILE, --control_file CONTROL_FILE
                        Path to control file.
  -o OUTPATH, --outpath OUTPATH
                        path to outdir
  -u UNIPROT_ID, --uniprot_id UNIPROT_ID
                        Uniprot Ascession ID
  -p PDB_ID, --pdb_id PDB_ID
                        PDB ID to process
  -c CHAIN, --chain CHAIN
                        Chain of the PDB to use
  -l LOG_FILE, --log_file LOG_FILE
                        Path to logging file
  -m LIPMS_FILES, --lipms_files LIPMS_FILES
                        Path to lipms files to use
  -s SKIP_CONTACT_LIB, --skip_contact_lib SKIP_CONTACT_LIB
                        True: skip contact lib generation | False: do not skip
```
The control file has the following format.

```
pdb_dir = path-to-slug/Native_Entanglements_in_PDBs/Representative_experimental_structures/Ecoli/PDBs/ 
mapping_dir = path-to-slug/Native_Entanglements_in_PDBs/Representative_experimental_structures/Ecoli/PDBs/
clustered_ent_files = path-to-slug/Native_Entanglements_in_PDBs/Entanglements/Ecoli/EXP/mapped_NoSlipKNots_clustered_GE/
scope_defs_file = data/cop-cla-defs.txt
scope_file = data/cop-cla-latest.txt
essential_genes_file = data/deg_annotation_p.csv
path_to_disprot = data/DisProt_release_2023_12_with_ambiguous_evidences.tsv
path_to_contact_pt = data/*_contact_potential.csv
```
| key | description |
|----------|----------|
| pdb_dir   | path to PDB directory   |
| mapping_dir    | path to mapping directory  |
| clustered_ent_files | path to clustered entanglments files |
| path_to_contact_pt | path to contact potentials |
| scope_defs_file | path to SCOP definition file |
| scope_file | path to SCOP data file |
| essential_genes_file | path to Deg database file |
| path_to_disprot | path to DisProt database file |


If you have the [SLUG] then you can use the command files located [here](src/command_lists/Gen_proteome_features_EXP_FLiPPR.cmds) to reproduce the features files used in this work for the experimental structures and [here](src/command_lists/Gen_proteome_features_AF_FLiPPR.cmds) for the AlphaFold structures. Please modify any other pathing as necessary. 

### output file descriptions
#### Protein feature files located in outpath/res_features_lib/
A "|" separated file with 1 row per residue in the structure

| Name                  | Description   |
|-----------------------|---------------|
| gene                  | Uniprot accession ID              |
| pdb                   | Protein databank ID              |
| chain                 | Chain ID from PDB file              |
| uniprot_length        | Canonical length of protein in AA              |
| essential             | Protein essential or not [[1]](#1)             |
| ent_present           | Native entanglement present              |
| pdb_resid             | PDB residue ID              |
| resname               | Residue type in 3 letter code             |
| AA                    | Residue type in 1 letter code              |
| nearest_neighbors     | Residue IDs of those residues with alpha carbons within 8A of this residue              |
| num_nearest_neighbors | # of nearest neighbors to this residue              |
| region                | Residue was in an entangled region              |
| ent_idx               | Index of the entangled region the residue was in              |
| res_sasa              | Solvant accessible surface area              |
| median_sasa           | Median solvant accessible surface area of this residue and its nearest neighbors             |
| NC                    | 0 = not a unique entanglement native contact <br> 1 = a unique entanglement native contact <br> 2 = within the +/- 3 buffer of the native contact              |
| crossing              | 0 = not a unique entanglement crossing residue <br> 1 = a unique entanglement crossing residue <br> 2 = within the +/- 3 buffer of the crossing residue              |
| mapped_resid          | Mapped residue ID in the canonical FASTA sequence              |
| secondary_struct      | Secondary structure determined by STRIDE              |
| SCOP_class            | SCOP classificiation              |
| IDR                   | In an IDR determined by DisProt              |
| cut_str               | A "/" separated string containing detailed information about which LiPMS experiments this residue was determined to have a signficant change in conformation              |
| cut_C_Rall            | A cut was observed at any timepoint in the cyto-serum only buffer LiPMS experiments              |
| cut_CD_Rall           | A cut was observed at any timepoint in the cyto-serum + DnaK buffer LiPMS experiments              |
| cut_CG_Rall           | A cut was observed at any timepoint in the cyto-serum + GroEL buffer LiPMS experiments              |
| pdb_coverage          | Percent of the canonical sequence resolved in the PDB              |
| unresolved_IDR        | Depriciated - IGNORE              |
| buried                | Depriciated - IGNORE              |

#### Contact files located in outpath/contact_type{n}_lib/
A "|" separated file with 1 contact per row. There are three different types of contacts calculated:
1. Alpha carbons with in 8A of eachother  
2. any heavy atoms within 4.5A  
3. side-chain heavy atoms within 4.5A  
For each contact an estimate of the potential is reported using three different scales commonly used in litterature.  
  
Betancourt-Thirumalai (BT) scale [[2]](#2)  
Kolinski-Godzik-Skolnick (KGS) scale [[3]](#3)  
Miyazawa-Jernigan (MJ) scale [[4]](#4)  
  
#### Unique entanglement feature files located in outpath/uent_features_lib/
A "|" separated file with 1 unique entanglement per row. The definitions of the unique entanglement features are defined [here](docs/entanglement_complexity_metrics.jpg). 
  

## Generate Gene Lists for all other analysis
Here we generate all lists of genes used in this analysis by applying all combinations of several thresholds to the Feature Files.
1. We threshold the Sum of the Peptide Abundance (SPA) for each protein observed in the original Limited Proteolysis Mass Spec (LiP-MS) experiments native samples in 10-percentile incriments along the CDFs for each buffer system: [Cyto-serum](Figures/SPA_CDFs/C_Rall_CDF_vs_spa.png) [Cyto-serum+DnaK](Figures/SPA_CDFs/CD_Rall_CDF_vs_spa.png) [Cyto-serum+GroEl](Figures/SPA_CDFs/CD_Rall_CDF_vs_spa.png)
2. We threshold the coverage observed in the experiments (COV) in 10-percentile incriments along the CDFs for each buffer system: [Cyto-serum](Figures/LiPMScov_CDFs/C_Rall_CDF_vs_cov.png) [Cyto-serum+DnaK](Figures/LiPMScov_CDFs/CD_Rall_CDF_vs_cov.png) [Cyto-serum+GroEl](Figures/LiPMScov_CDFs/CD_Rall_CDF_vs_cov.png)
3. We separate genes into essential and non-essential based on the knock-out dataset reported in the [DEG](http://origin.tubic.org/deg/public/index.php) database.
4. We separate genes into refoldable and non-refoldable based on if they have atleast 1 significant change in proteolysis suseptibility after refolding. We control for FDR using the BH method proteome wide and therefore any change in structure still observed is most likely not a false postive where as others have corrected on a per protein basis and thus use a higher threshold to determine non-refoldability. We also only consider half-tryptic peptides and neglect full-tryptic peptides entirely as the interpreation of changes in structure resulting from them is still poorly concieved in the feild. 
5. We remove proteins with known Knots from the dataset as reported in [KnotProt2.0](https://knotprot.cent.uw.edu.pl/) and [AlphaKnot](https://alphaknot.cent.uw.edu.pl/). The specific database files we used are located [here](data/KnotProt2.0.txt) for experimental structrues and [here](data/AlphaKnotProt.txt) for Alphafold structures. 

### Usage of [src/data/MakeGeneLists.py](src/data/MakeGeneLists.py)
```
usage: MakeGeneLists.py [-h] -f FEATURE_FILES -o OUTPATH -t TAG -e EXPRESS_FILE -l LIPMS_COVS_FILE -b BUFFER -th THRESHOLD -tp TIMEPOINT [-m MASK] -c LIPMS_COV_THRESHOLD -k KNOTS

Process user specified arguments

options:
  -h, --help            show this help message and exit
  -f FEATURE_FILES, --feature_files FEATURE_FILES
                        Path to residue feature files
  -o OUTPATH, --outpath OUTPATH
                        Path to output directory
  -t TAG, --tag TAG     Tag for output filenames
  -e EXPRESS_FILE, --express_file EXPRESS_FILE
                        path to expression control file or Total
  -l LIPMS_COVS_FILE, --lipms_covs_file LIPMS_COVS_FILE
                        path to lipms coverage file
  -b BUFFER, --buffer BUFFER
                        Buffer system to use: C, CD, CG
  -th THRESHOLD, --threshold THRESHOLD
                        threshold for either SPA: 0 - 90 in 10 int steps
  -tp TIMEPOINT, --timepoint TIMEPOINT
                        R1min, R5min, R2hr, Rall
  -m MASK, --mask MASK  Global mask for genes. if present no gene will pass unless in mask regardless of other filters
  -c LIPMS_COV_THRESHOLD, --lipms_cov_threshold LIPMS_COV_THRESHOLD
                        threshold for either LiPMS-COV: 0 - 90 in 10 int steps
  -k KNOTS, --knots KNOTS
                        List of genes to ignore that contain knots
```

If you have the [SLUG] then you can use the command files located [here](src/command_lists/MakeGeneLists_EXP_Rall.cmds) to reproduce the gene lists used in this work for the experimental structures and [here](src/command_lists/MakeGeneLists_EXP_Rall.cmds) for the AlphaFold structures. Please modify any other pathing as necessary. 

The count statistics of each gene list can be found [here](data/Combined_genelist_counts_EXP.csv) for the experimental dataset and [here](data/Combined_genelist_counts_AF.csv) for the Alphafold datasets. 


### References  

<a id="1">[1]</a>: Baba, T., Ara, T., Hasegawa, M., Takai, Y., Okumura, Y., Baba, M., Datsenko, K. A., Tomita, M., Wanner, B. L., & Mori, H. (2006). Construction of Escherichia coli K-12 in-frame, single-gene knockout mutants: The Keio collection. Molecular Systems Biology, 2. https://doi.org/10.1038/msb4100050

<a id="2">[2]</a>: Betancourt, M. R., & Thirumalai, D. (1999). Pair potentials for protein folding: choice of reference states and sensitivity of predicted native states to variations in the interaction schemes. Protein Science : A Publication of the Protein Society, 8(2), 361. https://doi.org/10.1110/PS.8.2.361

<a id="3">[3]</a>: Skolnick, J., Jaroszewski, L., Kolinski, A., & Godzik, A. (1997). Derivation and testing of pair potentials for protein folding. When is the quasichemical approximation correct? Protein Science, 6(3), 676–688. https://doi.org/10.1002/pro.5560060317

<a id="4">[4]</a>: Miyazawa, S., & Jernigan, R. L. (1996). Residue-Residue Potentials with a Favorable Contact Pair Term and an Unfavorable High Packing Density Term, for Simulation and Threading. In J. Mol. Biol (Vol. 256).