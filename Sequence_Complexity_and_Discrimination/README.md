# Sequence_Complexity_and_Discrimination  

## General workflow
```mermaid
graph TD
    A[Hydropathy Anal.] --> B[Plotting]
    C[Contact Enrich. Anal.] --> D[Plotting]
    E[OR Trend Anal.]
    click A "https://github.com/obrien-lab-psu/Failure-to-Form_Native_Entanglements/tree/main/Sequence_Complexity_and_Discrimination#hydropathy-analysis-of-loop-forming-contacts"
    click B "https://github.com/obrien-lab-psu/Failure-to-Form_Native_Entanglements/tree/main/Sequence_Complexity_and_Discrimination#usage-of-plot_hydropathy_resultspy"
    click C "https://github.com/obrien-lab-psu/Failure-to-Form_Native_Entanglements/tree/main/Sequence_Complexity_and_Discrimination#loop-forming-contact-enrichment"
    click D "https://github.com/obrien-lab-psu/Failure-to-Form_Native_Entanglements/tree/main/Sequence_Complexity_and_Discrimination#usage-of-plot_energetics_results_with_permutationpy"
    click E "https://github.com/obrien-lab-psu/Failure-to-Form_Native_Entanglements/tree/main/Sequence_Complexity_and_Discrimination#trend-analysis-in-loop-forming-contacts"
``` 

## Hydropathy Analysis of Loop Forming Contacts
| Class               | Amino Acid and Hydrophobicity Score |
|---------------------|-------------------------------------|
| Strong hydrophobic  | Ile 4.5 / Val 4.2 / Leu 3.8 / Phe 2.8 / Cys 2.5 / Met 1.9 / Ala 1.8 |
| Weak hydrophobic    | Gly -0.4 / Thr -0.7 / Ser -0.8 / Trp -0.9 / Tyr -1.3 / Pro -1.6     |
| Hydrophilic         | His -3.2 / Gln -3.5 / Glu -3.5 / Asn -3.5 / Asp -3.5 / Lys -3.9 / Arg -4.5 |  
  
[Citation](https://pubmed.ncbi.nlm.nih.gov/7108955/)  

### Usage of [compare_hydropathy.py](src/data/compare_hydropathy.py)
```
usage: compare_hydropathy.py [-h] -eg ESS_GENE_LIST -neg NONESS_GENE_LIST -l LOG_FILE -c CONTACT_FILES -e UENT_FILES -o OUTPATH -t TAG

Process user specified arguments

options:
  -h, --help            show this help message and exit
  -eg ESS_GENE_LIST, --Ess_gene_list ESS_GENE_LIST
                        path to Ess gene list
  -neg NONESS_GENE_LIST, --NonEss_gene_list NONESS_GENE_LIST
                        path to Ess gene list
  -l LOG_FILE, --log_file LOG_FILE
                        Path to logging file
  -c CONTACT_FILES, --contact_files CONTACT_FILES
                        path to native contact files
  -e UENT_FILES, --uent_files UENT_FILES
                        path to unique entanglement files
  -o OUTPATH, --outpath OUTPATH
                        path to output directory. will be made if doesnt exist
  -t TAG, --tag TAG     tag for output file
```

If you have the [SLUG] then you can use the command files located [here](src/command_lists/compare_hydropathy.cmds) to reproduce loop forming contact hydropathy analysis used in this work in the experimental data set and the AlphaFold structures. Please modify any other pathing as necessary. 
  

### Usage of [Plot_hydropathy_results.py](src/data/Plot_hydropathy_results.py)
```
usage: Plot_hydropathy_results.py [-h] -f INP_FILES -o OUTPATH

Process regression data and generate plots.

options:
  -h, --help            show this help message and exit
  -f INP_FILES, --inp_files INP_FILES
                        Input file pattern for hydropathy data.
  -o OUTPATH, --outpath OUTPATH
                        Path to output directory.
```

If you have the [SLUG] then you can use the command files located [here](src/command_lists/Plot_compare_hydropathy.cmds) to reproduce the results plots for the loop forming contact hydropathy analysis used in this work in the experimental data set and the AlphaFold structures. Please modify any other pathing as necessary. 

### Results of loop forming contact hydropathy analysis
For the set of experimentally derived structures and native entanglements  
![Experimental structure loop closing hydropathy](Figures/Hydropathy/EXP/hydropathy_results.png)  
The raw data for these plots can be found [here](data/Hydropathy/EXP/hydropathy_results.csv)  
  
For the set of Alphafold derived structures and native entanglements  
![Experimental structure loop closing hydropathy](Figures/Hydropathy/AF/hydropathy_results.png)  
The raw data for these plots can be found [here](data/Hydropathy/AF/hydropathy_results.csv)  
  
  
## Loop forming contact enrichment 

### Usage of [compare_energetics_results_with_permutation.py](src/data/compare_energetics_results_with_permutation.py)
```
usage: compare_energetics_results_with_permutation.py [-h] -f FASTA_FILES -g REP_GENE_LIST -Eg ESS_GENE_LIST -NEg NONESS_GENE_LIST -l LOG_FILE -c CONTACT_FILES -r
                                                      RESFEAT_FILES -e UENT_FILES -o OUTPATH -p NUM_PERMUTE --floor FLOOR

Process user specified arguments

options:
  -h, --help            show this help message and exit
  -f FASTA_FILES, --fasta_files FASTA_FILES
                        path to fasta files for all genes in dataset
  -g REP_GENE_LIST, --rep_gene_list REP_GENE_LIST
                        path to representative gene list used in the calculation of F_ab
  -Eg ESS_GENE_LIST, --Ess_gene_list ESS_GENE_LIST
                        path to Essential gene list used for mask in Fc_ab and FcG_ab calcs
  -NEg NONESS_GENE_LIST, --NonEss_gene_list NONESS_GENE_LIST
                        path to Non-Essential gene list used for mask in Fc_ab and FcG_ab calcs
  -l LOG_FILE, --log_file LOG_FILE
                        Path to logger file
  -c CONTACT_FILES, --contact_files CONTACT_FILES
                        path to native contact files
  -r RESFEAT_FILES, --resFeat_files RESFEAT_FILES
                        path to residue Feature files
  -e UENT_FILES, --uent_files UENT_FILES
                        path to unique entanglement files
  -o OUTPATH, --outpath OUTPATH
                        path to output directory. will be made if doesnt exist
  -p NUM_PERMUTE, --num_permute NUM_PERMUTE
                        Number of permutations
  --floor FLOOR         Floor of Fc_ab contact matrix
```

If you have the [SLUG] then you can use the command files located [here](src/command_lists/compare_energetics_with_permutation.cmds) to reproduce loop forming contact enrichement analysis used in this work in the experimental data set and the AlphaFold structures. Please modify any other pathing as necessary. 
  

### Usage of [Plot_energetics_results_with_permutation.py](src/data/Plot_energetics_results_with_permutation.py)
```
usage: Plot_energetics_results_with_permutation.py [-h] -En ENORM -EsdE ESS_DELTAE -EsdEp ESS_DELTAE_PVALUES -NEsdE NONESS_DELTAE -NEsdEp NONESS_DELTAE_PVALUES -dDE
                                                   DELTADELTAE -dDEp DELTADELTAE_PVALUES -o OUTPATH -t TAG

Process regression data and generate plots.

options:
  -h, --help            show this help message and exit
  -En ENORM, --Enorm ENORM
                        path to Enorm file to plot
  -EsdE ESS_DELTAE, --Ess_deltaE ESS_DELTAE
                        path to Essential gene GT deltaE file to plot
  -EsdEp ESS_DELTAE_PVALUES, --Ess_deltaE_pvalues ESS_DELTAE_PVALUES
                        path to fdr pvalues file
  -NEsdE NONESS_DELTAE, --NonEss_deltaE NONESS_DELTAE
                        path to NonEssential gene GT deltaE file to plot
  -NEsdEp NONESS_DELTAE_PVALUES, --NonEss_deltaE_pvalues NONESS_DELTAE_PVALUES
                        path to fdr pvalues file
  -dDE DELTADELTAE, --deltaDeltaE DELTADELTAE
                        path to deltaDeltaE file to plot
  -dDEp DELTADELTAE_PVALUES, --deltaDeltaE_pvalues DELTADELTAE_PVALUES
                        path to fdr pvalues file
  -o OUTPATH, --outpath OUTPATH
                        Path to output directory.
  -t TAG, --tag TAG     tag for outfile: deltaE, Enorm,e ct....
```

If you have the [SLUG] then you can use the command files located [here](src/command_lists/Plot_compare_energetics_with_permutation.cmds) to reproduce the results plots for the loop forming contact enrichment analysis used in this work in the experimental data set and the AlphaFold structures. Please modify any other pathing as necessary. 
  
  
## Results of loop forming contact enrichment analysis
### Experimental structures
For the cyto-serum only set of observal proteins using 100k permutations  
![Cyto-serum loop forming contact enrichment](Figures/Contact_Enrichment/EXP/deltaDeltaE_C_50_p100000_EXP.png)  
Raw data for plotting can be found [here](data/Contact_Enrichment/EXP/C_50_p100000/)  

For the cyto-serum + Dnak only set of observal proteins using 100k permutations  
![Cyto-serum + DnaK loop forming contact enrichment](Figures/Contact_Enrichment/EXP/deltaDeltaE_CD_50_p100000_EXP.png) 
Raw data for plotting can be found [here](data/Contact_Enrichment/EXP/CD_50_p100000/)  
  
For the cyto-serum + GroEL only set of observal proteins using 100k permutations  
![Cyto-serum + GroEL loop forming contact enrichment](Figures/Contact_Enrichment/EXP/deltaDeltaE_CG_50_p100000_EXP.png)  
Raw data for plotting can be found [here](data/Contact_Enrichment/EXP/CG_50_p100000/)  
  
  

### Alphafold structures
For the cyto-serum only set of observal proteins using 100k permutations  
![Cyto-serum loop forming contact enrichment](Figures/Contact_Enrichment/AF/deltaDeltaE_C_50_p100000_AF.png)  
Raw data for plotting can be found [here](data/Contact_Enrichment/AF/C_50_p100000/)  
  
For the cyto-serum + Dnak only set of observal proteins using 100k permutations  
![Cyto-serum + DnaK loop forming contact enrichment](Figures/Contact_Enrichment/AF/deltaDeltaE_CD_50_p100000_AF.png)  
Raw data for plotting can be found [here](data/Contact_Enrichment/AF/CD_50_p100000/)  
  
For the cyto-serum + GroEL only set of observal proteins using 100k permutations  
![Cyto-serum + GroEL loop forming contact enrichment](Figures/Contact_Enrichment/AF/deltaDeltaE_CG_50_p100000_AF.png)  
Raw data for plotting can be found [here](data/Contact_Enrichment/AF/CG_50_p100000/)  
  


## Results of loop forming contact enrichment analysis (OR instead of contact potential)
### Experimental structures
For the cyto-serum only set of observal proteins using 100k permutations  
![Cyto-serum loop forming contact enrichment](Figures/Contact_Enrichment_OR_anal/EXP/OR_C_50_p100000_EXP.png)  
Raw data for plotting can be found [here](data/Contact_Enrichment_OR_anal/EXP/C_50_p100000/)  

For the cyto-serum + Dnak only set of observal proteins using 100k permutations  
![Cyto-serum + DnaK loop forming contact enrichment](Figures/Contact_Enrichment_OR_anal/EXP/OR_CD_50_p100000_EXP.png) 
Raw data for plotting can be found [here](data/Contact_Enrichment_OR_anal/EXP/CD_50_p100000/)  
  
For the cyto-serum + GroEL only set of observal proteins using 100k permutations  
![Cyto-serum + GroEL loop forming contact enrichment](Figures/Contact_Enrichment_OR_anal/EXP/OR_CG_50_p100000_EXP.png)  
Raw data for plotting can be found [here](data/Contact_Enrichment/EXP/CG_50_p100000/)  
  
  
### Alphafold structures
For the cyto-serum only set of observal proteins using 100k permutations  
![Cyto-serum loop forming contact enrichment](Figures/Contact_Enrichment_OR_anal/AF/OR_C_50_p100000_AF.png)  
Raw data for plotting can be found [here](data/Contact_Enrichment_OR_anal/AF/C_50_p100000/)  
  
For the cyto-serum + Dnak only set of observal proteins using 100k permutations  
![Cyto-serum + DnaK loop forming contact enrichment](Figures/Contact_Enrichment_OR_anal/AF/OR_CD_50_p100000_AF.png)  
Raw data for plotting can be found [here](data/Contact_Enrichment_OR_anal/AF/CD_50_p100000/)  
  
For the cyto-serum + GroEL only set of observal proteins using 100k permutations  
![Cyto-serum + GroEL loop forming contact enrichment](Figures/Contact_Enrichment_OR_anal/AF/OR_CG_50_p100000_AF.png)  
Raw data for plotting can be found [here](data/Contact_Enrichment_OR_anal/AF/CG_50_p100000/)  
  
The command files used to run these analysis using the OR instead of the "energetics" can be found [here](src/command_lists/Compare_OR_with_permutation.cmds) and to plot [here](src/command_lists/Plot_compare_OR_with_permutation.cmds). 


## Trend analysis in loop forming contacts 
From the loop forming contact enrichement analysis we observe that F-F, F-Y, T-Y, S-Y, and H-C loop forming contacts are enriched in non-essential proteins across all three LiPMS experimental conditions. 

### Usage of [Gen_SigContact_df.py](src/data/Gen_SigContact_df.py)
This script creates a dataframe with counts of all the loop forming contacts that match F-F, F-Y, T-Y, S-Y, and H-C for ech gene. 

```
usage: Gen_SigContact_df.py [-h] -l LOG_FILE -c CONTACT_FILES -e UENT_FILES -o OUTPATH

Process user specified arguments

options:
  -h, --help            show this help message and exit
  -l LOG_FILE, --log_file LOG_FILE
                        Path to logging file
  -c CONTACT_FILES, --contact_files CONTACT_FILES
                        path to native contact files
  -e UENT_FILES, --uent_files UENT_FILES
                        path to unique entanglement files
  -o OUTPATH, --outpath OUTPATH
                        path to output directory. will be made if doesnt exist
```
If you have the [SLUG] then you can use the command files located [here](/src/command_lists/Gen_SigContact_df.cmds) to reproduce the per gene significant contact count dataframe used in this work in the experimental data set and the AlphaFold structures. Please modify any other pathing as necessary. 


### Usage of [Trend_analysis.py](src/data/Trend_analysis.py)

```
usage: Trend_analysis.py [-h] -g GENE_LIST -Eg ESS_GENE_LIST -NEg NONESS_GENE_LIST -r RESFEAT_FILES -c LOOPCONTACT_DF -o OUTPATH -t TAG -b BUFF -s SPA --LiPMScov LIPMSCOV -l LOG_FILE

Process user specified arguments

options:
  -h, --help            show this help message and exit
  -g GENE_LIST, --gene_list GENE_LIST
                        path to all ent gene list use
  -Eg ESS_GENE_LIST, --ess_gene_list ESS_GENE_LIST
                        path to essentail ent gene list use
  -NEg NONESS_GENE_LIST, --noness_gene_list NONESS_GENE_LIST
                        path to nonessentail ent gene list use
  -r RESFEAT_FILES, --resFeat_files RESFEAT_FILES
                        path to residue Feature files
  -c LOOPCONTACT_DF, --loopcontact_df LOOPCONTACT_DF
                        path to dataframe containing loop contact classifications
  -o OUTPATH, --outpath OUTPATH
                        path to output directory. will be made if doesnt exist
  -t TAG, --tag TAG     tag for final output image
  -b BUFF, --buff BUFF  buffer used C CD CG
  -s SPA, --spa SPA     spa used 0, 10, 20, 30, 40, 50, 60, 70, 80, 90
  --LiPMScov LIPMSCOV   LiPMS coverage used 0, 10, 20, 30, 40, 50, 60, 70, 80, 90
  -l LOG_FILE, --log_file LOG_FILE
                        Path to logging file
```

If you have the [SLUG] then you can use the command files located [here](src/command_lists/Trend_analysis.cmds) to reproduce the trend analysis for the significant loop forming contacts found in this work in the experimental data set and the AlphaFold structures. Please modify any other pathing as necessary. 

## Results of Trend analysis in number of significant loop forming contacts
### Experimental structures
Fraction of proteins with atleast 1 of a loop closing contact type (cyto-serum only) 
![Fraction of proteins with atleast 1 of a loop closing contact type](Figures/Trend_analysis/FractionGenes_w_LoopContacts/EXP/EssVSNonEss_stats_atleast_n1_LFC_EXP_C_spa50_C_spa50_LiPMScov50.png)  
Raw data for this plot can be found [here](data/Trend_analysis/FractionGenes_w_LoopContacts/EXP/EssVSNonEss_stats_atleast_n1_LFC_EXP_C_spa50_C_spa50_LiPMScov50.csv).  
Results for other buffer systems can be found here [+Dnak](Figures/Trend_analysis/FractionGenes_w_LoopContacts/EXP/EssVSNonEss_stats_atleast_n1_LFC_EXP_CD_spa50_CD_spa50_LiPMScov50.png) [+GroEL](Figures/Trend_analysis/FractionGenes_w_LoopContacts/EXP/EssVSNonEss_stats_atleast_n1_LFC_EXP_CG_spa50_CG_spa50_LiPMScov50.png)  

<p float="left">
  <img src="Figures/Trend_analysis/OR_trends/EXP/atleast_n_loopcontacts_ORdiff_linear_reg_FF-FY-SY-HC_EXP_C_spa50_C_spa50_LiPMScov50.png" width="300" />
  <img src="Figures/Trend_analysis/OR_trends/EXP/atleast_n_loopcontacts_ORdiff_linear_reg_FF-FY-SY-HC_EXP_CD_spa50_CD_spa50_LiPMScov50.png" width="300" />
  <img src="Figures/Trend_analysis/OR_trends/EXP/atleast_n_loopcontacts_ORdiff_linear_reg_FF-FY-SY-HC_EXP_CG_spa50_CG_spa50_LiPMScov50.png" width="300" />
</p>

<p float="left">
  <img src="Figures/Trend_analysis/OR_trends/EXP/n_loopcontacts_ORdiff_linear_reg_FF-FY-SY-HC_EXP_C_spa50_C_spa50_LiPMScov50.png" width="300" />
  <img src="Figures/Trend_analysis/OR_trends/EXP/n_loopcontacts_ORdiff_linear_reg_FF-FY-SY-HC_EXP_CD_spa50_CD_spa50_LiPMScov50.png" width="300" />
  <img src="Figures/Trend_analysis/OR_trends/EXP/n_loopcontacts_ORdiff_linear_reg_FF-FY-SY-HC_EXP_CG_spa50_CG_spa50_LiPMScov50.png" width="300" />
</p>  

Raw data for these plots can be found [here](data/Trend_analysis/OR_trends/EXP/).  

### Alphafold structures
Fraction of proteins with atleast 1 of a loop closing contact type (cyto-serum only)   
![Fraction of proteins with atleast 1 of a loop closing contact type](Figures/Trend_analysis/FractionGenes_w_LoopContacts/AF/EssVSNonEss_stats_atleast_n1_LFC_AF_C_spa50_C_spa50_LiPMScov50.png)  
Results for other buffer systems can be found here [+Dnak](Figures/Trend_analysis/FractionGenes_w_LoopContacts/AF/EssVSNonEss_stats_atleast_n1_LFC_AF_CD_spa50_CD_spa50_LiPMScov50.png) [+GroEL](Figures/Trend_analysis/FractionGenes_w_LoopContacts/AF/EssVSNonEss_stats_atleast_n1_LFC_AF_CG_spa50_CG_spa50_LiPMScov50.png)  

<p float="left">
  <img src="Figures/Trend_analysis/OR_trends/AF/atleast_n_loopcontacts_ORdiff_linear_reg_FF-FY-SY-HC_AF_C_spa50_C_spa50_LiPMScov50.png" width="300" />
  <img src="Figures/Trend_analysis/OR_trends/AF/atleast_n_loopcontacts_ORdiff_linear_reg_FF-FY-SY-HC_AF_CD_spa50_CD_spa50_LiPMScov50.png" width="300" />
  <img src="Figures/Trend_analysis/OR_trends/AF/atleast_n_loopcontacts_ORdiff_linear_reg_FF-FY-SY-HC_AF_CG_spa50_CG_spa50_LiPMScov50.png" width="300" />
</p>

<p float="left">
  <img src="Figures/Trend_analysis/OR_trends/AF/n_loopcontacts_ORdiff_linear_reg_FF-FY-SY-HC_AF_C_spa50_C_spa50_LiPMScov50.png" width="300" />
  <img src="Figures/Trend_analysis/OR_trends/AF/n_loopcontacts_ORdiff_linear_reg_FF-FY-SY-HC_AF_CD_spa50_CD_spa50_LiPMScov50.png" width="300" />
  <img src="Figures/Trend_analysis/OR_trends/AF/n_loopcontacts_ORdiff_linear_reg_FF-FY-SY-HC_AF_CG_spa50_CG_spa50_LiPMScov50.png" width="300" />
</p>  

Raw data for these plots can be found [here](data/Trend_analysis/OR_trends/AF/).  

## Results for Fraction of amino acid types in datasets
### Experimental structures  
_Fraction of General hydrophobic residues per protein in the cyto-serum only buffer at the 50th SPA and LiPMScov thresholds_
![Fraction of General hydrophobic residues per protein](Figures/CalcHydrophobicFraction/EXP/EXP_C_Rall_spa50_LiPMScov50_perProt_Fraction_General_hydrophobic.png)   
The raw stats can be found [here](data/CalcHydrophobicFraction/EXP/EXP_C_Rall_spa50_LiPMScov50_perProt_Fraction_General_hydrophobic_STATS.csv)  

_Fraction of All hydrophobic residues per protein in the cyto-serum only buffer at the 50th SPA and LiPMScov thresholds_
![Fraction of Total hydrophobic residues per protein](Figures/CalcHydrophobicFraction/EXP/EXP_C_Rall_spa50_LiPMScov50_perProt_Fraction_Total_hydrophobic.png)   
The raw stats can be found [here](data/CalcHydrophobicFraction/EXP/EXP_C_Rall_spa50_LiPMScov50_perProt_Fraction_Total_hydrophobic_STATS.csv)  
  
_Fraction of strong hydrophobic residues per protein in the cyto-serum only buffer at the 50th SPA and LiPMScov thresholds_
![Fraction of strong hydrophobic residues per protein](Figures/CalcHydrophobicFraction/EXP/EXP_C_Rall_spa50_LiPMScov50_perProt_Fraction_Strong_hydrophobic.png)   
The raw stats can be found [here](data/CalcHydrophobicFraction/EXP/EXP_C_Rall_spa50_LiPMScov50_perProt_Fraction_Strong_hydrophobic_STATS.csv)  

_Fraction of Weak hydrophobic residues per protein in the cyto-serum only buffer at the 50th SPA and LiPMScov thresholds_
![Fraction of Weak hydrophobic residues per protein](Figures/CalcHydrophobicFraction/EXP/EXP_C_Rall_spa50_LiPMScov50_perProt_Fraction_Weak_hydrophobic.png)   
The raw stats can be found [here](data/CalcHydrophobicFraction/EXP/EXP_C_Rall_spa50_LiPMScov50_perProt_Fraction_Weak_hydrophobic_STATS.csv) 

_Fraction of Total hydrophobic residues per protein in whole experimental dataset_ 
![Fraction of Total hydrophobic residues per protein](Figures/CalcHydrophobicFraction/EXP/EXP_Total_Rall_spaTotal_LiPMScovTotal_perProt_Fraction_Total_hydrophobic.png)  
The raw stats can be found [here](data/CalcHydrophobicFraction/EXP/EXP_Total_Rall_spaTotal_LiPMScovTotal_perProt_Fraction_Total_hydrophobic_STATS.csv)

_Fraction of strong hydrophobic residues per protein in whole experimental dataset_ 
![Fraction of strong hydrophobic residues per protein](Figures/CalcHydrophobicFraction/EXP/EXP_Total_Rall_spaTotal_LiPMScovTotal_perProt_Fraction_Strong_hydrophobic.png)  
The raw stats can be found [here](data/CalcHydrophobicFraction/EXP/EXP_Total_Rall_spaTotal_LiPMScovTotal_perProt_Fraction_Strong_hydrophobic_STATS.csv)

_Fraction of Weak hydrophobic residues per protein in whole experimental dataset_ 
![Fraction of Weak hydrophobic residues per protein](Figures/CalcHydrophobicFraction/EXP/EXP_Total_Rall_spaTotal_LiPMScovTotal_perProt_Fraction_Weak_hydrophobic.png)  
The raw stats can be found [here](data/CalcHydrophobicFraction/EXP/EXP_Total_Rall_spaTotal_LiPMScovTotal_perProt_Fraction_Weak_hydrophobic_STATS.csv)

### Alphafold structures
_Fraction of General hydrophobic residues per protein in the cyto-serum only buffer at the 50th SPA and LiPMScov thresholds_
![Fraction of General hydrophobic residues per protein](Figures/CalcHydrophobicFraction/AF/AF_C_Rall_spa50_LiPMScov50_perProt_Fraction_General_hydrophobic.png)   
The raw stats can be found [here](data/CalcHydrophobicFraction/AF/AF_C_Rall_spa50_LiPMScov50_perProt_Fraction_General_hydrophobic_STATS.csv)  

_Fraction of All hydrophobic residues per protein in the cyto-serum only buffer at the 50th SPA and LiPMScov thresholds_
![Fraction of Total hydrophobic residues per protein](Figures/CalcHydrophobicFraction/AF/AF_C_Rall_spa50_LiPMScov50_perProt_Fraction_Total_hydrophobic.png)   
The raw stats can be found [here](data/CalcHydrophobicFraction/AF/AF_C_Rall_spa50_LiPMScov50_perProt_Fraction_Total_hydrophobic_STATS.csv)  
  
_Fraction of strong hydrophobic residues per protein in the cyto-serum only buffer at the 50th SPA and LiPMScov thresholds_
![Fraction of strong hydrophobic residues per protein](Figures/CalcHydrophobicFraction/AF/AF_C_Rall_spa50_LiPMScov50_perProt_Fraction_Strong_hydrophobic.png)   
The raw stats can be found [here](data/CalcHydrophobicFraction/AF/AF_C_Rall_spa50_LiPMScov50_perProt_Fraction_Strong_hydrophobic_STATS.csv)  

_Fraction of Weak hydrophobic residues per protein in the cyto-serum only buffer at the 50th SPA and LiPMScov thresholds_
![Fraction of Weak hydrophobic residues per protein](Figures/CalcHydrophobicFraction/AF/AF_C_Rall_spa50_LiPMScov50_perProt_Fraction_Weak_hydrophobic.png)   
The raw stats can be found [here](data/CalcHydrophobicFraction/AF/AF_C_Rall_spa50_LiPMScov50_perProt_Fraction_Weak_hydrophobic_STATS.csv) 

_Fraction of Total hydrophobic residues per protein in whole AF dataset_ 
![Fraction of Total hydrophobic residues per protein](Figures/CalcHydrophobicFraction/AF/AF_Total_Rall_spaTotal_LiPMScovTotal_perProt_Fraction_Total_hydrophobic.png)  
The raw stats can be found [here](data/CalcHydrophobicFraction/AF/AF_Total_Rall_spaTotal_LiPMScovTotal_perProt_Fraction_Total_hydrophobic_STATS.csv)

_Fraction of strong hydrophobic residues per protein in whole AF dataset_ 
![Fraction of strong hydrophobic residues per protein](Figures/CalcHydrophobicFraction/AF/AF_Total_Rall_spaTotal_LiPMScovTotal_perProt_Fraction_Strong_hydrophobic.png)  
The raw stats can be found [here](data/CalcHydrophobicFraction/AF/AF_Total_Rall_spaTotal_LiPMScovTotal_perProt_Fraction_Strong_hydrophobic_STATS.csv)

_Fraction of Weak hydrophobic residues per protein in whole AF dataset_ 
![Fraction of Weak hydrophobic residues per protein](Figures/CalcHydrophobicFraction/AF/AF_Total_Rall_spaTotal_LiPMScovTotal_perProt_Fraction_Weak_hydrophobic.png)  
The raw stats can be found [here](data/CalcHydrophobicFraction/AF/AF_Total_Rall_spaTotal_LiPMScovTotal_perProt_Fraction_Weak_hydrophobic_STATS.csv)
