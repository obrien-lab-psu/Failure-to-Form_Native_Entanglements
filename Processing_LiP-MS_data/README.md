# Processing_LiP-MS_data  

## General workflow
```mermaid
graph TD
    A[FLiPPR] --> B[FDR] --> C[Thresholding]
    D[Timepoint Overlap]
    E[Sel. Refolded]
    click A "https://github.com/obrien-lab-psu/Failure-to-Form_Native_Entanglements/tree/main/Processing_LiP-MS_data#flippr"
    click B "https://github.com/obrien-lab-psu/Failure-to-Form_Native_Entanglements/tree/main/Processing_LiP-MS_data#false-discovery-rate-correction"
    click C "https://github.com/obrien-lab-psu/Failure-to-Form_Native_Entanglements/tree/main/Processing_LiP-MS_data#thresholding-observed-proteins-by-spa-and-lipms-coverage"
    click D "https://github.com/obrien-lab-psu/Failure-to-Form_Native_Entanglements/tree/main/Processing_LiP-MS_data#overlap-of-misfolded-genes-and-pk-sites-across-1min-5min-and-2hr-of-refolding-in-lip-ms-experiments"
    click E "https://github.com/obrien-lab-psu/Failure-to-Form_Native_Entanglements/tree/main/Processing_LiP-MS_data#selecting-refolded-candidates-for-simulation"

``` 

## FLiPPR
We analyze the raw proteome discover data using the [FLiPPR](https://pubs.acs.org/doi/full/10.1021/acs.jproteome.3c00887)[[1]](#1) protocol to ensure peptides are merged not only on the identity level but also the PK cut-site level.  

### A high level overview of the theory
Ions are initially separated into 4 classes:  
1. Invalid peptide - not enough abundance in one or both cases (will be ignored in downstream analysis)
2. No missing abundances across all 6 samples  
3. One abundance missing from either the Native or Refolded cases (missing value is dropped and statistical tests are preformed on the smaller sample set) 
4. All-or-None: One sample is missing all 3 abundances while the other has all 3 (The 3 missing values are filled by Gaussian imputation of a distribution with mean and variance meant to mimic the lower limit of detection of the Mass Spec).  

For all the ions that can inform on the susceptibility at a given cut-site, the (ratio, P-value) pairs for those ions are collectively considered. If they all agree in direction (i.e., the signs of the t test statistics are all the same), then the oerall ratio for the cut-site is calculated by taking the median of the ratios of all ions that map to it, and the P-values are combined with Fishers method to provide an updated (ratio, P-value) for the cut-site.

If there are two independent ions and they disagree (e.g., the ion is more abundant in the test condition in the 2+ charge state but more abundant in the control condition in the 3+ charge state), then a median is still calculated, but the P-value is set to 1, implying there is no confidence as to whether this cut-site was more susceptible in the test or control condition. These cut-sites are discounted from the tally of the total valid cut-sites.

If there are three ions, then a majority rules heuristic is applied: the disagreeing ion is disregarded, and the (ratio, P-value)s are only combined for the agreeing ions. In practice, it is relatively rare for more than three ions to be mapped to the same cut-site, but where this occurs, if a majority (or all) of the ions agree in direction, they are combined. If there is a tie, the P-value is set to 1.

The false discovery rate correction is handled differently in this [work](https://github.com/obrien-lab-psu/Failure-to-Form_Native_Entanglements/tree/main/Processing_LiP-MS_data#usage-of-fdr_correction_of_lipms_datapy).  

### Usage of [FLiPPR.py](src/data/FLiPPR.py)
```
usage: FLiPPR.py [-h] -i INPFILE -o OUTPATH

Process user specified arguments

options:
  -h, --help            show this help message and exit
  -i INPFILE, --inpfile INPFILE
                        Path to raw LiPMS file to process
  -o OUTPATH, --outpath OUTPATH
                        Path to output directory
```

If you have the [SLUG] then you can use the command files located [here](src/command_lists/FliPPR.cmds) to reproduce the FLiPPR analysis done in this work. Please modify any other pathing as necessary. 

## False discovery rate correction
We perform the benjamini hochberg FDR across the whole proteome wide set of half-tryptic peptides instead of the individual proteins.   

### Usage of [FDR_correction_of_LiPMS_data.py](src/data/FDR_correction_of_LiPMS_data.py)
```
usage: FDR_correction_of_LiPMS_data.py [-h] -f LIPMS_FILES -t THRESHOLD -o OUTPATH -a ALPHA

Process user specified arguments for LiPMS files correction.

options:
  -h, --help            show this help message and exit
  -f LIPMS_FILES, --lipms_files LIPMS_FILES
                        path to LiPMS files to correct
  -t THRESHOLD, --threshold THRESHOLD
                        abundance change threshold
  -o OUTPATH, --outpath OUTPATH
                        outpath
  -a ALPHA, --alpha ALPHA
                        family wise error rate: alpha
```

If you have the [SLUG] then you can use the command files located [here](src/command_lists/FDR.cmds) to reproduce the FDR corrected FLiPPR analysis done in this work. Please modify any other pathing as necessary. 


## Thresholding observed proteins by SPA and LiPMS-Coverage
For robustness of results we scan all analysis done in this work along the sum of the peptide abundance SPA percentiles. We also only consider those where atleast 50% of the uniprot canonical sequence was identified in the experiment.   

### Usage of [thresholding.py](src/data/thresholding.py)
```
usage: thresholding.py [-h] -i INPFILES -o OUTPATH -t TAG

Process user specified arguments

options:
  -h, --help            show this help message and exit
  -i INPFILES, --inpfiles INPFILES
                        Path to raw LiPMS file to process
  -o OUTPATH, --outpath OUTPATH
                        Path to output directory
  -t TAG, --tag TAG     Thresholding tag: spa or cov
```

If you have the [SLUG] then you can use the command files located [here](src/command_lists/thresholding.cmds) to reproduce the thresholding of genes along the SPA and or LiPMS-COV analysis done in this work. Please modify any other pathing as necessary. 

## Overlap of misfolded genes and PK sites across 1min, 5min, and 2hr of refolding in LiP-MS experiments
We pool the observed misfolding events across all refolding times in this work for two reasons: 
1. we only care about whether a misfolding event was observed and not necessiarly about its persistence on long time scales.
2. due to the spare nature of the signals resulting from the LiP-MS experiment we pool the data to increase statistical power. 

We therefore need to show that there is significant overlap between the timepoints in both genes observed to misfold and those misfolding sites.  

### Usage of [Overlap.py](src/data/Overlap.py)
```
usage: Overlap.py [-h] -i INPFILES -o OUTPATH -t TAG -g GENE_LIST

Process user specified arguments

options:
  -h, --help            show this help message and exit
  -i INPFILES, --inpfiles INPFILES
                        Path to FLiPPR processed files
  -o OUTPATH, --outpath OUTPATH
                        Path to output directory
  -t TAG, --tag TAG     Thresholding tag: spa or cov
  -g GENE_LIST, --gene_list GENE_LIST
                        gene list to mask with
```

If you have the [SLUG] then you can use the command files located [here](src/command_lists/Overlap.cmds) to reproduce the analysis of overlap of genes and PK-cutsites at various refolding times done in this work. Please modify any other pathing as necessary.  

### Results of overlap analysis
Here we analyze the overlap of those genes that had a native SPA greater thant the 50th percentile and atleast 50% of their canonical uniprot sequence resolved in the experiment in the absence of any chaperones. We observe a stastically signifcant overlap of 67.9% across all three refolding times (with 40.1% of the signficant changes in structure consistent) and 76.6% when considering the 5min and 2hr timepoints (with 57.7% of the signficant changes in structure consistent).   

#### Gene overlap (Experimental)
![EXP GeneOverlapVenn3Diagram_spa50_LiPMScov50_ent_genes_C](Figures/GeneOverlap/EXP/GeneOverlapVenn3Diagram_spa50_LiPMScov50_ent_genes_C.png)
![EXP GeneOverlapVenn2Diagram_spa50_LiPMScov50_ent_genes_C](Figures/GeneOverlap/EXP/GeneOverlapVenn2Diagram_spa50_LiPMScov50_ent_genes_C.png) 
#### Gene overlap (Alphafold) 
![AF GeneOverlapVenn3Diagram_spa50_LiPMScov50_ent_genes_C](Figures/GeneOverlap/AF/GeneOverlapVenn3Diagram_spa50_LiPMScov50_ent_genes_C.png)
![AF GeneOverlapVenn2Diagram_spa50_LiPMScov50_ent_genes_C](Figures/GeneOverlap/AF/GeneOverlapVenn2Diagram_spa50_LiPMScov50_ent_genes_C.png) 


#### PK site overlap (Experimental)
![EXP PKsiteOverlapVenn3Diagram_spa50_LiPMScov50_ent_genes_C](Figures/PKsiteOverlap/EXP/PKsiteOverlapVenn3Diagram_spa50_LiPMScov50_ent_genes_C.png)
![EXP PKsiteOverlapVenn2Diagram_spa50_LiPMScov50_ent_genes_C](Figures/PKsiteOverlap/EXP/PKsiteOverlapVenn2Diagram_spa50_LiPMScov50_ent_genes_C.png) 
#### PK site overlap (Alphafold)
![AF PKsiteOverlapVenn3Diagram_spa50_LiPMScov50_ent_genes_C](Figures/PKsiteOverlap/AF/PKsiteOverlapVenn3Diagram_spa50_LiPMScov50_ent_genes_C.png)
![AF PKsiteOverlapVenn2Diagram_spa50_LiPMScov50_ent_genes_C](Figures/PKsiteOverlap/AF/PKsiteOverlapVenn2Diagram_spa50_LiPMScov50_ent_genes_C.png) 

## Selecting refolded candidates for simulation
Refolded candidates for our simulations meet the following criteria:
1. have a SPA >= 50th percentile value  
2. have coverage in the LiP-MS experiments >= to 50% of the canonical uniprot sequence  
3. have no PK cut-sites with a significant change in abundance between treated and untreated samples in the cyto-serum buffer  
4. are refolded across all three refolding schedules (1min, 5min 2hrs)  
  
### Usage of [Refolded.py](src/data/Refolded.py)
```
usage: Refolded.py [-h] -f LIPMS_FILES -t THRESHOLD -o OUTPATH -a ALPHA -e ENT_GENES

Process user specified arguments for LiPMS files correction.

options:
  -h, --help            show this help message and exit
  -f LIPMS_FILES, --lipms_files LIPMS_FILES
                        path to LiPMS files to correct
  -t THRESHOLD, --threshold THRESHOLD
                        abundance change threshold
  -o OUTPATH, --outpath OUTPATH
                        outpath
  -a ALPHA, --alpha ALPHA
                        family wise error rate: alpha
  -e ENT_GENES, --ent_genes ENT_GENES
                        path to list of ent genes to mask over the set of refolded to get the final set of genes
```

If you have the [SLUG] then you can use the command files located [here](src/command_lists/Refolded.cmds) to reproduce the sets of refolded genes done in this work. Please modify any other pathing as necessary. The resulting lists of refolded genes can be found [here](data/Refolded/EXP_all/ALL_Refolded.csv) and for those with native entanglements [here](data/Refolded/EXP_all/ENT_Refolded.csv).  

### References  

<a id="1">[1]</a>: Manriquez-Sandoval, E., Brewer, J., Lule, G., Lopez, S., & Fried, S. D. (2024). FLiPPR: A Processor for Limited Proteolysis (LiP) Mass Spectrometry Data Sets Built on FragPipe. Journal of Proteome Research, 23(7), 2332–2342. https://doi.org/10.1021/acs.jproteome.3c00887 