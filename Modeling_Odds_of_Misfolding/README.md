# Modeling_Odds_of_Misfolding  
Here we model the log-odds of observing a change in proteolysis suseptibility (i.e. misfolding) as a function of the amino acid composition and the region of the protein it was observed in (i.e. entangled region or non-entangled region). We use logistic regression and fit the following model using the maximal likelihood optimiztion method.

$$
\text{log(odds)} \approx B_0 + B_1 \cdot X_1 + B_2 \cdot X_2
$$

## General workflow
```mermaid
graph TD
    A[Fit Model] --> B[Propensity Score Matching]
``` 

## Binomial logistic regression



### Usage of [src/data/Regression.py](src/data/Regression.py)
```
usage: Regression.py [-h] -f RESFEAT_FILES -o OUTPATH -g GENE_LIST -t TAG -l LOAD_STYLE -b BUFFER -s SPA -c COV -r REG_FORMULA -v HOLD_VAR

Process user specified arguments

options:
  -h, --help            show this help message and exit
  -f RESFEAT_FILES, --resFeat_files RESFEAT_FILES
                        Path to residue feature files
  -o OUTPATH, --outpath OUTPATH
                        Path to output directory
  -g GENE_LIST, --gene_list GENE_LIST
                        Path to gene list to use
  -t TAG, --tag TAG     Tag for output filenames
  -l LOAD_STYLE, --load_style LOAD_STYLE
                        Load style (True: load by gene, False: load a single file with all genes present)
  -b BUFFER, --buffer BUFFER
                        Buffer system to use
  -s SPA, --spa SPA     SPA threshold
  -c COV, --cov COV     LiPMS cov threshold
  -r REG_FORMULA, --reg_formula REG_FORMULA
                        Regression formula
  -v HOLD_VAR, --hold_var HOLD_VAR
                        Variable to hold constant while calculating odds
```


If you have the [SLUG] then you can use the command files located [here](src/command_lists/Regressions_whole_proteome.cmds) to reproduce the binomial logistic regressions used in this work to fit the experimental data and the AlphaFold structures. Please modify any other pathing as necessary. 

### Results 
When considering only those proteins with native entanglements we observe a statistically significant greater odds of misfolding in the entangled region of proteins rather than not. Even in the presence of DnaK and GroEL.  
![All proteins with native entanglements](Figures/Regressions/whole_proteome/EXP/ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.png)    
  
When we restrict ourselves to the set of essential proteins determined by knockout experiments reported in the [DEG](http://origin.tubic.org/deg/public/index.php) database we observe no statistical association in the presence of chaperones.  
![Essential proteins with native entanglements](Figures/Regressions/whole_proteome/EXP/essential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.png)  
  
While for the complimentary set of non-essential proteins we still observe the bias towards misfolding involving native entanglements in the presence of chaperones.  
![Non-Essential proteins with native entanglements](Figures/Regressions/whole_proteome/EXP/nonessential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.png)  
  
  
Below is a table of all our results across both the sets of experimentally derived structures and Alphafold structures as well as at each individual time point observed in the LiPMS experiments.  

| Set Type               | Timepoint | Structure Type | Entangled Only | Link to File                                                                                   |
|------------------------|-----------|----------------|-----------------|------------------------------------------------------------------------------------------------|
| all                    | R1min     | EXP            | No-ENT              | [all_genes_R1min](Figures/Regressions/whole_proteome/EXP/all_genes_R1min_binomial_regression_results_var-region_LiPMScov50.png)       |
| essential_ent           | R1min     | EXP            | Yes-ENT             | [essential_ent_genes_R1min](Figures/Regressions/whole_proteome/EXP/essential_ent_genes_R1min_binomial_regression_results_var-region_LiPMScov50.png) |
| nonessential_ent        | R1min     | EXP            | Yes-ENT             | [nonessential_ent_genes_R1min](Figures/Regressions/whole_proteome/EXP/nonessential_ent_genes_R1min_binomial_regression_results_var-region_LiPMScov50.png) |
| ent                    | R1min     | EXP            | Yes-ENT             | [ent_genes_R1min](Figures/Regressions/whole_proteome/EXP/ent_genes_R1min_binomial_regression_results_var-region_LiPMScov50.png)         |
| essential              | R1min     | EXP            | No-ENT              | [essential_genes_R1min](Figures/Regressions/whole_proteome/EXP/essential_genes_R1min_binomial_regression_results_var-region_LiPMScov50.png) |
| nonessential           | R1min     | EXP            | No-ENT              | [nonessential_genes_R1min](Figures/Regressions/whole_proteome/EXP/nonessential_genes_R1min_binomial_regression_results_var-region_LiPMScov50.png) |
| all                    | R5min     | EXP            | No-ENT              | [all_genes_R5min](Figures/Regressions/whole_proteome/EXP/all_genes_R5min_binomial_regression_results_var-region_LiPMScov50.png)     |
| essential_ent           | R5min     | EXP            | Yes-ENT             | [essential_ent_genes_R5min](Figures/Regressions/whole_proteome/EXP/essential_ent_genes_R5min_binomial_regression_results_var-region_LiPMScov50.png) |
| nonessential_ent        | R5min     | EXP            | Yes-ENT             | [nonessential_ent_genes_R5min](Figures/Regressions/whole_proteome/EXP/nonessential_ent_genes_R5min_binomial_regression_results_var-region_LiPMScov50.png) |
| ent                    | R5min     | EXP            | Yes-ENT             | [ent_genes_R5min](Figures/Regressions/whole_proteome/EXP/ent_genes_R5min_binomial_regression_results_var-region_LiPMScov50.png)         |
| essential              | R5min     | EXP            | No-ENT              | [essential_genes_R5min](Figures/Regressions/whole_proteome/EXP/essential_genes_R5min_binomial_regression_results_var-region_LiPMScov50.png) |
| nonessential           | R5min     | EXP            | No-ENT              | [nonessential_genes_R5min](Figures/Regressions/whole_proteome/EXP/nonessential_genes_R5min_binomial_regression_results_var-region_LiPMScov50.png) |
| all                    | R2hr      | EXP            | No-ENT              | [all_genes_R2hr](Figures/Regressions/whole_proteome/EXP/all_genes_R2hr_binomial_regression_results_var-region_LiPMScov50.png)       |
| essential_ent           | R2hr      | EXP            | Yes-ENT             | [essential_ent_genes_R2hr](Figures/Regressions/whole_proteome/EXP/essential_ent_genes_R2hr_binomial_regression_results_var-region_LiPMScov50.png) |
| nonessential_ent        | R2hr      | EXP            | Yes-ENT             | [nonessential_ent_genes_R2hr](Figures/Regressions/whole_proteome/EXP/nonessential_ent_genes_R2hr_binomial_regression_results_var-region_LiPMScov50.png) |
| ent                    | R2hr      | EXP            | Yes-ENT             | [ent_genes_R2hr](Figures/Regressions/whole_proteome/EXP/ent_genes_R2hr_binomial_regression_results_var-region_LiPMScov50.png)         |
| essential              | R2hr      | EXP            | No-ENT              | [essential_genes_R2hr](Figures/Regressions/whole_proteome/EXP/essential_genes_R2hr_binomial_regression_results_var-region_LiPMScov50.png) |
| nonessential           | R2hr      | EXP            | No-ENT              | [nonessential_genes_R2hr](Figures/Regressions/whole_proteome/EXP/nonessential_genes_R2hr_binomial_regression_results_var-region_LiPMScov50.png) |
| all                    | Rall      | EXP            | No-ENT              | [all_genes_Rall](Figures/Regressions/whole_proteome/EXP/all_genes_Rall_binomial_regression_results_var-region_LiPMScov50.png)        |
| essential_ent           | Rall      | EXP            | Yes-ENT             | [essential_ent_genes_Rall](Figures/Regressions/whole_proteome/EXP/essential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.png) |
| nonessential_ent        | Rall      | EXP            | Yes-ENT             | [nonessential_ent_genes_Rall](Figures/Regressions/whole_proteome/EXP/nonessential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.png) |
| ent                    | Rall      | EXP            | Yes-ENT             | [ent_genes_Rall](Figures/Regressions/whole_proteome/EXP/ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.png)         |
| essential              | Rall      | EXP            | No-ENT              | [essential_genes_Rall](Figures/Regressions/whole_proteome/EXP/essential_genes_Rall_binomial_regression_results_var-region_LiPMScov50.png) |
| nonessential           | Rall      | EXP            | No-ENT              | [nonessential_genes_Rall](Figures/Regressions/whole_proteome/EXP/nonessential_genes_Rall_binomial_regression_results_var-region_LiPMScov50.png) |
| all                    | R1min     | AF            | No-ENT              | [all_genes_R1min](Figures/Regressions/whole_proteome/AF/all_genes_R1min_binomial_regression_results_var-region_LiPMScov50.png)       |
| essential_ent           | R1min     | AF            | Yes-ENT             | [essential_ent_genes_R1min](Figures/Regressions/whole_proteome/AF/essential_ent_genes_R1min_binomial_regression_results_var-region_LiPMScov50.png) |
| nonessential_ent        | R1min     | AF            | Yes-ENT             | [nonessential_ent_genes_R1min](Figures/Regressions/whole_proteome/AF/nonessential_ent_genes_R1min_binomial_regression_results_var-region_LiPMScov50.png) |
| ent                    | R1min     | AF            | Yes-ENT             | [ent_genes_R1min](Figures/Regressions/whole_proteome/AF/ent_genes_R1min_binomial_regression_results_var-region_LiPMScov50.png)         |
| essential              | R1min     | AF            | No-ENT              | [essential_genes_R1min](Figures/Regressions/whole_proteome/AF/essential_genes_R1min_binomial_regression_results_var-region_LiPMScov50.png) |
| nonessential           | R1min     | AF            | No-ENT              | [nonessential_genes_R1min](Figures/Regressions/whole_proteome/AF/nonessential_genes_R1min_binomial_regression_results_var-region_LiPMScov50.png) |
| all                    | R5min     | AF            | No-ENT              | [all_genes_R5min](Figures/Regressions/whole_proteome/AF/all_genes_R5min_binomial_regression_results_var-region_LiPMScov50.png)     |
| essential_ent           | R5min     | AF            | Yes-ENT             | [essential_ent_genes_R5min](Figures/Regressions/whole_proteome/AF/essential_ent_genes_R5min_binomial_regression_results_var-region_LiPMScov50.png) |
| nonessential_ent        | R5min     | AF            | Yes-ENT             | [nonessential_ent_genes_R5min](Figures/Regressions/whole_proteome/AF/nonessential_ent_genes_R5min_binomial_regression_results_var-region_LiPMScov50.png) |
| ent                    | R5min     | AF            | Yes-ENT             | [ent_genes_R5min](Figures/Regressions/whole_proteome/AF/ent_genes_R5min_binomial_regression_results_var-region_LiPMScov50.png)         |
| essential              | R5min     | AF            | No-ENT              | [essential_genes_R5min](Figures/Regressions/whole_proteome/AF/essential_genes_R5min_binomial_regression_results_var-region_LiPMScov50.png) |
| nonessential           | R5min     | AF            | No-ENT              | [nonessential_genes_R5min](Figures/Regressions/whole_proteome/AF/nonessential_genes_R5min_binomial_regression_results_var-region_LiPMScov50.png) |
| all                    | R2hr      | AF            | No-ENT              | [all_genes_R2hr](Figures/Regressions/whole_proteome/AF/all_genes_R2hr_binomial_regression_results_var-region_LiPMScov50.png)       |
| essential_ent           | R2hr      | AF            | Yes-ENT             | [essential_ent_genes_R2hr](Figures/Regressions/whole_proteome/AF/essential_ent_genes_R2hr_binomial_regression_results_var-region_LiPMScov50.png) |
| nonessential_ent        | R2hr      | AF            | Yes-ENT             | [nonessential_ent_genes_R2hr](Figures/Regressions/whole_proteome/AF/nonessential_ent_genes_R2hr_binomial_regression_results_var-region_LiPMScov50.png) |
| ent                    | R2hr      | AF            | Yes-ENT             | [ent_genes_R2hr](Figures/Regressions/whole_proteome/AF/ent_genes_R2hr_binomial_regression_results_var-region_LiPMScov50.png)         |
| essential              | R2hr      | AF            | No-ENT              | [essential_genes_R2hr](Figures/Regressions/whole_proteome/AF/essential_genes_R2hr_binomial_regression_results_var-region_LiPMScov50.png) |
| nonessential           | R2hr      | AF            | No-ENT              | [nonessential_genes_R2hr](Figures/Regressions/whole_proteome/AF/nonessential_genes_R2hr_binomial_regression_results_var-region_LiPMScov50.png) |
| all                    | Rall      | AF            | No-ENT              | [all_genes_Rall](Figures/Regressions/whole_proteome/AF/all_genes_Rall_binomial_regression_results_var-region_LiPMScov50.png)        |
| essential_ent           | Rall      | AF            | Yes-ENT             | [essential_ent_genes_Rall](Figures/Regressions/whole_proteome/AF/essential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.png) |
| nonessential_ent        | Rall      | AF            | Yes-ENT             | [nonessential_ent_genes_Rall](Figures/Regressions/whole_proteome/AF/nonessential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.png) |
| ent                    | Rall      | AF            | Yes-ENT             | [ent_genes_Rall](Figures/Regressions/whole_proteome/AF/ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.png)         |
| essential              | Rall      | AF            | No-ENT              | [essential_genes_Rall](Figures/Regressions/whole_proteome/AF/essential_genes_Rall_binomial_regression_results_var-region_LiPMScov50.png) |
| nonessential           | Rall      | AF            | No-ENT              | [nonessential_genes_Rall](Figures/Regressions/whole_proteome/AF/nonessential_genes_Rall_binomial_regression_results_var-region_LiPMScov50.png) |


