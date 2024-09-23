# Chaperone_Client_Associations  
Here we asses the statistical association known DnaK and GroEL clients with gene essentiality and native entanglemnts. We collate together a set of client proteins from experiments that directly measure the association of a given protein with either DnaK or GroEL either via co-immunoprecipitaton or his-tag purification followed by proteomics methods for identification. 

## General workflow
```mermaid
graph TD
    A[Assoc. essentiality and client]
    B[DnaK binding site]
``` 

## Association of gene essentiality and clien status
We use the Fisher exact test to determine if there is any statistically significant associaiton between proteins being essential and a chaperone client. 
The contingency table is designed as follows:  
| | Client Yes | Client No 
|------|-----------|----------------|
| Essential & Ent. | | |
| Nonessential & Ent. | | |

### Usage of [Assoc_Client_n_Essential.py](src/data/Assoc_Client_n_Essential.py)
```
usage: Assoc_Client_n_Essential.py [-h] -Ag GENE_LIST -Eg ESS_GENE_LIST -Entg ENT_GENE_LIST -c CLIENT_GENE_LIST -l LOG_FILE -o OUTPATH -t TAG -b BUFF -s SPA --LiPMScov LIPMSCOV

Process user specified arguments

options:
  -h, --help            show this help message and exit
  -Ag GENE_LIST, --gene_list GENE_LIST
                        path to all gene list to use
  -Eg ESS_GENE_LIST, --ess_gene_list ESS_GENE_LIST
                        path to essential gene list to use
  -Entg ENT_GENE_LIST, --ent_gene_list ENT_GENE_LIST
                        path to ent gene list to use
  -c CLIENT_GENE_LIST, --client_gene_list CLIENT_GENE_LIST
                        path to list of genes that are known clients
  -l LOG_FILE, --log_file LOG_FILE
                        Path to logging file
  -o OUTPATH, --outpath OUTPATH
                        path to output directory. will be made if doesnt exist
  -t TAG, --tag TAG     tag for final output image
  -b BUFF, --buff BUFF  buffer used C CD CG
  -s SPA, --spa SPA     spa used 0, 10, 20, 30, 40, 50, 60, 70, 80, 90
  --LiPMScov LIPMSCOV   LiPMS coverage used 0, 10, 20, 30, 40, 50, 60, 70, 80, 90
```

If you have the [SLUG] then you can use the command files located [here](src/command_lists/Assoc_Client_n_Essential.cmds) to reproduce the Fisher exact test used in this work to determine the association between protein essentiality and client state in the experimental data set and the AlphaFold structures. Please modify any other pathing as necessary. 