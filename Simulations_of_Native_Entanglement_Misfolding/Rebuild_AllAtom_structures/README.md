# Rebuild missing residues in All-atom structures
PDB files often have missing residues in either long turn regions or with intrinsically disordered regions. As long as these regions are small it is not hard to rebuild them without impact the quality of the model. 
A list of our candidates to rebuild can be found [here](data/simulation_candidates_ids.csv).  

## Getting the PDB files for rebuilding
Here we preprocess the PDB files and standardize them against the canonical FASTA sequence before building the C-alpha coarse grained models. The script [Get_PDBs.py](src/data/Get_PDBs.py) is meant to simplify a lot of this processing but the output still requires careful inspection. A high level overview of the processing is as follows:  
1. Each unique *gene* *pdb* *chain* in the candidates.csv file provided will be downloaded and only that chain will be saved.  
  
2. The PDB sequence of the chain from the PDB will be read using the *Biopython* module.  
  
3. The canonical FASTA seuqnece will be downloaded directly from [Uniprot](https://www.uniprot.org/) using the *requests* module and the rest API call.   
   For example [https://rest.uniprot.org/uniprotkb/P0AD61.fasta](https://rest.uniprot.org/uniprotkb/P0AD61.fasta).  
     
4. The two sequences are aligned using [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) and the best alignment is choosen by the maximal [BLAST score](https://www.nlm.nih.gov/ncbi/workshops/2023-08_BLAST_evol/blast_score.html). 
   
5. Missing residues and mutations are found but importantly this code does not handel insertions or true deletions!  
    NOTE: There are no insertions or true deletions in the dataset of proteins used in this work.  
      
6. Mutated residues are removed and remodeled with the canonical sequence.  
    CAUTION: should be used here as if there are extensive mutations that differ in physicochemical characteristics from the canonical residue then another PDB maybe a better choice than remodeling. If you wish to avoid rebuilding mutated residues change the variable *non_missing_mutant_residues_df* in [Get_PDBs.py](src/data/Get_PDBs.py) to include the MUTATION column in making the resid_mapping object.   
      
7. The missing and mutated residues are rebuilt using the [MODELLER](https://salilab.org/modeller/download_installation.html) which does require a free license to operate. To prevent non-realistic bonds between the rebuilt residues and the remainder of the protein (which remains fixed during the rebuilding processes) we include the residues on either side of the missing residues in the MODELLER template.  
    NOTE: if you are rebuilding long sections of the protein such as with IDR tails then modeller may need a little help to prevent these rebuilt sections from threading gaps in the protein.  See the *special_restraints* section of the *MyModel* class inside the *rebuild_missing_residues* function for examples on how to defined lowerbound distance restraints to prevent the rebuilt terminal tails from coming to close too the protein.  
          
    For a list of the missing and mutated residues rebuilt in this project please see the log file located [here](data/Get_Rebuilt_res.log).  
      
8.  Finally the [CATH](https://www.cathdb.info/) domain information is pulled from the [boundaries-seqreschopping](data/cath-domain-boundaries-seqreschopping.txt) and [domain-list](data/cath-domain-list.txt) files provided by the user and mapped to the fasta sequence. 
      
    | CATH code | CLASS | Description                  |
    |-----------|-------|------------------------------|
    | 1         | a     | Mainly alpha-helical         |
    | 2         | b     | Mainly beta-strand           |
    | 3         | c     | Alpha-helical and beta-strand|
    | -         | n     | not defined                  |
    
    By definition missing residues not resolved in the PDB are not present in the CATH domain data and thus a simple approach is employeed where any missing residue is assigned to the next already assigned CATH domain ahead of it in the sequence. Any residue still not mapped to a domain after this first pass is then assigned to the previous already assigned CATH domain. 
      
    This approach assumes the following: (1) There are only short sections of missing residues within already defined domains and thus no matter which direction you scan the previous or next official CATH domain will be the same. (2) Most long stretches or missing residues are located on the terminal tails. 
      
    Therefore, if you have a structure that violates these assumptions you may need to manually curate the domain.txt file created. 
      
    Finally if no CATH data is present then the [STRIDE](https://webclu.bio.wzw.tum.de/stride/) program is used to estimate the secondary structure content of the domain and classify it based on a 20% threshold. If the protein has both alpha-helical and beta-strand content above 20% then the structure is classified as mixed a/b and if neither meets the threshold then it is assigned a value of *n* and should be manually curated by the user.  
      
    Two of the candidates (P0AES0 and P37747) had CATH data available for a portion of the resolved PDB structure but was missing annotations for another significant portion of the structure. Again we used the STRIDE method to approximate the structural class of these missing domains.  


### Usage of [Get_PDBs.py](src/data/Get_PDBs.py)
You will need the following external software installed and accessible in your PATH to run this code: 
1.  [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)  
2.  [STRIDE](https://webclu.bio.wzw.tum.de/stride/)  
3.  [MODELLER](https://salilab.org/modeller/download_installation.html) python package which does require a free license to operate. (does not need to be in your PATH but does need to be installed in your current conda env.)   
4.  Download the most up-to-date CATH files located [here](https://www.cathdb.info/wiki?id=data:index). You will need the cath-domain-boundaries-seqreschopping.txt and cath-domain-list.txt files.  

```
usage: Get_PDBs.py [-h] --candidates CANDIDATES --outpath OUTPATH --log LOG --cath CATH --cath_domain_list CATH_DOMAIN_LIST

Process user specified arguments

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
                        Path to candidates file
  --outpath OUTPATH     Path to output directory
  --log LOG             Path to logging file
  --cath CATH           Path to latest cath-domain-boundaries-seqreschopping.txt file
  --cath_domain_list CATH_DOMAIN_LIST
                        Path to latest cath-domain-list.txt file
```

If you have the [SLUG] then you can use the command files located [here](src/command_files/Get_PDBs.cmd) to processes and rebuild the PDBs in this work. Please modify any other pathing as necessary. 


## Results
Below are the results of the rebuilding process (red) overlayed ontop of the original PDB (blue). While non-rebuilt residues are held fixed during the rebuilding process newly added residues can change how nglview renders cartoon secondary structure so there maybe changes in the representation that are unexpected. Make sure you check each structure in the *ball-and-stick* model in VMD or Pymol to valid the rebuild.  

data/gifs/P0A6B4_4WR3_A_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P0A6B4_4WR3_A_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P0A6E6_1AQT_A_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P0A6E6_1AQT_A_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P0A6I0_1KDO_A_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P0A6I0_1KDO_A_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P0A6L2_1DHP_A_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P0A6L2_1DHP_A_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P0A6T9_3A7L_A_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P0A6T9_3A7L_A_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P0A763_2HUR_C_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P0A763_2HUR_C_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P0A790_4CRZ_A_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P0A790_4CRZ_A_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P0A7J3_6XZ7_H_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P0A7J3_6XZ7_H_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P0A7N9_6XZ7_b_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P0A7N9_6XZ7_b_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P0A8I5_3DXX_A_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P0A8I5_3DXX_A_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P0AA25_6LUR_D_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P0AA25_6LUR_D_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P0ACP7_1QPZ_A_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P0ACP7_1QPZ_A_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P0ADY3_6XZ7_K_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P0ADY3_6XZ7_K_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P0ADZ0_6QUL_U_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P0ADZ0_6QUL_U_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P0AES0_2IO9_B_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P0AES0_2IO9_B_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P0AG51_6XZ7_Z_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P0AG51_6XZ7_Z_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P21599_6K0K_A_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P21599_6K0K_A_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P31142_1URH_A_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P31142_1URH_A_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P37747_1I8T_A_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P37747_1I8T_A_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P45748_1HRU_A_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P45748_1HRU_A_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P60438_6PCR_N_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P60438_6PCR_N_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P60546_2ANB_A_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P60546_2ANB_A_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P61175_6XZ7_S_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P61175_6XZ7_S_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P65556_2FKB_A_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P65556_2FKB_A_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P77214_1ZEQ_X_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P77214_1ZEQ_X_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/P77754_6BIE_A_rebuilt.gif
<p align='center'>
  <img src=data/gifs/P77754_6BIE_A_rebuilt.gif alt='My GIF'>
</p>  
  
data/gifs/Q46856_1OJ7_D_rebuilt.gif
<p align='center'>
  <img src=data/gifs/Q46856_1OJ7_D_rebuilt.gif alt='My GIF'>
</p>  


### Usage of [Get_Rebuilt_res.py](src/data/Get_Rebuilt_res.py)
This script just checks for all those residues that were rebuilt in this study and for any Insertions or Deletions reported in the resid_mapping.txt files provided by our previous [publication](https://pubmed.ncbi.nlm.nih.gov/38296158/).  
```
usage: Get_Rebuilt_res.py [-h] --inpfiles INPFILES --mapping MAPPING --outpath OUTPATH

Process user specified arguments

options:
  -h, --help           show this help message and exit
  --inpfiles INPFILES  Path to candidates file
  --mapping MAPPING    Path to mapping files to check for insertions or deletions
  --outpath OUTPATH    Path to output directory
```
The resulting log file with details can be viewed [here](data/Get_Rebuilt_res.log).  
  
If you have the [SLUG] then you can use the command files located [here](src/command_files/Get_Rebuilt_res.cmds) to reproduce the set of missing residues rebuilt in this work. Please modify any other pathing as necessary. 