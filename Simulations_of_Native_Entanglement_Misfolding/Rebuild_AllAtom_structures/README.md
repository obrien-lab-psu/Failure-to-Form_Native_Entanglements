# Rebuild missing residues in All-atom structures
PDB files often have missing residues in either long turn regions or with intrinsically disordered regions. As long as these regions are small it is not hard to rebuild them without impact the quality of the model. 
A list of our candidates to rebuild can be found [here](data/simulation_candidates_ids.csv).  

## Getting the PDB files for rebuilding

### Usage of [Get_PDBs.py](src/data/Get_PDBs.py)
```
usage: Get_PDBs.py [-h] --candidates CANDIDATES --outpath OUTPATH --log LOG

Process user specified arguments

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
                        Path to candidates file
  --outpath OUTPATH     Path to output directory
  --log LOG             Path to logging file
```

If you have the [SLUG] then you can use the command files located [here](src/command_files/Benchmarks_GPU_4YNG.cmds) to reproduce the benchmarks simulations. Please modify any other pathing as necessary. 

## Results
  
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