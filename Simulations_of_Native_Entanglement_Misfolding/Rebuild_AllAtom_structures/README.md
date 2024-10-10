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

P0A6B4_4WR3_A  
<p align="center">
  <img src="data/gifs/P0A6B4_4WR3_A_rebuilt.gif" alt="My GIF">
</p>  