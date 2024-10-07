# Benchmarking C-alpha Coarse grained molecular dynamics sims on PSU ROAR GPUs

### Usage of [dynamics.py](src/data/dynamics.py)
```
usage: dynamics.py [-h] --psffile PSFFILE --corfile CORFILE --prmfile PRMFILE --temp TEMP --outpath OUTPATH --outname OUTNAME --log LOG --steps STEPS --GPU GPU [--restart RESTART]

Process user specified arguments

options:
  -h, --help         show this help message and exit
  --psffile PSFFILE  Path to protein structure file
  --corfile CORFILE  Path to protein coordinate file
  --prmfile PRMFILE  Path to force field file .xlm format
  --temp TEMP        Temperature in K to run the simulation at
  --outpath OUTPATH  Path to output directory
  --outname OUTNAME  name for output files
  --log LOG          Path to logging file
  --steps STEPS      Number of simulation steps to run
  --GPU GPU          True use GPU | False use CPUs
  --restart RESTART  True to restart | False to start from scratch (optional - if not specified assumed to be False)
```

If you have the [SLUG] then you can use the command files located [here](src/command_files/Benchmarks_GPU_4YNG.cmds) to reproduce the benchmarks simulations. Please modify any other pathing as necessary. 

## Results
After running 30 trajectories on the mgc-nih GPUs for 105ns. (really ran 10 trajectories and restarted them each twice to test restart capabilites). I obtained the run times located [here](data/benchmark_times.txt) that have a mean run time of 28.53 (28.50, 28.57) min for a rate of ~4.53 hr/us for the 470 residue model.  
