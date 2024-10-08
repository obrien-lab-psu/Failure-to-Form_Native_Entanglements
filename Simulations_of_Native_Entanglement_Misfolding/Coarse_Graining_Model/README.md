# Coarse_graing_protein_models

## C-alpha coarse graining of PDB all-atom models


### Usage of [opt_nscal.py](src/data/opt_nscal.py)
```

  Usage: python opt_nscal.py
                --input | -i <input.pdb> for CG model creation
                --domain | -d <domain.dat> for domain defination
                [--temp | -t] <Temperature> in Kelvin
                [--ppn | -n] <number of CPUs> for each trajectory. Default 1.
                        if 0, use a single GPU for each trajectory.
                [--tpn | -p] <total number of CPUs>. Default 10.
                [--ntraj | -j] <number of trajectories>. Default 10.
                [--restart | -r] <0 or 1> restart optimization. Default 0, not restart.
                [--nscal | -s] <nscal_level.dat> for nscal levels. Default values were
                               obtained from a training set of 18 small single-domain
                               proteins.
                [--casm | -c] <0 or 1> CG model type. Default 0, C-alpha model.
                              1, C-alpha side chain model.
                [--help | -h]

  Example domain.dat:
  1:96 302:350 a #Domain 1 is from resid 1 to 96 and 302 to 350, and in
                 #alpha-helix class
  97:155 b #Domain 2 is from resid 97 to 155 and in beta-sheet class
  156:301 c #Domain 3 is from resid 156 to 301 and in alpha-beta class

  Example nscal_level.dat:
  a 1.1954 1.4704 1.7453 2.0322 2.5044 1.7453
  b 1.4732 1.8120 2.1508 2.5044 2.5044 2.1508
  c 1.1556 1.4213 1.6871 1.9644 2.5044 1.6871
  i 1.2747 1.5679 1.8611 2.1670 2.5044 1.8611

```

If you have the [SLUG] then you can use the command files located [here](src/command_files/opt_nscal_cmds.sh) to reproduce the coarse grained protein files. Please modify any other pathing as necessary. 

## Results

