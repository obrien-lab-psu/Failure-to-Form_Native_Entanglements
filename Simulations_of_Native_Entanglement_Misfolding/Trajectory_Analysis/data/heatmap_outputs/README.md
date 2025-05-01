## Test for convergence of Calpha coarse grained models used in study  
There are three sub directories for each of the three datasets analyzed.  
[setID1](data/heatmap_outputs/setID1) is the set of refoldable proteins   
[setID2](data/heatmap_outputs/setID2) is the set of non-refoldable proteins size matched to setID1  
[setID3](data/heatmap_outputs/setID3) is the set of non-refoldable proteins that reflects the whole proteome size distribution  

Within each folder is a figure for each of the proteins simulated in datasets 1 and 3.  Each figure has 8 subplots cooresponding to the running average and running standard deviation of the four key order parameters used in this work. Fraction of native contacts (Q), Fraction of native contacts with a change in entanglement (G), Fraction of secondary structures that packed with correct chirality (K), relative difference in solvant accessible surface area (Z, or zeta).  

For data sets 1 and 2 the full trajectories are shown and the white dashed line indicates where we truncated the data to examine just the last steady state portion of each trajectory.  For dataset 3 