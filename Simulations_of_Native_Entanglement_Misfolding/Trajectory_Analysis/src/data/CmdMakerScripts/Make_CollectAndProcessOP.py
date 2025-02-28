import os
import numpy as np
import glob
import pandas as pd



cmds = {1:[], 2:[], 3:[]}

candidates = pd.read_csv('data/simulation_candidates_ids.csv')
for rowi, row in candidates.iterrows():
    gene = row['gene']
    pdb = row['pdb']
    chain = row['chain']
    tag = f'{gene}_{pdb}_{chain}'
    setID = row['set']


    cmd = f'python src/data/CollectAndProcessOP.py --outpath ../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CollectAndProcessOP/setID{setID}/ --candidates data/simulation_candidates_ids.csv --toplevel ../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/ --outname {tag} --Mirrorfile data/annotated_mirrors.csv --setID {setID} --gene {gene}'
    cmds[setID] += [cmd]


for setID in [1, 2, 3]:
    if len(cmds[setID]) != 0:
        cmd_file = f'src/command_files/CollectAndProcessOP_setID{setID}.cmds'
        np.savetxt(cmd_file, cmds[setID], fmt='%s')
        print(f'SAVED: {cmd_file} {len(cmds[setID])}')
