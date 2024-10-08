import os
import glob
import pandas as pd
import numpy as np

df = pd.read_csv('data/simulation_candidates_ids.csv')
print(df)

array = np.asarray(df[['gene', 'pdb', 'chain']].values, dtype=str)
array = np.unique(array, axis=0)
print(array, array.shape)

cmds = []
for gene, pdb, chain in array:
    pre_pdb = f'data/pre_rebuilt/{gene}_{pdb}_{chain}.pdb'
    post_pdb = f'data/post_rebuilt/{gene}_{pdb}_{chain}_rebuilt.pdb'

    if os.path.exists(pre_pdb) and os.path.exists(post_pdb):

        cmd = f'python src/data/Make_gif.py --prebuilt_PDB {pre_pdb} --postbuilt_PDB {post_pdb} --outpath data/gifs/'
        if cmd not in cmds:
            print(cmd)
            cmds += [cmd]