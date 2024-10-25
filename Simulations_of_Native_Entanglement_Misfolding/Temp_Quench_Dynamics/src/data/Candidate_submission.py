import os
import numpy as np

cmd_files = np.loadtxt('temp_finished_CG.txt', dtype=str)
for cmdf in cmd_files:
    #print(cmdf)

    pdb = cmdf.split('/')[-1].split('_')[1]
    cmd = f'python src/data/submit_cmds.py {cmdf} data/temp.slurm {pdb} 50'
    print(cmd)

    os.system(cmd)
