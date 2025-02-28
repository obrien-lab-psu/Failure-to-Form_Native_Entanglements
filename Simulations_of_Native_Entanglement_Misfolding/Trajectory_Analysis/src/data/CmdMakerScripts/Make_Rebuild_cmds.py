import os
import numpy as np
import glob
import pandas as pd

top_level = '/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/'
tags = os.listdir(top_level)
print(tags)


candidates = pd.read_csv('data/simulation_candidates_ids.csv')
candidates = candidates[candidates['set'].isin([0,3])]
print(candidates)

#########################################################################################################
cmd_file = f'src/command_files/Rebuild_last_frames.cmds'
cmds = []

for gene, pdb, chain, set in candidates[['gene', 'pdb', 'chain', 'set']].values:
    tag = f'{gene}_{pdb}_{chain}'
    print(gene, pdb, chain, tag)

    for t in np.arange(0,50):
        print(f'Making backmapping commands for {tag} {t}')
        
        script = f"python src/data/backmap.py" 
        AA = f'--aa_pdb ../Rebuild_AllAtom_structures/data/post_rebuilt/{tag}_rebuilt.pdb'
        CG = glob.glob(f'../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/{tag}/Quenching/{tag}_t{t}_quench_finalframe*.pdb')[0]
        if not os.path.exists(CG):
            print(f'CG file for {tag} {t} does NOT exist!')
            quit()
        CG = f"--cg_pdb {CG}"
        misc = '-p 1'
        cmd = ' '.join([script, AA, CG, misc])
        #print(cmd)
        cmds += [cmd]
        
## save cmd file if theere are any to save
if len(cmds) != 0:
    np.savetxt(cmd_file, cmds, fmt='%s')
    print(f'SAVED: {cmd_file} {len(cmds)}')
else:
    print(f'No commands made for to save')
###############################################