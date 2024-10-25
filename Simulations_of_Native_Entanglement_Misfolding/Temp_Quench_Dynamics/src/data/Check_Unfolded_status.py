import glob
import pandas as pd
import MDAnalysis as mda
import os
import shutil
import glob

list_dirs = os.listdir('../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/')
print(list_dirs)

check = {'index':[], 'gene':[], 'pdb':[], 'chain':[], 'Num_unfolded':[]}
cmds = []
for i, tag in enumerate(list_dirs):
    gene, pdb, chain = tag.split('_')
    #print(i, tag, gene, pdb, chain)

    # get list of all final frames
    finalF = glob.glob(f'../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/{tag}/Unfolding/*final*.pdb')
    #print(f'Number of final frames: {len(finalF)}')
    completed_traj = []
    for t in range(0,50):
        f = [x for x in finalF if f't{t}_' in x]

        # if the final frame is not even present get the command and increase the steps just for the fun of it!
        if len(f) != 1:
            print(f'No final frame found for traj {t} in {tag}')
            #get the command from the command files and double to unfolding time
            cmd = list(os.popen(f'grep {tag}_t{t}_ src/command_files/{tag}_unfolding.cmds'))[0]
            cmd = cmd.replace('66666667', '266666668').replace('\n', '') + ' --restart True\n'
            #print(cmd)  
            cmds += [cmd]
           
        # if the final frame is present check if it was at the last step and if so increase it by double again
        else:
            f = f[0]
            frame = f.split('/')[-1].split('_')[-1].replace('finalframe', '').replace('.pdb','')
            traj = f.split('/')[-1].split('_')[3].replace('t','')
            #print(f, frame, traj)
            if frame == '13334':
                print(f'Final frame found for traj {t} in {tag} but was at 1us indicating the protein did not unfold')
                print(f, frame, traj)

                #get the command from the command files and double to unfolding time
                cmd = list(os.popen(f'grep {tag}_t{t}_ src/command_files/{tag}_unfolding.cmds'))[0]
                cmd = cmd.replace('66666667', '266666668').replace('\n', '') + ' --restart True\n'                #print(cmd)  
                cmds += [cmd]

            else:
                completed_traj += [traj]
        #print(f'Number of completed final frames: {len(completed_traj)}')

    check['index'] += [i]
    check['gene'] += [gene]
    check['pdb'] += [pdb]
    check['chain'] += [chain]
    check['Num_unfolded'] += [len(completed_traj)]

check = pd.DataFrame(check)
print(check.to_string())
check.to_csv('data/Check_Unfolding.csv', index=False)
print(f'SAVED: data/Check_Unfolding.csv')

print(len(cmds))
with open('restart.cmds', 'w') as fh:
    for cmd in cmds:
        fh.write(cmd)
print(f'SAVED: restart.cmds')
    