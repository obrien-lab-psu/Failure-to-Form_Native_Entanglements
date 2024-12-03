import os
import numpy as np
import glob

top_level = '/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/'
tags = os.listdir(top_level)
print(tags)

cmd_file = f'src/command_files/Cluster_ChangesInEnt.cmds'
cmds = []

for tag in tags:
    print(f'Making Cluster_ChangesInEnt for {tag}')

    for t in np.arange(0,50):
        
        # check if EntInfo input file exists and if this clustered output file already exists
        EntInfo = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/G/{tag}_t{t}.EntInfo'
        outfile = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Cluster_ChangesInEnt/{tag}_t{t}_clustered.EntInfo'
        if os.path.exists(outfile):
            #print(f'{outfile} ALREADY EXISTS')
            continue

        if os.path.exists(EntInfo):

            script = f'python src/data/Cluster_ChangesInEnt.py'
            outdir = f'--outpath /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Cluster_ChangesInEnt/'
            EntInfo = f'--EntInfofile /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/G/{tag}_t{t}.EntInfo'
            misc = f'--outname {tag}_t{t}'
            cmd = ' '.join([script, outdir, EntInfo, misc])
            #print(cmd)
            cmds += [cmd]
        

## save non ITS command file
if len(cmds) != 0:
    np.savetxt(cmd_file, cmds, fmt='%s')
    print(f'SAVED: {cmd_file} {len(cmds)}')
else:
    print(f'No commands made for to save')

print('Normal Termination')
#python src/data/Cluster_ChangesInEnt.py --outpath ../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/Cluster_ChangesInEnt/ --EntInfofile ../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/P0A7J3_6XZ7_H/G/P0A7J3_6XZ7_H_t2.EntInfo --outname P0A7J3_6XZ7_H_t2
#python src/data/Cluster_ChangesInEnt.py --outpath ../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/Cluster_ChangesInEnt/ --EntInfofile ../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/P0A7J3_6XZ7_H/G/P0A7J3_6XZ7_H_t2.EntInfo --outname P0A7J3_6XZ7_H_t2
#python src/data/Cluster_ChangesInEnt.py --outpath ../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/Cluster_ChangesInEnt/ --EntInfofile ../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/P0A6E6_1AQT_A/G/P0A6E6_1AQT_A_t16.EntInfo --outname P0A6E6_1AQT_A_t16

## summary of files already present
print(f'Already made summary')
for tag in tags:
    counter = 0
    ent_counter = 0
    for t in np.arange(0,50):
        # check if EntInfo input file exists and if this clustered output file already exists
        EntInfo = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/G/{tag}_t{t}.EntInfo'
        if os.path.exists(EntInfo):
            ent_counter += 1

        outfile = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Cluster_ChangesInEnt/{tag}_t{t}_clustered.EntInfo'
        if os.path.exists(outfile):
            counter += 1

    print(f'{tag} {ent_counter} {counter}')
