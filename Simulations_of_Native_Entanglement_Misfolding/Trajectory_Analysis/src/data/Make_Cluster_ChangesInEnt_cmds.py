import os
import numpy as np
import glob

top_level = '/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/'
tags = os.listdir(top_level)
print(tags)


#########################################################################################################
cmd_file = f'src/command_files/Cluster_ChangesInEnt.cmds'
cmds = []

for tag in tags:
    print(f'Making Cluster_ChangesInEnt for {tag}')
    gene, _, _ = tag.split('_')

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
            EntInfo = f'--EntInfofile {EntInfo}'
            misc = f'--outname {tag}_t{t}'
            cmd = ' '.join([script, outdir, EntInfo, misc])
            #print(cmd)
            cmds += [cmd]
        
## save cmd file if theere are any to save
if len(cmds) != 0:
    np.savetxt(cmd_file, cmds, fmt='%s')
    print(f'SAVED: {cmd_file} {len(cmds)}')
else:
    print(f'No commands made for to save')
#########################################################################################################


#########################################################################################################
#### make the native clustering files
cmd_file = f'src/command_files/Cluster_Native_ChangesInEnt.cmds'
cmds = []

for tag in tags:
    print(f'Making Cluster_ChangesInEnt for {tag}')
    gene, _, _ = tag.split('_')

    for t in np.arange(0,50):
        
        # check if EntInfo input file exists and if this clustered output file already exists
        EntInfo = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Native/G/{tag}_t{t}.EntInfo'
        outfile = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Native/Cluster_ChangesInEnt/{tag}_t{t}_clustered.EntInfo'
        if os.path.exists(outfile):
            #print(f'{outfile} ALREADY EXISTS')
            continue

        if os.path.exists(EntInfo):

            script = f'python src/data/Cluster_ChangesInEnt.py'
            outdir = f'--outpath /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Native/Cluster_ChangesInEnt/'
            EntInfo = f'--EntInfofile {EntInfo}'
            misc = f'--outname {tag}_t{t}'
            cmd = ' '.join([script, outdir, EntInfo, misc])
            #print(cmd)
            cmds += [cmd]
        
## save cmd file if theere are any to save
if len(cmds) != 0:
    np.savetxt(cmd_file, cmds, fmt='%s')
    print(f'SAVED: {cmd_file} {len(cmds)}')
else:
    print(f'No commands made for to save')
#########################################################################################################

## summary of files already present
print(f'Already made summary')
for tag in tags:
    counter = 0
    ent_counter = 0
    Ncounter = 0
    Nent_counter = 0
    for t in np.arange(0,50):
        # check if EntInfo input file exists and if this clustered output file already exists
        EntInfo = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/G/{tag}_t{t}.EntInfo'
        if os.path.exists(EntInfo):
            ent_counter += 1
        nEntInfo = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Native/G/{tag}_t{t}.EntInfo'
        if os.path.exists(nEntInfo):
            Nent_counter += 1

        outfile = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Cluster_ChangesInEnt/{tag}_t{t}_clustered.EntInfo'
        if os.path.exists(outfile):
            counter += 1
        noutfile = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Native/Cluster_ChangesInEnt/{tag}_t{t}_clustered.EntInfo'
        if os.path.exists(noutfile):
            Ncounter += 1

    print(f'{tag} | Quench {ent_counter} {counter} | Native {Nent_counter} {Ncounter}')
