import os
import numpy as np
import glob

top_level = '/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/'
tags = os.listdir(top_level)
print(tags)

cmd_file = f'src/command_files/BuildKineticModels.cmds'
cmds = []

ITScmd_file = f'src/command_files/BuildKineticModels_ITS.cmds'
ITScmds = []

for tag in tags:
    print(f'Making BuildKinticModel for {tag}')

    MSMmapping = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/BuildKineticModel/{tag}_MSMmapping.csv'
    if os.path.exists(MSMmapping):
        continue
    ## check that all 50 G Q files are present
    counter = 0
    for t in np.arange(0,50):
        
        # check if G and Q file alread exist
        G = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/G/{tag}_t{t}.G'
        Q = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Q/{tag}_t{t}.Q'
        #print(G)
        #print(Q)
        if os.path.exists(G) and os.path.exists(Q):
            #print(f'GQ files FOUND for {tag} traj {t}:\n{Q}\n{G}')
            counter += 1
            continue

    print(f'Number of GQ files {counter}')
    if counter == 50:

        script = f'python src/data/BuildKineticModel.py'
        out = f' --outpath /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/BuildKineticModel/'
        OPpath = f' --OPpath /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/'
        misc = f'--lagtime 10 --outname {tag}'
        cmd = ' '.join([script, out, OPpath, misc])
        #print(cmd)
        cmds += [cmd]

        misc = misc + f'_ITS --ITS True'
        cmd = ' '.join([script, out, OPpath, misc])
        #print(cmd)
        ITScmds += [cmd]
        

## save non ITS command file
if len(cmds) != 0:
    np.savetxt(cmd_file, cmds, fmt='%s')
    print(f'SAVED: {cmd_file} {len(cmds)}')
else:
    print(f'No commands made for {tag} to save')

## save ITS command file
if len(ITScmds) != 0:
    np.savetxt(ITScmd_file, ITScmds, fmt='%s')
    print(f'SAVED: {ITScmd_file} {len(ITScmds)}')
else:
    print(f'No ITS commands made for {tag} to save')



## summary of files already existing
print(f'Summary of files already existing')
for tag in tags:

    MSMmapping = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/BuildKineticModel/{tag}_MSMmapping.csv'
    if os.path.exists(MSMmapping):
        Mcount = 1
    else:
        Mcount = 0

    ## check that all 50 G Q files are present
    Qcount = 0
    Gcount = 0
    for t in np.arange(0,50):
        
        # check if G and Q file alread exist
        G = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/G/{tag}_t{t}.G'
        Q = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Q/{tag}_t{t}.Q'
        MSMmapping = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/BuildKineticModel/{tag}_MSMmapping.csv'

        if os.path.exists(Q):
            Qcount += 1
            
        if os.path.exists(G):
            Gcount += 1

        
    print(f'{tag}: Qcount {Qcount} | Gcount {Gcount} | Mcount {Mcount} ')