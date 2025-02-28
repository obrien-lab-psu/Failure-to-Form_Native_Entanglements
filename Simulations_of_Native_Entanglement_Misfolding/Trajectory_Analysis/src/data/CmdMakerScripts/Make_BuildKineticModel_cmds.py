import os
import numpy as np
import glob
import pandas as pd

top_level = '/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/'
tags = os.listdir(top_level)
print(tags)

cmd_file = f'src/command_files/BuildKineticModels.cmds'
cmds = []

candidates = pd.read_csv('data/simulation_candidates_ids.csv')
for rowi, row in candidates.iterrows():
    gene = row['gene']
    pdb = row['pdb']
    chain = row['chain']
    tag = f'{gene}_{pdb}_{chain}'
    setID = row['set']
    print(f'Making BuildKinticModel for {tag}')
    print(setID, type(setID))

    MSMmapping = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/BuildKineticModel/setID{setID}/{tag}_MSMmapping.csv'
    #if os.path.exists(MSMmapping):
    #    continue
    ## check that all 50 G Q files are present
    Qcounter = 0
    Gcounter = 0
    for t in np.arange(0,50):
        
        # check if G and Q file alread exist
        G = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Cluster_ChangesInEnt/{tag}_t{t}_clustered.G'
        Q = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Q/{tag}_t{t}.Q'
        #print(G)
        #print(Q)
        if os.path.exists(Q):
            Qcounter += 1
        else:
            print(Q)

        if os.path.exists(G):
            Gcounter += 1
        else:
            print(G)
            

    print(f'Number of Q files {Qcounter} | Number of G files {Gcounter}')
    if Gcounter == 50 and Qcounter == 50:

        script = f'python src/data/BuildKineticModel.py'
        out = f' --outpath /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/BuildKineticModel/setID{setID}/'
        OPpath = f' --OPpath /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/'
        misc = f'--lagtime 10 --outname {tag}'
        start_end = f'--start 24000 --end 26666'
        
        if setID == 3:
            cmd = ' '.join([script, out, OPpath, misc, start_end])
        else:
            cmd = ' '.join([script, out, OPpath, misc])
        #print(cmd)
        cmds += [cmd]
        

## save non ITS command file
if len(cmds) != 0:
    np.savetxt(cmd_file, cmds, fmt='%s')
    print(f'SAVED: {cmd_file} {len(cmds)}')
else:
    print(f'No commands made for {tag} to save')




## summary of files already existing
print(f'Summary of files already existing')
for rowi, row in candidates.iterrows():
    gene = row['gene']
    pdb = row['pdb']
    chain = row['chain']
    tag = f'{gene}_{pdb}_{chain}'
    setID = row['set']

    MSMmapping = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/BuildKineticModel/setID{setID}/{tag}_MSMmapping.csv'
    if os.path.exists(MSMmapping):
        Mcount = 1
    else:
        Mcount = 0

    ## check that all 50 G Q files are present
    Qcount = 0
    Gcount = 0
    for t in np.arange(0,50):
        
        # check if G and Q file alread exist
        G = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Cluster_ChangesInEnt/{tag}_t{t}_clustered.G'
        Q = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Q/{tag}_t{t}.Q'

        if os.path.exists(Q):
            Qcount += 1
            
        if os.path.exists(G):
            Gcount += 1

        
    print(f'{tag} setID {setID}: Qcount {Qcount} | Gcount {Gcount} | Mcount {Mcount} ')