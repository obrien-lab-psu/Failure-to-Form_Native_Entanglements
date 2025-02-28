import os
import numpy as np
import glob

top_level = '/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/'
tags = os.listdir(top_level)
#print(tags)

for tag in tags:
    tag_file = f'src/command_files/{tag}_GQ.cmds'
    cmds = []
    for t in np.arange(0,50):
        
        # check if G and Q file alread exist
        G = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/G/{tag}_t{t}.G'
        Q = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Q/{tag}_t{t}.Q'
        found = False
        if os.path.exists(G):
            #print(f'G files FOUND for {tag} traj {t}: {G}')
            if os.path.exists(Q):
                #print(f'Q files FOUND for {tag} traj {t}: {Q}')
                found = True

        if found:
            continue


        # find final PDB file
        #print(os.path.join(top_level, f'{tag}/Unfolding/*_t{t}_unfolding_finalframe*.pdb'))
        dcd = glob.glob(os.path.join(top_level, f'{tag}/Quenching/{tag}_t{t}_quench.dcd'))
        #print(dcd)
        if len(dcd) > 1:
            raise ValueError(f'There was more than 1 DCD found for trajectory {t} in {tag}:\n{dcd}')

        elif len(dcd) == 1:
            dcd = dcd[0]

            script = f'python src/data/GQ.py'
            psf = f'--psf /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/{tag}/setup/{tag}_rebuilt_clean_ca.psf'
            cor = f'--cor /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/{tag}/setup/{tag}_rebuilt_clean_ca.cor'
            dcdstr = f'--dcd {dcd}'
            out = f' --outpath /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/'
            sec = f'--sec_elements /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/{tag}/setup/secondary_struc_defs.txt'
            misc = f'--outname {tag}_t{t} --start -2667'
            cmd = ' '.join([script, psf, cor, dcdstr, sec, out, misc])
            #print(cmd)

            cmds += [cmd]


    if len(cmds) != 0:
        np.savetxt(tag_file, cmds, fmt='%s')
        print(f'SAVED: {tag_file} {len(cmds)}')
    else:
        print(f'No commands made for {tag} to save')

## summarize the finished files
for i, tag in enumerate(tags):
    Q_count = 0
    G_count = 0

    for t in np.arange(0,50):
        # check if G and Q file alread exist
        G = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/G/{tag}_t{t}.G'
        Q = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Q/{tag}_t{t}.Q'
        if os.path.exists(G):
            G_count += 1
        if os.path.exists(Q):
            Q_count += 1
              
    print(f'{i} {tag} {G_count} {Q_count}')

