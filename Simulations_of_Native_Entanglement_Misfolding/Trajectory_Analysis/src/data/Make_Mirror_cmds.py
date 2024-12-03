import os
import numpy as np
import glob

top_level = '/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/'
tags = os.listdir(top_level)
#print(tags)

for tag in tags:
    tag_file = f'src/command_files/{tag}_Mirror.cmds'
    tag_native_file = f'src/command_files/{tag}_Native_Mirror.cmds'

    cmds = []
    native_cmds = []
    for t in np.arange(0,50):
        
        ######################################################################
        ## find quench dcd
        dcd = glob.glob(os.path.join(top_level, f'{tag}/Quenching/{tag}_t{t}_quench.dcd'))
        #print(f'Quenching dcd: {dcd}')

        ## check if output file already exists
        Kfile = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Mirror/Quench/K_{tag}_t{t}_quench.dat'
        #print(f'{Kfile}', os.path.exists(Kfile))
        if not os.path.exists(Kfile):
            ## make the quench commands
            if len(dcd) > 1:
                raise ValueError(f'There was more than 1 DCD found for trajectory {t} in {tag}:\n{dcd}')

            elif len(dcd) == 1:
                dcd = dcd[0]
                #calc_chirality_number.pl -i ./setup/1zmr_model_clean_ca.cor -d ./setup/domain_def.dat -s ./setup/secondary_struc_defs.txt -t ./prod_traj_reduce_saving/1_prod.dcd -o ./analysis/chirality/
                script = f'perl src/data/calc_chirality_number.pl'
                cor = f'-i /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/{tag}/setup/{tag}_rebuilt_clean_ca.cor'
                dcdstr = f'-t {dcd}'
                out = f' -o /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Mirror/Quench/'
                sec = f'-s /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/{tag}/setup/secondary_struc_defs.txt'
                dom = f'-d /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Simulations_of_Native_Entanglement_Misfolding/Rebuild_AllAtom_structures/data/domains/{tag}.txt'
                cmd = ' '.join([script, cor, dcdstr, sec, dom, out])
                #print(cmd)

                cmds += [cmd]

                ## check if output directory exists. If not make it
                if not os.path.exists(f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Mirror/Quench/'):
                    os.makedirs(f'/storage/group/epo2/default/ims86/git_slugs/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Mirror/Quench/')
                    #print(f'MADE: /storage/group/epo2/default/ims86/git_slugs/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Mirror/Quench/')
                

        ######################################################################
        ## find native dcd
        # find top level round of CG optimization
        CG_dir = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Coarse_Graining_Model/'
        rounds = os.listdir(f'{CG_dir}{tag}/')
        #print(rounds)
        round_numbers = [int(d.split('_')[1]) for d in rounds if d.startswith('round')]
        max_round = max(round_numbers)
        #print(tag, round_numbers, max_round)
        native_dcd = glob.glob(os.path.join(CG_dir, f'{tag}/round_{max_round}/{t}.dcd'))
        #print(f'native_dcd: {native_dcd} {len(native_dcd)}')

        ## check if output file already exists
        Kfile = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Mirror/Native/K_{t}.dat'
        #print(f'{Kfile}', os.path.exists(Kfile))
        if not os.path.exists(Kfile):
            ## make the native commandas
            if len(native_dcd) == 0:
                #print(f'No native dcd found for traj {t}')
                continue

            elif len(native_dcd) == 1:
                dcd = native_dcd[0]
                #calc_chirality_number.pl -i ./setup/1zmr_model_clean_ca.cor -d ./setup/domain_def.dat -s ./setup/secondary_struc_defs.txt -t ./prod_traj_reduce_saving/1_prod.dcd -o ./analysis/chirality/
                script = f'perl src/data/calc_chirality_number.pl'
                cor = f'-i /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/{tag}/setup/{tag}_rebuilt_clean_ca.cor'
                dcdstr = f'-t {dcd}'
                out = f' -o /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Mirror/Native/'
                sec = f'-s /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/{tag}/setup/secondary_struc_defs.txt'
                dom = f'-d /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Simulations_of_Native_Entanglement_Misfolding/Rebuild_AllAtom_structures/data/domains/{tag}.txt'
                cmd = ' '.join([script, cor, dcdstr, sec, dom, out])
                #print(cmd)

                native_cmds += [cmd]

                ## check if output directory exists. If not make it
                if not os.path.exists(f'/storage/group/epo2/default/ims86/git_slugs/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Mirror/Native/'):
                    os.makedirs(f'/storage/group/epo2/default/ims86/git_slugs/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Mirror/Native/')
                    #print(f'MADE: /storage/group/epo2/default/ims86/git_slugs/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Mirror/Native/')

    if len(cmds) != 0:
        np.savetxt(tag_file, cmds, fmt='%s')
        print(f'SAVED: {tag_file} {len(cmds)}')
    else:
        print(f'No commands made for {tag} to save')


    if len(native_cmds) != 0:
        np.savetxt(tag_native_file, native_cmds, fmt='%s')
        print(f'SAVED: {tag_native_file} {len(native_cmds)}')
    else:
        print(f'No commands made for Native {tag} to save')

## summarize the finished files
for i, tag in enumerate(tags):
    Quench_count = 0
    for t in np.arange(0,50):
        Kfile = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Mirror/Quench/K_{tag}_t{t}_quench.dat'
        if os.path.exists(Kfile):
            Quench_count += 1

    Native_count = 0
    for t in np.arange(0,50):
        Kfile = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/Mirror/Native/K_{t}.dat'
        if os.path.exists(Kfile):
            Native_count += 1
    print(f'{i} {tag} {Quench_count} {Native_count}')