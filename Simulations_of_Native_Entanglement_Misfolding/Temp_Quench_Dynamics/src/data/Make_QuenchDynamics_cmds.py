import os
import numpy as np
import glob

top_level = '/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/'
tags = os.listdir(top_level)
print(tags)

for tag in tags:
    tag_file = f'src/command_files/{tag}_Quenching.cmds'
    cmds = []
    for t in np.arange(0,1):
        
        # find final PDB file
        #print(os.path.join(top_level, f'{tag}/Unfolding/*_t{t}_unfolding_finalframe*.pdb'))
        final_cor = glob.glob(os.path.join(top_level, f'{tag}/Unfolding/*_t{t}_unfolding_finalframe*.pdb'))
        if len(final_cor) > 1:
            raise ValueError(f'There was more than 1 final frame found for trajectory {t} in {tag}:\n{final_cor}')
        elif len(final_cor) == 1:
            final_cor = final_cor[0]
            #print(f'final_cor: {final_cor}')
            script = f'python src/data/Dynamics.py'
            psf = f'--psffile ../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/{tag}/setup/{tag}_rebuilt_clean_ca.psf'
            cor = f'--corfile {final_cor}'
            prm = f'--prmfile ../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/{tag}/setup/{tag}_rebuilt_clean_nscal1_fnn1_go_bt.xml'
            temp = f'--temp 310' 
            out = f' --outpath ../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/{tag}/Quenching/'
            misc = f'--outname {tag}_t{t}_quench --steps 133333334 --GPU True'
            cmd = ' '.join([script, psf, cor, prm, temp, out, misc])
            #print(cmd)

            cmds += [cmd]
        elif len(final_cor) == 0:
            print(f'Error: There was no final frame for  trajectory {t} in {tag}:\n{final_cor}')

    if len(cmds) != 0:
        np.savetxt(tag_file, cmds, fmt='%s')
        print(f'SAVED: {tag_file}')
    else:
        print(f'No commands made for {tag}')

    ## submit commands
    submission = f'python src/data/submit_cmds.py {tag_file} data/temp.slurm {tag.split("_")[1]} 1'
    print(submission)
    #os.system(submission)
