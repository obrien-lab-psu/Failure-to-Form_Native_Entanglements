import os
import numpy as np

tags = os.listdir('/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Coarse_Graining_Model/')
print(tags)

for tag in tags:
    tag_file = f'src/command_files/{tag}_unfolding.cmds'
    cmds = []
    for t in np.arange(0,50):
        script = f'python src/data/UnfoldingDynamics.py'
        psf = f'--psffile ../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/{tag}/setup/{tag}_rebuilt_clean_ca.psf'
        cor = f'--corfile ../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/{tag}/setup/{tag}_rebuilt_clean_ca.cor'
        prm = f'--prmfile ../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/{tag}/setup/{tag}_rebuilt_clean_nscal1_fnn1_go_bt.xml'
        temp = f'--temp 800' 
        sec = f'--sec_elements ../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/{tag}/setup/secondary_struc_defs.txt'
        out = f' --outpath ../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/{tag}/Unfolding/'
        misc = f'--outname {tag}_t{t}_unfolding --steps 66666667 --GPU True'
        cmd = ' '.join([script, psf, cor, prm, temp, sec, out, misc])
        #print(cmd)

        cmds += [cmd]
    np.savetxt(tag_file, cmds, fmt='%s')
    print(f'SAVED: {tag_file}')

