import os
import numpy as np
import glob

top_level = '/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/'
tags = os.listdir(top_level)
CG_dir = f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Coarse_Graining_Model/'
#print(tags)
cmds = []
outfile = f'src/command_files/Zeta.cmds'
for tag in tags:

    # find top level round of CG optimization
    rounds = os.listdir(f'{CG_dir}{tag}/')
    round_numbers = [int(d.split('_')[1]) for d in rounds if d.startswith('round')]
    max_round = max(round_numbers)
    print(tag, round_numbers, max_round)


    script = f'python src/data/SASA.py'
    psf = f'--psf /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/{tag}/setup/{tag}_rebuilt_clean_ca.psf'
    cor = f'--cor /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/{tag}/setup/{tag}_rebuilt_clean_ca.cor'
    quench_dcdstr = f'--quench_dcds /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/{tag}/Quenching/{tag}_t\*_quench.dcd'
    native_dcdstr = f'--native_dcds /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Coarse_Graining_Model/{tag}/round_{max_round}/\*.dcd'
    out = f' --outpath /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{tag}/'
    misc = f'--outname {tag}'
    fasta = f'--FASTApath /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Rebuild_AllAtom_structures/FASTA/{tag}.fasta'
    cmd = ' '.join([script, psf, cor, quench_dcdstr, native_dcdstr, out, misc, fasta])
    #print(cmd)

    cmds += [cmd]


if len(cmds) != 0:
    np.savetxt(outfile, cmds, fmt='%s')
    print(f'SAVED: {outfile} {len(cmds)}')


