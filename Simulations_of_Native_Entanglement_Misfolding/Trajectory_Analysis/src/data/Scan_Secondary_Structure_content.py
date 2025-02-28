import pandas as pd
import numpy  as np
import os, sys
import glob
import mdtraj as md

#########################################################################################################
def get_stride_residues(pdbfile):
    stride_cmd = f'stride -o {pdbfile}'
    #print(f'CALL: {stride_cmd}')
    results = os.popen(stride_cmd)
    SS_res = []
    for resline in results:
        if resline.startswith('LOC'):
            resline = resline.split()
            SS = resline[1]
            i = int(resline[3])
            j = int(resline[6])
            if SS in ['AlphaHelix', 'Strand', '310Helix']:
                #print(resline, SS, i, j)
                SS_res += [np.arange(i, j + 1)]
    SS_res = set(np.hstack(SS_res))
    #print(f'SS_res: {SS_res}')
            
    return SS_res
#########################################################################################################


#########################################################################################################
top = '/storage/group/epo2/default/ims86/'
#AA = f'--aa_pdb ../Rebuild_AllAtom_structures/data/post_rebuilt/{tag}_rebuilt.pdb'

rebuilt_files = glob.glob(f'highQ_Misfoled_Candidates_rebuilt/*.pdb')
Fraction_SS_df = {'gene':[], 'traj':[], 'frame':[], 'FracSS':[]}
for rebuilt in rebuilt_files:
    gene = rebuilt.split('/')[-1].split('_')[0]
    traj = rebuilt.split('/')[-1].split('_')[1].replace('t', '')
    frame = rebuilt.split('/')[-1].split('_')[2].replace('f', '')
    print(rebuilt, gene, traj, frame)

    AAfile = glob.glob(f'../Rebuild_AllAtom_structures/data/post_rebuilt/{gene}*_rebuilt.pdb')[0]
    print(f'AAfile: {AAfile}')

    AA_SS = get_stride_residues(AAfile)
    rebuilt_SS = get_stride_residues(rebuilt)

    overlap = AA_SS.intersection(rebuilt_SS)
    FracSS = len(overlap) / len(AA_SS)
    print(f'AA_SS: {len(AA_SS)} | rebuilt_SS: {len(rebuilt_SS)} | overlap: {len(overlap)} | FracSS: {FracSS}')

    Fraction_SS_df['gene'] += [gene]
    Fraction_SS_df['traj'] += [traj]
    Fraction_SS_df['frame'] += [frame]
    Fraction_SS_df['FracSS'] += [FracSS]

Fraction_SS_df = pd.DataFrame(Fraction_SS_df)
Fraction_SS_df = Fraction_SS_df.sort_values(by=['FracSS'])
print(f'Fraction_SS_df:\n{Fraction_SS_df.to_string()}')
Fraction_SS_df.to_csv(f'data/LossGain_candidate_SS_setID3.csv', index=False)
print(f'SAVE: data/LossGain_candidate_SS_setID3.csv')
##############################################################################################
print(f'NORMAL TERMINATION')