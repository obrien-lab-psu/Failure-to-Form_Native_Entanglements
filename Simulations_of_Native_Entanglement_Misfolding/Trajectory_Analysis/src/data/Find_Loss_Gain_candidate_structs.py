import pandas as pd
import numpy  as np
import os, sys
import glob
import mdtraj as md

setID = 2
#########################################################################################################
def extract_frame_from_trajectory(psf_path, dcd_path, frame_index, output_pdb_path):
    """
    Extract a specific frame from a DCD trajectory and save it to a new PDB file.

    Args:
        psf_path (str): Path to the .psf file (structure file).
        dcd_path (str): Path to the .dcd file (trajectory file).
        frame_index (int): The index of the frame to extract (zero-based).
        output_pdb_path (str): Path to save the extracted frame as a .pdb file.
    """
    # Load the trajectory using MDTraj
    #print(f"Loading trajectory from {dcd_path} with structure from {psf_path}...")
    traj = md.load(dcd_path, top=psf_path)

    # Extract the requested frame (MDTraj uses zero-based indexing)
    #print(f"Extracting frame {frame_index}...")
    frame = traj[frame_index]
    #print(frame)

    # Save the extracted frame to a new PDB file
    #print(f"Saving extracted frame to {output_pdb_path}...")
    frame.save_pdb(output_pdb_path)

    print(f"Frame {frame_index} successfully saved to {output_pdb_path}.")

#########################################################################################################

#########################################################################################################
def make_single_rebuild_cmd(gene, traj, output_pdb_path):
    top_level = '/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/'
    tags = os.listdir(top_level)
    tag = [tag for tag in tags if gene in tag][0]
    print(f'tag: {tag}')

    # check if this is a loss or gain candidate
    script = f"python src/data/backmap.py" 
    AA = f'--aa_pdb ../Rebuild_AllAtom_structures/data/post_rebuilt/{tag}_rebuilt.pdb'
    CG = output_pdb_path
    if not os.path.exists(CG):
        print(f'CG file for {tag} {traj} does NOT exist!')
        quit()
    CG = f"--cg_pdb {CG}"
    misc = '-p 1'
    cmd = ' '.join([script, AA, CG, misc])
    print(cmd)
    return cmd
#########################################################################################################

#########################################################################################################
def get_OP_summary(gene, traj):
    ## get the Q file
    #../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/P31142_1URH_A/Q/P31142_1URH_A_t40.Q
    Qfile = glob.glob(f'../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{gene}*/Q/{gene}*_t{traj}.Q')[0]
    print(f'Qfile: {Qfile}')
    Q = pd.read_csv(Qfile)
    Q = Q['Q'].values[-2600:]
    maxQpos = np.argmax(Q)
    revpos = 2600 - maxQpos
    print(Q.shape, maxQpos, Q[maxQpos], revpos, Q[-revpos])

    ## Get DCD and PSF
    #../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/P31142_1URH_A/setup/P31142_1URH_A_rebuilt_clean_ca.cor
    #../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/P31142_1URH_A/setup/P31142_1URH_A_rebuilt_clean_ca.psf
    #../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/P31142_1URH_A/Quenching/P31142_1URH_A_t40_quench.dcd
    dcdfile = glob.glob(f'../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/{gene}*/Quenching/{gene}_*_t{traj}_quench.dcd')[0]
    psffile = glob.glob(f'../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/{gene}*/setup/{gene}_*_rebuilt_clean_ca.psf')[0]
    print(f'dcdfile: {dcdfile}')
    print(f'psffile: {psffile}')
    output_pdb_path = os.path.join(f'./highQ_Misfoled_Candidates/{gene}_t{traj}_f{-revpos}_Q{Q[-revpos]:.3f}.pdb')
    extract_frame_from_trajectory(psffile, dcdfile, -revpos, output_pdb_path)
    
    EntInfo_f = glob.glob(f'{top}git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/{gene}_*/Cluster_ChangesInEnt/{gene}_*_t{traj}_clustered.EntInfo')[0]
    EntInfo = pd.read_csv(EntInfo_f, low_memory=False)
    #print(EntInfo)
    frames = np.unique(EntInfo['Frame'].values)
    frame = frames[-revpos]
    print(f'frame: {frame}')
    
    frame_EntInfo = EntInfo[EntInfo['Frame'] == frame]
    print(f'frame_EntInfo:\n{frame_EntInfo.to_string()}')

    cmd = make_single_rebuild_cmd(gene, traj, output_pdb_path)
    return cmd
#########################################################################################################

#########################################################################################################
top = '/storage/group/epo2/default/ims86/'
## load the Misfolding propensity and <Q> rankings
MQ_rankings = pd.read_csv(f'{top}git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingMechanism/setID{setID}/DATA/Combined_and_Processed_threshold_metrics_setID{setID}.csv')
MQ_rankings = MQ_rankings[(MQ_rankings['<Q>'] >= 0.75) & (MQ_rankings['MisfoldingProp'] >= 0.7) & (MQ_rankings['num_misfolded_frames'] >= 2000)]
MQ_rankings = MQ_rankings.sort_values(by=['<Q>'], ascending=False)
print(f'MQ_rankings:\n{MQ_rankings.to_string()}')


LossGain_fractions = pd.read_csv(f'{top}git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingMechanism/setID{setID}/DATA/PerTraj_Mechanism_summary_setID{setID}.csv')
#print(LossGain_fractions)
# Drop the first two columns by selecting column indices
#LossGain_fractions = LossGain_fractions.drop(LossGain_fractions.columns[[0, 1]], axis=1)
# parse those within the thresholded MQ_rankings
parsed_fractions = []
for gene, traj in MQ_rankings[['gene', 'traj']].values:
    loc = LossGain_fractions[(LossGain_fractions['gene'] == gene) & (LossGain_fractions['traj'] == traj)]
    #print(loc)
    parsed_fractions += [loc]
LossGain_fractions = pd.concat(parsed_fractions)
print(f'LossGain_fractions:\n{LossGain_fractions.to_string()}')
##############################################################################################


##############################################################################################
print(f'\n{"#"*100}\n{"#"*100}\n')
## get Loss candidates
Loss_candidates = LossGain_fractions[LossGain_fractions['FracLoss'] > 0.8]
Loss_candidates = Loss_candidates.sort_values(by=['FracLoss'], ascending=False)
print(f'Loss_candidates:\n{Loss_candidates.to_string()} {len(Loss_candidates)}')

# generate rebuild commands for these candidates
cmds = []
cmd_file = f'src/command_files/Rebuild_last_frames_Loss_Candidates_setID{setID}.cmds'
for rowi, row in Loss_candidates.iterrows():
    gene, traj = row['gene'], row['traj']
    print(f'\n{"#"*100}\nLOSS {gene} {traj}')

    cmd = get_OP_summary(gene, traj)
    cmds += [cmd]
## save cmd file if theere are any to save
if len(cmds) != 0:
    np.savetxt(cmd_file, cmds, fmt='%s')
    print(f'SAVED: {cmd_file} {len(cmds)}')
else:
    print(f'No commands made for to save')
##############################################################################################


##############################################################################################
print(f'\n{"#"*100}\n{"#"*100}\n')
## get Gain candidates
Gain_candidates = LossGain_fractions[LossGain_fractions['FracGain'] > 0.8]
Gain_candidates = Gain_candidates.sort_values(by=['FracGain'], ascending=False)
print(f'Gain_candidates:\n{Gain_candidates.to_string()} {len(Gain_candidates)}')

# generate rebuild commands for these candidates
cmds = []
cmd_file = f'src/command_files/Rebuild_last_frames_Gain_Candidates_setID{setID}.cmds'
for rowi, row in Gain_candidates.iterrows():
    gene, traj = row['gene'], row['traj']
    print(f'\n{"#"*100}\nGAIN {gene} {traj}')

    cmd = get_OP_summary(gene, traj)
    cmds += [cmd]
## save cmd file if theere are any to save
if len(cmds) != 0:
    np.savetxt(cmd_file, cmds, fmt='%s')
    print(f'SAVED: {cmd_file} {len(cmds)}')
else:
    print(f'No commands made for to save')
##############################################################################################
print(f'NORMAL TERMINATION')