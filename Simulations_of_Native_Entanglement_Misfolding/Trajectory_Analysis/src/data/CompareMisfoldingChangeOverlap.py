#!/usr/bin/env python3
import logging, os, sys
import time
import argparse
import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.colors as mcolors
import seaborn as sns
from scipy.stats import mode, ttest_1samp, ttest_ind
import pickle
#pd.set_option('display.max_rows', 5000)

class Analysis:
    """
    A class to perform analysis on Molecular dynamics trajecotry data and plot the results. 
    """
    #######################################################################################
    def __init__(self, args):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
        ("--outpath", type=str, required=True, help="Path to output directory")
        ("--candidates", type=str, required=True, help="A file containing two columns. The candidate tag and the groupID")
        ("--CollectedOPpath", type=str, required=True, help="path to the CollectAndProcessesOP DATA/ folder")
        ("--outname", type=str, required=True, help="base name for output files")
        """

        # parse the parameters 
        self.candidates = pd.read_csv(args.candidates)
        logging.info(f'candidates:\n{self.candidates}')
        print(f'candidates:\n{self.candidates}')

        self.outpath = args.outpath
        logging.info(f'outpath: {self.outpath}')
        print(f'outpath: {self.outpath}')

        self.outname = args.outname
        logging.info(f'outname: {self.outname}')
        print(f'outname: {self.outname}')

        self.CollectedOPpath = args.CollectedOPpath
        logging.info(f'CollectedOPpath: {self.CollectedOPpath}')
        print(f'CollectedOPpath: {self.CollectedOPpath}')


        ## make Plots dir
        self.plot_path = os.path.join(self.outpath, 'Plots')
        if not os.path.exists(self.plot_path):
            os.makedirs(self.plot_path)
            print(f'Made directory: {self.plot_path}')  
        print(f'plot_path: {self.plot_path}') 
        
        ## make data dir
        self.data_path = os.path.join(self.outpath, 'DATA')
        if not os.path.exists(self.data_path):
            os.makedirs(self.data_path)
            print(f'Made directory: {self.data_path}')
        print(f'data_path: {self.data_path}')    
    #######################################################################################

    #######################################################################################
    def load_OP(self,):
        """
        Load the following files
        1. Collected G, Q, K and native frame reference df
        """
        # 1. Collected G, Q, K and native frame reference df
        OPFiles = glob.glob(os.path.join(self.CollectedOPpath, f'CollectedGQK.csv'))
        print(f'OPFiles: {OPFiles}')
        if len(OPFiles) != 1:
            raise ValueError(f'There should only be 1 G,Q,K collection file named CollectedGQK.csv not {len(OPFiles)}')
        self.OP_df = pd.read_csv(OPFiles[0])
        print(f'OP_df:\n{self.OP_df}')  
    #######################################################################################    

    #######################################################################################
    def CalcChangeOverlap(self, NativeBy='MSS', scope='full'):
        print(f'Calc Change Overlap for each protein and traj using NativeBy {NativeBy} and scope {scope}')
        """
        Calculates for each allowed from the percent of changes that have overlap with any other another nature in the frame

        NativeBy: can be MSM, Ref, or None
            MSM = using the metastable state with the highest <Q> and lowest <G> as the native state. Will examine all frames outside this native meta stable state.
            Ref = using <Qref> - 3*sigma and <Gref> + 3*sigma as thresholds for selecting native frames. Will examine all frames with Q greater than the threshold and G less than. 
            None = No native defined frames. examine all frames. 
        scope is the level at which data was analyzed. all chunks of the trajectory were analyzed (full) or just the last 10% (last10)
        """
        #################################################
        ## Determine the scope of the analysis
        ChangeOverlap_df = {'tag':[], 'setID':[], 'traj':[], 'NativeBy':[], 'Metric':[], 'Label':[], 'Frame':[], 'Ent_cID':[], 'CompEnt_cID':[], 'ChangeNature':[], 'ChangeNatureInfo':[], 'Overlap':[], 'ChangeOverlap':[]}
        if scope in ['full', 'last10']:
            if scope == 'full':
                scope_dict = {1:'0-10%', 2:'10-20%', 3:'20-30%', 4:'30-40%', 5:'40-50%', 6:'50-60%', 7:'60-70%', 8:'70-80%', 9:'80-90%', 10:'90-100%'}
            else:
                scope_dict = {10:'90-100%'}
        else:
            raise ValueError(f'Scope can only be full or last10')

        #####################################################
        ## Get the frame chunks to split the analysis on
        FrameEnd = 26666
        print(f'FrameEnd: {FrameEnd}')
        # Split the array into 10 equal chunks
        frame_chunks = [np.arange(0, FrameEnd + 1)] + np.array_split(np.arange(0, FrameEnd + 1), 10)
        print(f'frame_chunks: {frame_chunks}')

        ########################################################################################
        ## Check if the files are already made. If so warn user and load files     
        ## If not then calculate the TotalChange, TotalLoss, TotalGain, and TotalPure values and save the csv   
        ChangeOverlap_outfile = os.path.join(self.data_path, f'ChangeOverlap_Scope-{scope}_NativeBy{NativeBy}_DAD.csv')

        if os.path.exists(ChangeOverlap_outfile):

            logging.info(f'Files already exists and will be loaded')

            self.ChangeOverlap_df = pd.read_csv(ChangeOverlap_outfile)
            logging.info(f'Loaded: {ChangeOverlap_outfile}')


        else:    
            ########################################################################################
            # Get the ChangeOverlap 
            for gene, pdb, chain, setID, set_name in self.candidates.values:

                tag = f'{gene}_{pdb}_{chain}'

                # Find the EntInfo file
                EntInfofile = os.path.join(self.CollectedOPpath, f'CollectedEntInfo_{tag}.csv')
                print(f'EntInfofile: {EntInfofile}')

                ## Quality check that there is atleast 1 file for each made
                if not os.path.exists(EntInfofile):
                    raise ValueError(f'Could not find EntInfo for this traj {traj}')

                ## Load the entanglement info file for this traj
                EntInfo_df = pd.read_csv(EntInfofile, low_memory=False)
                print(f'EntInfo_df:\n{EntInfo_df}')

                ## Load GQK data for this tag
                loc_GQK_df = self.OP_df[self.OP_df['gene'] == gene]
                print(f'loc_GQK_df:\n{loc_GQK_df}')

                for traj in range(0,50):
                    print(gene, pdb, chain, setID, set_name, tag, traj)

                    ## Step 0: Update PLoss and OddsLoss dataframes with header info for this row            

                    ########################################################################################
                    ## Step 1: get the G. EntInfo, and Q files for the Quenching
                    ## Load the entanglement info file for this traj
                    EntInfo_traj_df = EntInfo_df[EntInfo_df['traj'] == traj]
                    print(f'EntInfo_traj_df:\n{EntInfo_traj_df}')

                    # Get the reference state EntInfo
                    RefEntInfo_df = EntInfo_traj_df[EntInfo_traj_df['Frame'] == -1]
                    #print(f'RefEntInfo_df:\n{RefEntInfo_df}')
                
                    ########################################################################################
                    # Step 2: Get list of misfolded frames for this traj
                    traj_msm = loc_GQK_df[loc_GQK_df['traj'] == traj]
                    #print(traj_msm)
                    if NativeBy in ['Ref', 'MSM']:
                        traj_msm = traj_msm[traj_msm[f'NativeBy{NativeBy}'] == False]
                    #print(f'traj_msm:\n{traj_msm}')

                    misfolded_frames = traj_msm['frame'].values
                    #print(f'misfolded_frames: {misfolded_frames}')

                    EntInfo_traj_df = EntInfo_traj_df[EntInfo_traj_df['Frame'].isin(misfolded_frames)]
                    print(f'Misfolded EntInfo_traj_df:\n{EntInfo_traj_df}')

                    ########################################################################################
                    # Step 3: Get the Entanglement info for those frames and count the unique elements of the NchangeType CchangeType
                    for i,label in scope_dict.items():

                        ## Get the upper and lower bounds of frames to analyze 
                        lower_bound = frame_chunks[i][0]
                        upper_bound = frame_chunks[i][-1]
                        lower_bound, upper_bound = int(lower_bound), int(upper_bound)
                        #print(i, label, lower_bound, upper_bound)
                        
                        ## Filter the DataFrame for rows within this percentile range
                        EntInfo_chunk = EntInfo_traj_df[(EntInfo_traj_df['Frame'] >= lower_bound) & (EntInfo_traj_df['Frame'] < upper_bound)]
                        #print(EntInfo_chunk)
        
                        ## If the chunk is empty then move on
                        if len(EntInfo_chunk) == 0:
                            continue
                            
                        else:
                            
                            ###################################################################################
                            ## Loop through each frame of the chunk
                            for frame, frame_df in EntInfo_chunk.groupby('Frame'):
                                #print(f'\nframe: {frame}:\n{frame_df}')

                                ## For each frame split out the changes in entanglement into a dictionary where the keys are an index for that entanglement and the values are a dictionary with info about it
                                # This separation is independant of terminal tail. So if a give loop has changes in both the N and C terminus then each of those changes gets its own entry. 
                                frame_ents = {}
                                entIDX = 0

                                for cID, c_df in frame_df.groupby('cID'):
                                    #print(c_df)
                                  
                                    cID_seti = set()
                                    cID_setj = set()
                                    cID_Nsetcrossings = set()
                                    cID_Csetcrossings = set()

                                    for rowi, row in c_df.iterrows():
                                        seti = np.arange(row['i'] - 3, row['i'] + 4)
                                        setj = np.arange(row['j'] - 3, row['j'] + 4)

                                        ## check for N terminal change. if loss get crossing from reference. 
                                        Nsetcrossings = np.asarray([])
                                        if 'Loss' in row['NchangeType']:
                                            ref_cross = RefEntInfo_df[(RefEntInfo_df['i'] == row['i']) & (RefEntInfo_df['j'] == row['j'])]
                                            cross = ref_cross['crossingsN'].values[0]                             
                                        elif 'Gain' in row['NchangeType'] or 'Pure' in row['NchangeType']:
                                            cross = row['crossingsN']
                                        else:
                                            cross = np.asarray([])

                                        if isinstance(cross, str):
                                            Nsetcrossings = np.hstack([np.arange(abs(float(c.replace('*',''))) - 3, abs(float(c.replace('*',''))) + 4) for c in cross.split(',')]).astype(int)
                                        elif isinstance(cross, float):
                                            cross = abs(cross)
                                            Nsetcrossings = np.arange(cross - 3, cross + 4)

                                        
                                        ## check for C terminal change. if loss get crossing from reference. 
                                        Csetcrossings = np.asarray([])
                                        if 'Loss' in row['CchangeType']:
                                            ref_cross = RefEntInfo_df[(RefEntInfo_df['i'] == row['i']) & (RefEntInfo_df['j'] == row['j'])]
                                            cross = ref_cross['crossingsC'].values[0]                                        
                                        elif 'Gain' in row['CchangeType'] or 'Pure' in row['CchangeType']:
                                            cross = row['crossingsC']
                                        else:
                                            cross = np.asarray([])

                                        if isinstance(cross, str):
                                            Csetcrossings = np.hstack([np.arange(abs(float(c.replace('*',''))) - 3, abs(float(c.replace('*',''))) + 4) for c in cross.split(',')]).astype(int)
                                        elif isinstance(cross, float):
                                            cross = abs(cross)
                                            Csetcrossings = np.arange(cross - 3, cross + 4)

                                        #print(cID, seti, setj, Nsetcrossings, Csetcrossings)
                                        cID_seti = cID_seti.union(seti)
                                        cID_setj = cID_setj.union(setj)
                                        cID_Nsetcrossings = cID_Nsetcrossings.union(Nsetcrossings)
                                        cID_Csetcrossings = cID_Csetcrossings.union(Csetcrossings)
                                    
                                    if row['NchangeType'] != 'NoChange':
                                        frame_ents[entIDX] = {'loop':cID_seti.union(cID_setj), 'changeType':row['NchangeType'], 'setcrossings':cID_Nsetcrossings, 'Termini':'N'}
                                        entIDX += 1
                                    if row['CchangeType'] != 'NoChange':
                                        frame_ents[entIDX] = {'loop':cID_seti.union(cID_setj), 'changeType':row['CchangeType'], 'setcrossings':cID_Csetcrossings, 'Termini':'C'}
                                        entIDX += 1 
                                ###################################################################################

                                ######################################################################
                                # Calculate the number of unique ent with any overlap with another unique ent of a different nature 
                                #print(f'\n### {tag} TRAJ: {traj} FRAME: {frame}\nFrameDF: {frame_df.to_string()}')
                                for entIDX in frame_ents.keys():                             
                                    entIDX_info = frame_ents[entIDX]
                                    #print(f'\nChecking cID {entIDX} {entIDX_info} for overlap with:')
                                    ChangeOverlap = False

                                    for testID, info in frame_ents.items():
                                        if entIDX == testID:
                                            continue
                                        #print(testID, info)

                                        # Determine if there was a change in Nature
                                        ChangeNature = False
                                        if 'Loss' in entIDX_info['changeType'] and 'Loss' not in info['changeType']:
                                            ChangeNature = True
                                        elif 'Gain' in entIDX_info['changeType'] and 'Gain' not in info['changeType']:
                                            ChangeNature = True
                                        elif 'Pure' in entIDX_info['changeType'] and 'Pure' not in info['changeType']:
                                            ChangeNature = True
                                        #print(f'ChangeNature: {ChangeNature}')
                                        ChangeNature_str = entIDX_info['changeType'] + '->' + info['changeType']

                                        #Determine if there was any overlap
                                        set1 = entIDX_info['loop'].union(entIDX_info['setcrossings'])
                                        #print(f'set1: {set1}')
                                        set2 = info['loop'].union(info['setcrossings'])
                                        #print(f'set2: {set2}')
                                        Overlap = set1.intersection(set2)
                                        Overlap_str = ','.join([str(s) for s in Overlap])
                                        #print(f'Overlap: {Overlap}')
                                        #print(f'{testID} {info}\nChangeNature: {ChangeNature} Overlap: {Overlap}')

                                        if Overlap and ChangeNature:
                                            #print(f'Both Overlap and ChangeNature conditions met: {Overlap} {ChangeNature}')
                                            ChangeOverlap = True

                                        ChangeOverlap_df['tag'] += [tag]
                                        ChangeOverlap_df['setID'] += [setID]
                                        ChangeOverlap_df['traj'] += [traj]
                                        ChangeOverlap_df['NativeBy'] += [NativeBy]
                                        ChangeOverlap_df['Metric'] += ['ChangeOverlap']
                                        ChangeOverlap_df['Label'] += [label]
                                        ChangeOverlap_df['Frame'] += [frame]
                                        ChangeOverlap_df['Ent_cID'] += [entIDX]
                                        ChangeOverlap_df['CompEnt_cID'] += [testID]
                                        ChangeOverlap_df['ChangeNature'] += [ChangeNature]
                                        ChangeOverlap_df['ChangeNatureInfo'] += [ChangeNature_str]
                                        ChangeOverlap_df['Overlap'] += [Overlap_str]
                                        ChangeOverlap_df['ChangeOverlap'] += [ChangeOverlap]
                                ###################################################################################
                            
            ## Save the PLoss, OddsLoss, and component metric dataframes
            for k,v in ChangeOverlap_df.items():
                print(k, len(v))
            self.ChangeOverlap_df = pd.DataFrame(ChangeOverlap_df)           
            print(f'ChangeOverlap_df:\n{self.ChangeOverlap_df}')
            self.ChangeOverlap_df.to_csv(ChangeOverlap_outfile, index=False)
            logging.info(f'SAVED: {ChangeOverlap_outfile}')

        return self.ChangeOverlap_df
    #######################################################################################
    
    #######################################################################################
    def ChangeStats(self, df, setID=2, outfiletag=f'None', OnlyGainLoss=False, label=f'None'):
        """
        Takes the ChangeOverlap dataframe created by CalcChangeOverlap and for each frame calculates the fraction of changes that switched their nature and had overlap.
        """
        print(f'ChangeStats: setID {setID} | OnlyGainLoss: {OnlyGainLoss} | label: {label}')
        outfile_csv = os.path.join(self.data_path, f'ChangeOverlapStats_{outfiletag}_OnlyGainLoss-{OnlyGainLoss}_label-{label}.csv')

        # select only those rows for the requested setID and OnlyGainLoss arguments
        # Select rows where 'ChangeNatureInfo' contains only 'Gain' and or 'Loss'. No pure
        if OnlyGainLoss:
            filtered_df = df[df["ChangeNatureInfo"].str.contains("Gain") | df["ChangeNatureInfo"].str.contains("Loss", na=False)]
            filtered_df = filtered_df[filtered_df['setID']==setID]
            filtered_df = filtered_df[filtered_df['Label']==label]
        else:
            filtered_df = df[df['setID']==setID]
            filtered_df = filtered_df[filtered_df['Label']==label]


        ## loop through the tag groups, then the traj groups, the the frame groups, then finally by the Ent_cID
        Prob = 0
        marProbLoss = 0 # P(L->G)
        numLoss = 0
        marProbGain = 0
        numGain = 0
        marProbPure = 0
        numPure = 0
        NumChanges = 0

        perTrajOverlap = []
        perTrajOverlapCondLoss = []
        perTrajOverlapCondGain = []
        perTrajOverlapCondPure = []

        perFrameOverlap = []
        perFrameOverlapCondLoss = []
        perFrameOverlapCondGain = []
        perFrameOverlapCondPure = []
        start_time = time.time()
        for tag, tag_df in filtered_df.groupby('tag'):
            for traj, traj_df in tag_df.groupby('traj'):
                perTraj_NumChanges = 0
                perTraj_NumOverlap = 0
                perTraj_NumLoss = 0
                perTraj_NumLossOverlap = 0
                perTraj_NumGain = 0
                perTraj_NumGainOverlap = 0
                perTraj_NumPure = 0
                perTraj_NumPureOverlap = 0

                for frame, frame_df in traj_df.groupby('Frame'):
                    #print(f'frame_df:\n{frame_df.to_string()}')
                    unique_Ent_cID = np.unique(frame_df['Ent_cID'].values)
                    #if len(unique_Ent_cID) == 1:
                    #        continue
                    perFrame_NumChanges = 0
                    perFrame_NumOverlap = 0
                    perFrame_NumLoss = 0
                    perFrame_NumLossOverlap = 0
                    perFrame_NumGain = 0
                    perFrame_NumGainOverlap = 0
                    perFrame_NumPure = 0
                    perFrame_NumPureOverlap = 0
                    for Ent_cID, Ent_cID_df in frame_df.groupby('Ent_cID'):
                        perTraj_NumChanges += 1 # keep track of the number of unique entanglements
                        perFrame_NumChanges += 1
                        NumChanges += 1

                        if Ent_cID_df['ChangeOverlap'].any():
                            perTraj_NumOverlap += 1
                            perFrame_NumOverlap += 1
                            Prob += 1
                        
                        if Ent_cID_df['ChangeNatureInfo'].str.startswith('Loss').all():
                            perTraj_NumLoss += 1
                            perFrame_NumLoss += 1
                            numLoss += 1
                        if Ent_cID_df['ChangeNatureInfo'].str.startswith('Gain').all():
                            perTraj_NumGain += 1  
                            perFrame_NumGain += 1  
                            numGain += 1              
                        if Ent_cID_df['ChangeNatureInfo'].str.startswith('Pure').all():
                            perTraj_NumPure += 1    
                            perFrame_NumPure += 1  
                            numPure += 1

                        if Ent_cID_df['ChangeNatureInfo'].str.startswith('Loss').all() and Ent_cID_df['ChangeOverlap'].any():
                            perTraj_NumLossOverlap += 1
                            perFrame_NumLossOverlap += 1
                            marProbLoss += 1
                        
                        if Ent_cID_df['ChangeNatureInfo'].str.startswith('Gain').all() and Ent_cID_df['ChangeOverlap'].any():
                            perTraj_NumGainOverlap += 1
                            perFrame_NumGainOverlap += 1
                            marProbGain += 1

                        if Ent_cID_df['ChangeNatureInfo'].str.startswith('Pure').all() and Ent_cID_df['ChangeOverlap'].any():
                            perTraj_NumPureOverlap += 1
                            perFrame_NumPureOverlap += 1
                            marProbPure += 1

                    perFrameOverlap += [perFrame_NumOverlap/perFrame_NumChanges]   
                    #print(f'Frame: {frame} | perFrameOverlap: {perFrame_NumOverlap/perFrame_NumChanges} | perFrame_NumOverlap {perFrame_NumOverlap} | perFrame_NumChanges: {perFrame_NumChanges}')
 
                    if perFrame_NumLoss != 0:
                        perFrameOverlapCondLoss += [perFrame_NumLossOverlap/perFrame_NumLoss]
                        #print(f'Frame: {frame} | perFrameOverlapCondLoss: {perFrame_NumLossOverlap/perFrame_NumLoss} | perFrame_NumLossOverlap {perFrame_NumLossOverlap} | perFrame_NumLoss: {perFrame_NumLoss}')

                    if perFrame_NumGain != 0:
                        perFrameOverlapCondGain += [perFrame_NumGainOverlap/perFrame_NumGain]
                        #print(f'Frame: {frame} | perFrameOverlapCondGain: {perFrame_NumGainOverlap/perFrame_NumGain} | perFrame_NumGainOverlap {perFrame_NumGainOverlap} | perFrame_NumGain: {perFrame_NumGain}')

                    if perFrame_NumPure != 0:
                        perFrameOverlapCondPure += [perFrame_NumPureOverlap/perFrame_NumPure]
                        #print(f'Frame: {frame} | perFrameOverlapCondPure: {perFrame_NumPureOverlap/perFrame_NumPure} | perFrame_NumPureOverlap {perFrame_NumPureOverlap} | perFrame_NumPure: {perFrame_NumPure}')

                perTrajOverlap += [perTraj_NumOverlap/perTraj_NumChanges]    
                if perTraj_NumLoss != 0:
                    perTrajOverlapCondLoss += [perTraj_NumLossOverlap/perTraj_NumLoss]
                if perTraj_NumGain != 0:
                    perTrajOverlapCondGain += [perTraj_NumGainOverlap/perTraj_NumGain]
                if perTraj_NumPure != 0:
                    perTrajOverlapCondPure += [perTraj_NumPureOverlap/perTraj_NumPure]

        FrameOverlap_stats = get_stats(perFrameOverlap)
        print(f'\nFrameOverlap_stats: {FrameOverlap_stats}')
        FrameOverlapCondLoss_stats = get_stats(perFrameOverlapCondLoss)
        print(f'FrameOverlapCondLoss_stats: {FrameOverlapCondLoss_stats}')
        FrameOverlapCondGain_stats = get_stats(perFrameOverlapCondGain)
        print(f'FrameOverlapCondGain_stats: {FrameOverlapCondGain_stats}')
        FrameOverlapCondPure_stats = get_stats(perFrameOverlapCondPure)
        print(f'FrameOverlapCondPure_stats: {FrameOverlapCondPure_stats}')


        TrajOverlap_stats = get_stats(perTrajOverlap)
        print(f'\nTrajOverlap_stats: {TrajOverlap_stats}')
        TrajOverlapCondLoss_stats = get_stats(perTrajOverlapCondLoss)
        print(f'TrajOverlapCondLoss_stats: {TrajOverlapCondLoss_stats}')
        TrajOverlapCondGain_stats = get_stats(perTrajOverlapCondGain)
        print(f'TrajOverlapCondGain_stats: {TrajOverlapCondGain_stats}')
        TrajOverlapCondPure_stats = get_stats(perTrajOverlapCondPure)
        print(f'TrajOverlapCondPure_stats: {TrajOverlapCondPure_stats}')


        print(f'\nOverall probability and marginals')
        print(f'Prob overlap: {Prob/NumChanges}')
        print(f'Prob overlap given Loss: {marProbLoss/numLoss}')
        print(f'Prob overlap given Gain: {marProbGain/numGain}')
        print(f'Prob overlap given Pure: {marProbPure/numPure}')

        end_time = time.time()
        print(f'delta time: {end_time - start_time}')
    
    #######################################################################################

#######################################################################################
def bootstrap(data):
    boot_means = []
    for b in range(10000):
        boot_samp = np.random.choice(data, size=len(data))
        boot_mean = np.mean(boot_samp)
        #print(b, boot_mean)
        boot_means += [boot_mean]

    lb = np.percentile(boot_means, 2.5)
    ub = np.percentile(boot_means, 97.5)
    return (lb, ub)
#######################################################################################

#######################################################################################
def get_stats(arr):
    """
    Get the <> and 95% ci for a stats array
    """
    mean = np.mean(arr)
    std = np.std(arr)
    lb, ub = bootstrap(arr)
    median = np.median(arr)
    return (median, mean, std, lb, ub)
#######################################################################################

#######################################################################################
def statistic(x, y, axis):
    return np.mean(x, axis=axis) - np.mean(y, axis=axis)
#######################################################################################

#######################################################################################
def permutation_test(arr1, arr2, num_permutations=10000):
    """
    Perform a permutation test to compare two arrays.

    Parameters:
    - arr1: First array of data.
    - arr2: Second array of data.
    - num_permutations: Number of permutations to perform (default is 10,000).

    Returns:
    - p_value: The p-value of the test.
    """
    # Calculate the observed difference in means (or some other statistic you care about)
    observed_diff = np.mean(arr1) - np.mean(arr2)

    # Combine both arrays for random permutation
    combined = np.concatenate([arr1, arr2])

    # Perform permutations and calculate the difference in means for each permutation
    perm_diffs = []
    for _ in range(num_permutations):
        np.random.shuffle(combined)  # Shuffle the combined array
        perm_arr1 = combined[:len(arr1)]  # First array after shuffle
        perm_arr2 = combined[len(arr1):]  # Second array after shuffle
        perm_diff = np.mean(perm_arr1) - np.mean(perm_arr2)  # Difference in means
        perm_diffs.append(perm_diff)

    # Convert perm_diffs to a numpy array for easier calculation
    perm_diffs = np.array(perm_diffs)

    # Calculate p-value: the proportion of permutations where the absolute difference is >= observed difference
    p_value = np.mean(np.abs(perm_diffs) >= np.abs(observed_diff))
    if p_value == 0.0:
        p_value = 1/num_permutations

    return observed_diff, p_value
    #######################################################################################

############## MAIN #################
def main():
    
    script_name = f'CompareMisfoldingMechanism'
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("--candidates", type=str, required=True, help="A file containing two columns. The candidate tag and the groupID")
    parser.add_argument("--CollectedOPpath", type=str, required=True, help="path to the CollectAndProcessesOP DATA/ folder")
    parser.add_argument("--outname", type=str, required=True, help="base name for output files")
    parser.add_argument("--setID", type=int, required=True, help="setID to analyze (2 or 3)")
    args = parser.parse_args()

    ## make output folder
    if not os.path.exists(args.outpath):
        os.makedirs(args.outpath)
        print(f'Made directory: {args.outpath}')

    ## make logging dir
    logs = os.path.join(args.outpath, 'logs')
    if not os.path.exists(logs):
        os.makedirs(logs)
        print(f'Made directory: {logs}')    
    
    ## make DATA dir
    data_path = os.path.join(args.outpath, 'DATA')
    if not os.path.exists(data_path):
        os.makedirs(data_path)
        print(f'Made directory: {data_path}')  
   

    # Setup logging configuration
    logfile = os.path.join(logs, f'{args.outname}.log')
    print(f'logfile: {logfile}')
    logging.basicConfig(filename=logfile, level=logging.INFO, format='%(asctime)s %(message)s')
    logging.info(f'{"#"*50}NEW RUN {script_name}{"#"*50}')

    
    # Step 0: initialize the simulation object 
    anal = Analysis(args)
    anal.load_OP()
    
    # Step 4: Calculate metrics, get stats, and plot
    for NativeBy in ['Ref', 'MSS']:

        ###############################################################################
        # Change Overlap analysis to determine for each frame the fraction of changes that have overlap with a change of a different nature
        ChangeOverlap_df = anal.CalcChangeOverlap(NativeBy=NativeBy, scope='last10')
        print(f'ChangeOverlap_df:\n{ChangeOverlap_df}')
        ChangeStats_df = anal.ChangeStats(ChangeOverlap_df, setID=args.setID, outfiletag=f'setID{args.setID}_NativeBy{NativeBy}', OnlyGainLoss=False, label='90-100%')
        ChangeStats_df = anal.ChangeStats(ChangeOverlap_df, setID=args.setID, outfiletag=f'setID{args.setID}_NativeBy{NativeBy}', OnlyGainLoss=True, label='90-100%')
       
        ###############################################################################


    print(f'logfile: {logfile}')

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()

print(f'NORMAL TERMINATION: {time.time() - start_time}')
logging.info(f'NORMAL TERMINATION: {time.time() - start_time}')