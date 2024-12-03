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
from scipy.stats import mode, permutation_test, ttest_1samp, ttest_ind, mannwhitneyu
import pickle
from multiprocessing import Pool, cpu_count
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
    def CalcChangeMetrics(self, NativeBy='MSS', Metric='None', scope='full'):
        print(f'Calc {Metric} for each protein and traj')
        """
        Calculates the requested metric over the requested scope of the trajectory data provided

        NativeBy: can be MSS, Ref, or None
            MSS = using the metastable state with the highest <Q> and lowest <G> as the native state. Will examine all frames outside this native meta stable state.
            Ref = using <Qref> - 3*sigma and <Gref> + 3*sigma as thresholds for selecting native frames. Will examine all frames with Q greater than the threshold and G less than. 
            None = No native defined frames. examine all frames. 
        Metric: can be PLoss, PGain, PPure, OddsLoss, OddsGain, OddsPure, PLossNoPure, OddsLossNoPure, PGainNoPure, OddsGainNoPure
        scope is the level at which data was analyzed. all chunks of the trajectory were analyzed (full) or just the last 10% (last10)
        """

        # Quality check to ensure the Metric requested is in the list
        if Metric not in 'PLoss, PGain, PPure, OddsLoss, OddsGain, OddsPure, PLossNoPure, OddsLossNoPure, PGainNoPure, OddsGainNoPure':
            raise ValueError(f'The Metric {Metric} is not in the allowed list: PLoss, PGain, PPure, OddsLoss, OddsGain, OddsPure, PLossNoPure, OddsLossNoPure, PGainNoPure, OddsGainNoPure')

        #################################################
        ## Determine the scope of the analysis
        if scope in ['full', 'last10']:
            if scope == 'full':
                scope_dict = {0:'full', 1:'0-10%', 2:'10-20%', 3:'20-30%', 4:'30-40%', 5:'40-50%', 6:'50-60%', 7:'60-70%', 8:'70-80%', 9:'80-90%', 10:'90-100%'}
                Change_df = {'tag':[], 'setID':[], 'traj':[], 'NativeBy':[], 'Metric':[],
                    'full':[], '0-10%':[], '10-20%':[], '20-30%':[], '30-40%':[], 
                    '40-50%':[], '50-60%':[], '60-70%':[], '70-80%':[], '80-90%':[], '90-100%':[]}
                df = {k:[] for k,v in Change_df.items()}

            elif scope == 'last10':
                scope_dict = {0:'full', 10:'90-100%'}
                Change_df = {'tag':[], 'setID':[], 'traj':[], 'NativeBy':[], 'Metric':[],
                    'full':[], '90-100%':[]}
                df = {k:[] for k,v in Change_df.items()}
        else:
            raise ValueError(f'Scope can only be full or last10')

        #####################################################
        ## Get the frame chunks to split the analysis on
        FrameEnd = 26666
        #print(f'FrameEnd: {FrameEnd}')
        # Split the array into 10 equal chunks
        frame_chunks = [np.arange(0, FrameEnd + 1)] + np.array_split(np.arange(0, FrameEnd + 1), 10)
        #print(f'frame_chunks: {frame_chunks}')

        ########################################################################################
        ## Check if the files are already made. If so warn user and load files     
        ## If not then calculate the TotalChange, TotalLoss, TotalGain, and TotalPure values and save the csv   
        TotalChange_outfile = os.path.join(self.data_path, f'TotalChange_Scope-{scope}_NativeBy{NativeBy}.csv')

        if os.path.exists(TotalChange_outfile):

            logging.info(f'Files already exists and will be loaded')

        else:    
            ########################################################################################
            # Get the TotalChange, TotalLoss, TotalGain, and TotalPure for each traj and save it
            for gene, pdb, chain, setID, set_name in self.candidates.values:
                tag = f'{gene}_{pdb}_{chain}'
                print(tag, gene, pdb, chain)

                ########################################################################################
                ## Step 1: get the G. EntInfo, and Q files for the Quenching
                #'Q-NativeFiles':[], 'Q-QuenchFiles':[], 'G-NativeFiles':[], 'EntInfo-NativeFiles':[], 'G-QuenchFiles':[], 'EntInfo-QuenchFiles':[]
                EntInfofile = os.path.join(self.CollectedOPpath, f'CollectedEntInfo_{tag}.csv')
                print(f'EntInfofile: {EntInfofile}')

                ## Quality check that there is atleast 1 file for each made
                if not os.path.exists(EntInfofile):
                    raise ValueError(f'Could not find EntInfo for this traj {traj}')

                ## Load the entanglement info file for this traj
                EntInfo_df = pd.read_csv(EntInfofile, low_memory=False)
                #print(f'EntInfo_df:\n{EntInfo_df}')

                ## Load GQK data for this tag
                loc_GQK_df = self.OP_df[self.OP_df['gene'] == gene]
                #print(f'loc_GQK_df:\n{loc_GQK_df}')

                for traj in range(0,50):
                    print(gene, pdb, chain, setID, set_name, tag, traj)

                    ## Step 0: Update PLoss and OddsLoss dataframes with header info for this row
                    Change_df['tag'] += [tag]
                    Change_df['setID'] += [setID]
                    Change_df['traj'] += [traj]
                    Change_df['NativeBy'] += [NativeBy]
                    Change_df['Metric'] += ['TotalChange']

                    Change_df['tag'] += [tag]
                    Change_df['setID'] += [setID]
                    Change_df['traj'] += [traj]
                    Change_df['NativeBy'] += [NativeBy]
                    Change_df['Metric'] += ['TotalGain']

                    Change_df['tag'] += [tag]
                    Change_df['setID'] += [setID]
                    Change_df['traj'] += [traj]
                    Change_df['NativeBy'] += [NativeBy]
                    Change_df['Metric'] += ['TotalLoss']

                    Change_df['tag'] += [tag]
                    Change_df['setID'] += [setID]
                    Change_df['traj'] += [traj]
                    Change_df['NativeBy'] += [NativeBy]   
                    Change_df['Metric'] += ['TotalPure']               

                    ## Load the entanglement info file for this traj
                    EntInfo_traj_df = EntInfo_df[EntInfo_df['traj'] == traj]
                    #print(f'EntInfo_traj_df:\n{EntInfo_traj_df}')

                    ########################################################################################
                    # Step 2: Get list of misfolded frames for this traj
                    loc_GQK_traj_df = loc_GQK_df[loc_GQK_df['traj'] == traj]
                    if NativeBy in ['Ref', 'MSS']:
                        loc_GQK_traj_df = loc_GQK_traj_df[loc_GQK_traj_df[f'NativeBy{NativeBy}'] == False]
                    #print(f'loc_GQK_traj_df:\n{loc_GQK_traj_df}')

                    misfolded_frames = loc_GQK_traj_df['frame'].values
                    #print(f'misfolded_frames: {misfolded_frames}')

                    EntInfo_traj_df = EntInfo_traj_df[EntInfo_traj_df['Frame'].isin(misfolded_frames)]
                    #print(f'Misfolded EntInfo_traj_df:\n{EntInfo_traj_df}')

                    ########################################################################################
                    # Step 3: Get the Entanglement info for those frames and count the unique elements of the NchangeType CchangeType
                    for i,label in scope_dict.items():
                        # Calculate the percentile range
                        lower_bound = frame_chunks[i][0]
                        upper_bound = frame_chunks[i][-1]
                        lower_bound, upper_bound = int(lower_bound), int(upper_bound)
                        #print(i, label, lower_bound, upper_bound)
                        
                        # Filter the DataFrame for rows within this percentile range
                        EntInfo_chunk = EntInfo_traj_df[(EntInfo_traj_df['Frame'] >= lower_bound) & (EntInfo_traj_df['Frame'] < upper_bound)]
                        #print(EntInfo_chunk)
                   
                        if len(EntInfo_chunk) == 0:
                            TotalChange = 0
                            TotalGain = 0
                            TotalLoss = 0
                            TotalPure = 0

                        else:
                            EntInfo_chunk = EntInfo_chunk[['Time(ns)', 'Frame', 'Gn', 'Gc', 'NchangeType', 'CchangeType', 'len_crossingsN', 'len_crossingsC', 'cID']].drop_duplicates()

                            TotalChange = 0
                            TotalGain = 0
                            TotalLoss = 0
                            TotalPure = 0

                            Nterm_counts = EntInfo_chunk['NchangeType'].value_counts()
                            #print(f'Nterm_counts:\n{Nterm_counts}')
                            for key, value in Nterm_counts.items():
                                if key != 'NoChange':
                                    TotalChange += value
                                    if 'Gain' in key:
                                        TotalGain += value
                                    if 'Loss' in key:
                                        TotalLoss += value
                                    if 'Pure' in key:
                                        TotalPure += value

                            Cterm_counts = EntInfo_chunk['CchangeType'].value_counts()
                            #print(f'Cterm_counts:\n{Cterm_counts}')
                            for key, value in Cterm_counts.items():
                                if key != 'NoChange':
                                    TotalChange += value
                                    if 'Gain' in key:
                                        TotalGain += value
                                    if 'Loss' in key:
                                        TotalLoss += value
                                    if 'Pure' in key:
                                        TotalPure += value

                        Change_df[f'{label}'] += [TotalChange]
                        Change_df[f'{label}'] += [TotalGain]
                        Change_df[f'{label}'] += [TotalLoss]
                        Change_df[f'{label}'] += [TotalPure]
            

            ## Save the PLoss, OddsLoss, and component metric dataframes
            Change_df = pd.DataFrame(Change_df)           
            print(f'Change_df:\n{Change_df}')
            Change_df.to_csv(TotalChange_outfile, index=False)
            logging.info(f'SAVED: {TotalChange_outfile}')
        del Change_df

        # Load the saved Change_df
        self.Change_df = pd.read_csv(TotalChange_outfile)
        logging.info(f'Loaded: {TotalChange_outfile}')
        print(f'Loaded: {TotalChange_outfile}')

        ##############################################################
        ## Calculate the actual metrics
        # PLoss, PGain, PPure, OddsLoss, OddsGain, OddsPure, PLossNoPure, OddsLossNoPure, PGainNoPure, OddsGainNoPure
        print(f'Calculating metric: {Metric}')
        for gene, pdb, chain, setID, set_name in self.candidates.values:

            tag = f'{gene}_{pdb}_{chain}'
            # Only look at set 2 for the moment as set 3 is still quenching

            for traj in range(0,50):

                df['tag'] += [tag]
                df['setID'] += [setID]
                df['traj'] += [traj]
                df['NativeBy'] += [NativeBy]   
                df['Metric'] += [Metric] 

                for i,label in scope_dict.items():
                    #print(gene, pdb, chain, setID, set_name, tag, traj, label)
                    #print(self.Change_df[(self.Change_df['tag'] == tag) & (self.Change_df['setID'] == setID) & (self.Change_df['traj'] == traj) & (self.Change_df['Metric'] == 'TotalChange')])
                    TotalChange = self.Change_df[(self.Change_df['tag'] == tag) & (self.Change_df['setID'] == setID) & (self.Change_df['traj'] == traj) & (self.Change_df['Metric'] == 'TotalChange')][label].values[0]
                    TotalLoss = self.Change_df[(self.Change_df['tag'] == tag) & (self.Change_df['setID'] == setID) & (self.Change_df['traj'] == traj) & (self.Change_df['Metric'] == 'TotalLoss')][label].values[0]
                    TotalGain = self.Change_df[(self.Change_df['tag'] == tag) & (self.Change_df['setID'] == setID) & (self.Change_df['traj'] == traj) & (self.Change_df['Metric'] == 'TotalGain')][label].values[0]
                    TotalPure = self.Change_df[(self.Change_df['tag'] == tag) & (self.Change_df['setID'] == setID) & (self.Change_df['traj'] == traj) & (self.Change_df['Metric'] == 'TotalPure')][label].values[0]
                    #print(TotalChange, TotalLoss, TotalGain, TotalPure)
             
                    if TotalChange == 0:
                        M = -1

                    elif Metric == 'PLoss':
                        M = TotalLoss/TotalChange
                        if M > 1:
                            raise ValueError(f'PLoss > 1: {M}')

                    elif Metric == 'PGain':
                        M = TotalGain/TotalChange

                    elif Metric == 'PPure':
                        M = TotalPure/TotalChange

                    elif Metric == 'PLossNoPure':
                        M = TotalLoss/(TotalLoss + TotalGain)

                    elif Metric == 'PGainNoPure':
                        M = TotalGain/(TotalLoss + TotalGain)

                    elif Metric == 'OddsLoss':
                        P = TotalLoss/TotalChange
                        if P != 1:
                            M = (P) / (1 - P)
                        else:
                            M = np.inf

                    elif Metric == 'OddsGain':
                        P = TotalGain/TotalChange
                        if P != 1:
                            M = (P) / (1 - P)
                        else:
                            M = np.inf

                    elif Metric == 'OddsPure':
                        P = TotalPure/TotalChange
                        if P != 1:
                            M = (P) / (1 - P)
                        else:
                            M = np.inf

                    elif Metric == 'OddsLossNoPure':
                        P = TotalLoss/(TotalLoss + TotalGain)
                        if P != 1:
                            M = (P) / (1 - P)
                        else:
                            M = np.inf

                    elif Metric == 'OddsGainNoPure':
                        P = TotalGain/(TotalLoss + TotalGain)
                        if P != 1:
                            M = (P) / (1 - P)
                        else:
                            M = np.inf

                    df[f'{label}'] += [M]

        df = pd.DataFrame(df)           
        print(f'df:\n{df}')
        df_outfile = os.path.join(self.data_path, f'{Metric}_Scope-{scope}_NativeBy{NativeBy}.csv')
        df.to_csv(df_outfile, index=False)
        logging.info(f'SAVED: {df_outfile}')
        print(f'SAVED: {df_outfile}')
        return df
    #######################################################################################

    #######################################################################################
    def OneSampleStats(self, df, setID=2, NativeBy='MSS', test_mean=1, label='None', scope='full'):
        """
        Calculate the OneSampleStats for any metric: Avg, Median, std, 95% ci of the mean, 1sample ttest stat and pvalue
        setID is an integer to be used to select the data
        df is a metric dataframe 
        NativeBy: can be MSS, Ref, or None
            MSS = using the metastable state with the highest <Q> and lowest <G> as the native state. Will examine all frames outside this native meta stable state.
            Ref = using <Qref> - 3*sigma and <Gref> + 3*sigma as thresholds for selecting native frames. Will examine all frames with Q greater than the threshold and G less than. 
            None = No native defined frames. examine all frames. 
        test_mean is the population mean used in the 1sample ttest
        outfiletag is just a string appended to the image file name
        scope is the level at which data was analyzed. all chunks of the trajectory were analyzed (full) or just the last 10% (last10)
        """
        logging.info(f'Performing 1 sample T-Test on {label} data for setID: {setID} and NativeBy: {NativeBy}')
        print(f'Performing 1 sample T-Test on {label} data for setID: {setID} and NativeBy: {NativeBy}')
        df = df[(df['setID'] == setID) & (df['NativeBy']==NativeBy)]
        print(f'{label}_df:\n{df}')
    
        #################################################
        ## Determine the scope of the analysis
        if scope in ['full', 'last10']:
            if scope == 'full':
                scope_list = [f'full', f'0-10%', f'10-20%', f'20-30%', f'30-40%', f'40-50%', f'50-60%', f'60-70%', f'70-80%', f'80-90%', f'90-100%']
            else:
                scope_list = [f'90-100%']

        else:
            raise ValueError(f'Scope can only be full or last10')

        ########################################################################################
        ## for each protein calculate the <R> and confidence intervals 
        df_stats = {'setID':[], 'NativeBy':[], 'Label':[], 'Metric':[], 'median':[], 'avg':[], 'std':[], 'lb':[], 'ub':[], 'ttest':[], 'pvalue':[]}

        for key in scope_list:
            df_stats['setID'] += [setID]
            df_stats['NativeBy'] += [NativeBy]
            df_stats['Label'] += [key]
            df_stats['Metric'] += [label]
            
            #drop -1 values and replace np.inf values with the largest OddsLoss
            data = df[key].values
            data = data[np.where(data != -1)]
            max_finite_value = np.max(data[np.isfinite(data)])
            data[data == np.inf] = max_finite_value
            #print(data, max_finite_value)

            median, mean, std, lb, ub = get_stats(data)
            # Perform the one-sample t-test
            if test_mean == None:
                t_statistic, p_value = np.nan, np.nan
            else:
                t_statistic, p_value = ttest_1samp(data, test_mean)

            #print(f'setID {setID} | NativeBy {NativeBy} | {key} | <{label}> = {mean} +/- ({lb}, {ub}) | simga({label}) = {std} | t_statistic: {t_statistic} w/ p_value: {p_value}')
            df_stats[f'median'] += [median]
            df_stats[f'avg'] += [mean]
            df_stats[f'std'] += [std]
            df_stats[f'lb'] += [lb]
            df_stats[f'ub'] += [ub]
            df_stats[f'ttest'] += [t_statistic]
            df_stats[f'pvalue'] += [p_value]

        df_stats = pd.DataFrame(df_stats)
        print(f'{label} df_stats:\n{df_stats}')
        outfile = os.path.join(self.data_path, f'{label}_setID{setID}_OneSampleStats_{scope}_NativeBy{NativeBy}.csv')
        df_stats.to_csv(outfile, index=False)
        logging.info(f'SAVED: {outfile}') 
        print(f'SAVED: {outfile}') 
        return df_stats   
    #######################################################################################

    #######################################################################################
    def TwoSampleStats(self, df1, df2, NativeBy='MSS', setID=2, test='permutation', labels=['df1', 'df2'], scope='full'):
        """
        Calculate the TwoSampleStats for any metric between two sets of data using either permutation or ttest
        setIDs is a tuple of two setIDs to be used to select the data
        df is a metric dataframe 
        NativeBy: can be MSS, Ref, or None
            MSS = using the metastable state with the highest <Q> and lowest <G> as the native state. Will examine all frames outside this native meta stable state.
            Ref = using <Qref> - 3*sigma and <Gref> + 3*sigma as thresholds for selecting native frames. Will examine all frames with Q greater than the threshold and G less than. 
            None = No native defined frames. examine all frames. 
        test is permutation or ttest depending on which test you want to perform
        outfiletag is just a string appended to the image file name
        scope is the level at which data was analyzed. all chunks of the trajectory were analyzed (full) or just the last 10% (last10)
        """
        print(f'Performing {test} on {labels} data for and NativeBy: {NativeBy}')
        df1 = df1[(df1['setID'].isin([setID])) & (df1['NativeBy']==NativeBy)]
        df2 = df2[(df2['setID'].isin([setID])) & (df2['NativeBy']==NativeBy)]
        print(f'{labels[0]}_df1:\n{df1}')
        print(f'{labels[1]}_df2:\n{df2}')

        #################################################
        ## Determine the scope of the analysis
        if scope in ['full', 'last10']:
            if scope == 'full':
                scope_list = [f'full', f'0-10%', f'10-20%', f'20-30%', f'30-40%', f'40-50%', f'50-60%', f'60-70%', f'70-80%', f'80-90%', f'90-100%']
            else:
                scope_list = [f'90-100%']

        else:
            raise ValueError(f'Scope can only be full or last10')

        ########################################################################################
        ## for each protein calculate the <R> and confidence intervals 
        df_stats = {'setID':[], 'NativeBy':[], 'Label':[], 'test':[], 'statistic':[], 'pvalue':[]}

        for key in scope_list:

            data1 = df1[key].values
            data1 = data1[np.where(data1 != -1)]
            max_finite_value = np.max(data1[np.isfinite(data1)])
            data1[data1 == np.inf] = max_finite_value
            median1, mean1, std1, lb1, ub1 = get_stats(data1)
            stats_str = f'{labels[0]} stats\n<>: {mean1:.3} ({lb1:.3}, {ub1:.3})\nMedian: {median1:.3}\nSigma: {std1:.3}\n'

            data2 = df2[key].values
            data2 = data2[np.where(data2 != -1)]
            max_finite_value = np.max(data2[np.isfinite(data2)])
            data2[data2 == np.inf] = max_finite_value
            median2, mean2, std2, lb2, ub2 = get_stats(data2)
            stats_str += f'\n{labels[1]} stats\n<>: {mean2:.3f} ({lb2:.3}, {ub2:.3})\nMedian: {median2:.3}\nSigma: {std2:.3}\n'

            combined = np.hstack([data1, data2])
            
            if test == 'permutation':
                if mean1 > mean2:
                    res = permutation_test((data1, data2), statistic, vectorized=True, alternative='greater')
                    U1, p = mannwhitneyu(data1, data2, alternative='greater')
                elif mean1 < mean2:
                    res = permutation_test((data1, data2), statistic, vectorized=True, alternative='less')
                    U1, p = mannwhitneyu(data1, data2, alternative='less')
                elif mean1 == mean2:
                    res = permutation_test((data1, data2), statistic, vectorized=True, alternative='two-sided')
                    U1, p = mannwhitneyu(data1, data2, alternative='two-sided')

            elif test == 'ttest':
                res = ttest_ind(data1, data2, alternative='greater', equal_var=False)
            statistic_val = res.statistic
            pvalue = res.pvalue

            ## make stats string
            stats_str += f'\npermut stat: {statistic_val:.3f} | pvalue: {pvalue:.3e}'

            ## MWU test
            stats_str += f'\nMWU stat: {U1:.3f} | pvalue: {p:.3e}'
            print(stats_str)

            df_stats['setID'] += [setID]
            df_stats['NativeBy'] += [NativeBy]
            df_stats['Label'] += [key]
            df_stats['test'] += [test]
            df_stats['statistic'] += [statistic_val]
            df_stats['pvalue'] += [pvalue]

            ## Plot histogram
            plt.hist(data1, bins=30, range=(min(combined), max(combined)), label=labels[0], density=True, alpha=0.5)
            plt.hist(data2, bins=30, range=(min(combined), max(combined)), label=labels[1], density=True, alpha=0.5)

            # Add custom text below the legend
            legend = plt.legend(loc='best')

            #plt.text(0.1, 0.75, stats_str, transform=plt.gca().transAxes)
            plt.text(1.05, 1.0, stats_str, transform=plt.gca().transAxes, fontsize=8, va='top', ha='left')
             # Add text to the last line of the legend
            plt.xlabel(f'Prob(ChangeType)')
            plt.ylabel(f'PDF')
            plt.tight_layout()
            outfile = os.path.join(self.plot_path, f'{labels[0]}Vs{labels[1]}_{test}_TwoSampleStats_Scope-{key}_NativeBy{NativeBy}.png')
            plt.savefig(outfile)
            plt.close()
            plt.clf()
            print(f'SAVED: {outfile}')

        df_stats = pd.DataFrame(df_stats)
        print(df_stats)
        outfile = os.path.join(self.data_path, f'{labels[0]}Vs{labels[1]}_{test}_TwoSampleStats_Scope-{scope}_NativeBy{NativeBy}.csv')
        df_stats.to_csv(outfile, index=False)
        logging.info(f'SAVED: {outfile}')  
        print(f'SAVED: {outfile}')  

        return df_stats   
    #######################################################################################

    #######################################################################################
    def Plot_OneSampleStats(self, df, outfiletag='temp', scope='full'):
        """
        Plots the two sample stats for two sets of data.
        df is a OneSampleStats dataframe with the mean, median, std, and conf. intervals 
        outfiletag is just a string appended to the image file name
        scope is the level at which data was analyzed. all chunks of the trajectory were analyzed (full) or just the last 10% (last10)
        """     
        #print(f'Plotting:\n{df}')
        outfile = os.path.join(self.plot_path, f'{outfiletag}_Stats_Scope-{scope}_plot.png')
        outfile_csv = os.path.join(self.plot_path, f'{outfiletag}_Stats_Scope-{scope}_plot.csv')

        #################################################
        ## Determine the scope of the analysis
        if scope in ['full', 'last10']:
            if scope == 'full':
                fig, axes = plt.subplots(3, 1, figsize=(4, 6), sharex=True)
            else:
                print(f'SCOPE: {scope}')
                fig, axes = plt.subplots(3, 1, figsize=(2, 6), sharex=True)
        else:
            raise ValueError(f'Scope can only be full or last10')
        
        # First plot: avg with error bars (lb and ub)
        axes[0].errorbar(df['Label'], df['avg'], yerr=[df['avg'] - df['lb'], df['ub'] - df['avg']], 
                        fmt='o', capsize=5, label='Average')
        axes[0].set_ylabel('Average (avg)')
        if 'Odds' not in outfiletag:
            #axes[0].set_ylim(bottom=0)
            axes[0].axhline(0.5)
        else:
            axes[0].axhline(1)
        #axes[0].set_title('Average with Error Bars')
        
        # Second plot: median
        axes[1].plot(df['Label'], df['median'], 'o', color='orange', label='Median')
        axes[1].set_ylabel('Median')
        #axes[1].set_title('Median')

        # Third plot: pvalue with log scale
        axes[2].plot(df['Label'], df['pvalue'], 'o', color='red', label='P-Value')
        axes[2].set_yscale('log')
        axes[2].set_ylabel('P-Value (log scale)')
        #axes[2].set_title('P-Value (Log Scale)')
        axes[2].axhline(0.05)
        
        # Set shared x-axis labels
        plt.xticks(rotation=45)
        #axes[2].set_xlabel('Label')
        # Adjust alignment of x-tick labels to be anchored at the upper right
        for tick in axes[2].get_xticklabels():
            tick.set_ha('right')  # Set horizontal alignment to right
            tick.set_va('top')    # Set vertical alignment to top
        
        # Layout adjustments
        plt.suptitle(outfiletag, fontsize=10)
        # Get current xlim and add 0.5 to either side
        current_xlim = plt.gca().get_xlim()
        new_xlim = (current_xlim[0] - 0.5, current_xlim[1] + 0.5)
        plt.xlim(new_xlim)
        plt.tight_layout()
        plt.savefig(outfile)
        plt.close()
        plt.clf()
        logging.info(f'SAVED: {outfile}')
        print(f'SAVED: {outfile}')

        df.to_csv(outfile_csv, index=False)
        print(f'SAVED: {outfile_csv}')
    #######################################################################################

    #######################################################################################
    def ProbObserPureSwitch(self, NativeBy='MSS', scope='full'):
        print(f'Calc ProbObserPureSwitch for each protein and traj')
        """
        Calculates the requested metric over the requested scope of the trajectory data provided

        NativeBy: can be MSS, Ref, or None
            MSS = using the metastable state with the highest <Q> and lowest <G> as the native state. Will examine all frames outside this native meta stable state.
            Ref = using <Qref> - 3*sigma and <Gref> + 3*sigma as thresholds for selecting native frames. Will examine all frames with Q greater than the threshold and G less than. 
            None = No native defined frames. examine all frames. 
        Metric: can be PLoss, PGain, PPure, OddsLoss, OddsGain, OddsPure, PLossNoPure, OddsLossNoPure, PGainNoPure, OddsGainNoPure
        scope is the level at which data was analyzed. all chunks of the trajectory were analyzed (full) or just the last 10% (last10)
        """

        ## Determine the scope of the analysis
        if scope in ['full', 'last10']:
            if scope == 'full':
                scope_dict = {0:'full', 1:'0-10%', 2:'10-20%', 3:'20-30%', 4:'30-40%', 5:'40-50%', 6:'50-60%', 7:'60-70%', 8:'70-80%', 9:'80-90%', 10:'90-100%'}

            elif scope == 'last10':
                scope_dict = {0:'full', 10:'90-100%'}
        else:
            raise ValueError(f'Scope can only be full or last10')

        #####################################################
        ## Get the frame chunks to split the analysis on
        FrameEnd = 26666
        #print(f'FrameEnd: {FrameEnd}')
        # Split the array into 10 equal chunks
        frame_chunks = [np.arange(0, FrameEnd + 1)] + np.array_split(np.arange(0, FrameEnd + 1), 10)
        #print(f'frame_chunks: {frame_chunks}')


        ########################################################################################
        # Get the TotalChange, TotalLoss, TotalGain, and TotalPure for each traj and save it
        num_frames = 0
        num_pureObs = 0
        for gene, pdb, chain, setID, set_name in self.candidates.values:
            tag = f'{gene}_{pdb}_{chain}'
            #print(tag, gene, pdb, chain)

            ########################################################################################
            ## Step 1: get the G. EntInfo, and Q files for the Quenching
            #'Q-NativeFiles':[], 'Q-QuenchFiles':[], 'G-NativeFiles':[], 'EntInfo-NativeFiles':[], 'G-QuenchFiles':[], 'EntInfo-QuenchFiles':[]
            EntInfofile = os.path.join(self.CollectedOPpath, f'CollectedEntInfo_{tag}.csv')
            #print(f'EntInfofile: {EntInfofile}')

            ## Quality check that there is atleast 1 file for each made
            if not os.path.exists(EntInfofile):
                raise ValueError(f'Could not find EntInfo for this traj {traj}')

            ## Load the entanglement info file for this traj
            EntInfo_df = pd.read_csv(EntInfofile, low_memory=False)
            #print(f'EntInfo_df:\n{EntInfo_df}')

            ## Load GQK data for this tag
            loc_GQK_df = self.OP_df[self.OP_df['gene'] == gene]
            #print(f'loc_GQK_df:\n{loc_GQK_df}')

            for traj in range(0,50):
                #print(gene, pdb, chain, setID, set_name, tag, traj)      

                ## Load the entanglement info file for this traj
                EntInfo_traj_df = EntInfo_df[EntInfo_df['traj'] == traj]
                #print(f'EntInfo_traj_df:\n{EntInfo_traj_df}')

                ########################################################################################
                # Step 2: Get list of misfolded frames for this traj
                loc_GQK_traj_df = loc_GQK_df[loc_GQK_df['traj'] == traj]
                if NativeBy in ['Ref', 'MSS']:
                    loc_GQK_traj_df = loc_GQK_traj_df[loc_GQK_traj_df[f'NativeBy{NativeBy}'] == False]
                #print(f'loc_GQK_traj_df:\n{loc_GQK_traj_df}')

                misfolded_frames = loc_GQK_traj_df['frame'].values
                #print(f'misfolded_frames: {misfolded_frames}')

                EntInfo_traj_df = EntInfo_traj_df[EntInfo_traj_df['Frame'].isin(misfolded_frames)]
                #print(f'Misfolded EntInfo_traj_df:\n{EntInfo_traj_df}')

                ########################################################################################
                # Step 3: Get the Entanglement info for those frames and count the unique elements of the NchangeType CchangeType
                for i,label in scope_dict.items():
                    # Calculate the percentile range
                    lower_bound = frame_chunks[i][0]
                    upper_bound = frame_chunks[i][-1]
                    lower_bound, upper_bound = int(lower_bound), int(upper_bound)
                    #print(i, label, lower_bound, upper_bound)
                    
                    # Filter the DataFrame for rows within this percentile range
                    EntInfo_chunk = EntInfo_traj_df[(EntInfo_traj_df['Frame'] >= lower_bound) & (EntInfo_traj_df['Frame'] < upper_bound)]
                    #print(EntInfo_chunk)

                    for frame, frame_df in EntInfo_chunk.groupby('Frame'):
                        #print(frame_df)
                        num_frames += 1 # counter to keep track of all frames analyzed

                        if frame_df['NchangeType'].str.contains('Pure').any() or frame_df['CchangeType'].str.contains('Pure').any():
                            num_pureObs += 1
                            

        P = num_pureObs/num_frames
        print(f'Prob(Pure) = {P} = {num_pureObs}/{num_frames}')



    #######################################################################################

#######################################################################################
def bootstrap(data, num_bootstraps=10000):
    """
    Computes bootstrap confidence intervals using multiprocessing.
    
    Parameters:
        data (array-like): The dataset to bootstrap.
        num_bootstraps (int): Number of bootstrap samples to generate.
        
    Returns:
        lb (float): Lower bound of the 95% confidence interval.
        ub (float): Upper bound of the 95% confidence interval.
    """
    size = len(data)
    num_workers = cpu_count()  # Get the number of CPU cores
    print(f'num_workers: {num_workers}')
    # Create a multiprocessing pool
    with Pool(num_workers) as pool:
        # Map the bootstrap_worker function to the pool
        boot_means = pool.starmap(bootstrap_worker, [(data, size)] * num_bootstraps)

    # Calculate confidence intervals
    lb = np.percentile(boot_means, 2.5)
    ub = np.percentile(boot_means, 97.5)
    
    return lb, ub

def bootstrap_worker(data, size):
    """
    Worker function to compute a single bootstrap mean.
    """
    boot_samp = np.random.choice(data, size=size)
    return np.mean(boot_samp)
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

    setID = args.setID
    # Step 0: initialize the simulation object 
    anal = Analysis(args)
    anal.load_OP()

    # Step 5: Calculate metrics, get stats, and plot
    for NativeBy in ['Ref', 'MSS']:

        ###############################################################################
        # Probability of observing a pure switch
        #anal.ProbObserPureSwitch(NativeBy=NativeBy, scope='last10')  
        #quit()

        ############################################################################### 
        # probability of Loss considering loss, gain, and pure switch changes in ent
        PLoss_df = anal.CalcChangeMetrics(NativeBy=NativeBy, Metric='PLoss', scope='last10')
        print(f'PLoss_df:\n{PLoss_df}')
        PLoss_stats = anal.OneSampleStats(PLoss_df, setID=setID, NativeBy=NativeBy, test_mean=0.5, label='PLoss', scope='last10')
        anal.Plot_OneSampleStats(PLoss_stats, outfiletag=f'PLoss_NativeBy{NativeBy}_{args.outname}', scope='last10')
       
        # probability of Gain considering loss, gain, and pure switch changes in ent
        PGain_df = anal.CalcChangeMetrics(NativeBy=NativeBy, Metric='PGain', scope='last10')
        print(f'PGain_df:\n{PGain_df}')
        PGain_stats = anal.OneSampleStats(PGain_df, setID=setID, NativeBy=NativeBy, test_mean=0.5, label='PGain', scope='last10')
        anal.Plot_OneSampleStats(PGain_stats, outfiletag=f'PGain_NativeBy{NativeBy}_{args.outname}', scope='last10')

        # probability of Pure considering loss, gain, and pure switch changes in ent
        PPure_df = anal.CalcChangeMetrics(NativeBy=NativeBy, Metric='PPure', scope='last10')
        print(f'PPure_df:\n{PPure_df}')
        PPure_stats = anal.OneSampleStats(PPure_df, setID=setID, NativeBy=NativeBy, test_mean=0.5, label='PPure', scope='last10')
        anal.Plot_OneSampleStats(PPure_stats, outfiletag=f'PPure_NativeBy{NativeBy}_{args.outname}', scope='last10')
        
        # combined probability stats for PLoss, PGain, and PPure
        combined_df = pd.concat([PLoss_stats, PGain_stats, PPure_stats])
        combined_df['Label'] = combined_df['Metric']
        print(f'CombinedProb:\n{combined_df}')
        anal.Plot_OneSampleStats(combined_df, outfiletag=f'CombinedProb_NativeBy{NativeBy}_{args.outname}', scope='last10')
      
        # get two sample stats for PLoss versus PGain
        anal.TwoSampleStats(PLoss_df, PGain_df, NativeBy=NativeBy, setID=setID, test='permutation', labels=['PLoss', 'PGain'], scope='last10')
        
        # probability of Loss considering only loss, gain changes in ent
        PLossNoPure_df = anal.CalcChangeMetrics(NativeBy=NativeBy, Metric='PLossNoPure', scope='last10')
        print(f'PLossNoPure_df:\n{PLossNoPure_df}')
        PLossNoPure_stats = anal.OneSampleStats(PLossNoPure_df, setID=setID, NativeBy=NativeBy, test_mean=0.5, label='PLossNoPure', scope='last10')
        anal.Plot_OneSampleStats(PLossNoPure_stats, outfiletag=f'PLossNoPure_NativeBy{NativeBy}_{args.outname}', scope='last10')

        # probability of Gain considering only loss, gain changes in ent
        PGainNoPure_df = anal.CalcChangeMetrics(NativeBy=NativeBy, Metric='PGainNoPure', scope='last10')
        print(f'PGainNoPure_df:\n{PGainNoPure_df}')
        PGainNoPure_stats = anal.OneSampleStats(PGainNoPure_df, setID=setID, NativeBy=NativeBy, test_mean=0.5, label='PGainNoPure', scope='last10')
        anal.Plot_OneSampleStats(PGainNoPure_stats, outfiletag=f'PGainNoPure_NativeBy{NativeBy}_{args.outname}', scope='last10')

        # combined probability stats for PLoss, PGain, and PPure
        combined_df = pd.concat([PLossNoPure_stats, PGainNoPure_stats])
        combined_df['Label'] = combined_df['Metric']
        print(f'CombinedProbNoPure:\n{combined_df}')
        anal.Plot_OneSampleStats(combined_df, outfiletag=f'CombinedProbNoPure_NativeBy{NativeBy}_{args.outname}', scope='last10')

        # get two sample stats for PLoss versus PGain
        anal.TwoSampleStats(PLossNoPure_df, PGainNoPure_df, NativeBy=NativeBy, setID=setID, test='permutation', labels=['PLossNoPure', 'PGainNoPure'], scope='last10')
        ###############################################################################


        ###############################################################################
        # Odds of Loss considering loss, gain, and pure switch changes in ent
        OddsLoss_df = anal.CalcChangeMetrics(NativeBy=NativeBy, Metric='OddsLoss', scope='last10')
        print(f'OddsLoss_df:\n{OddsLoss_df}')
        OddsLoss_stats = anal.OneSampleStats(OddsLoss_df, setID=setID, NativeBy=NativeBy, test_mean=0.5, label='OddsLoss', scope='last10')
        anal.Plot_OneSampleStats(OddsLoss_stats, outfiletag=f'OddsLoss_NativeBy{NativeBy}_{args.outname}', scope='last10')

        # Odds of Gain considering loss, gain, and pure switch changes in ent
        OddsGain_df = anal.CalcChangeMetrics(NativeBy=NativeBy, Metric='OddsGain', scope='last10')
        print(f'OddsGain_df:\n{OddsGain_df}')
        OddsGain_stats = anal.OneSampleStats(OddsGain_df, setID=setID, NativeBy=NativeBy, test_mean=0.5, label='OddsGain', scope='last10')
        anal.Plot_OneSampleStats(OddsGain_stats, outfiletag=f'OddsGain_NativeBy{NativeBy}_{args.outname}', scope='last10')

        # Odds of Pure considering loss, gain, and pure switch changes in ent
        OddsPure_df = anal.CalcChangeMetrics(NativeBy=NativeBy, Metric='OddsPure', scope='last10')
        print(f'OddsPure_df:\n{OddsPure_df}')
        OddsPure_stats = anal.OneSampleStats(OddsPure_df, setID=setID, NativeBy=NativeBy, test_mean=0.5, label='OddsPure', scope='last10')
        anal.Plot_OneSampleStats(OddsPure_stats, outfiletag=f'OddsPure_NativeBy{NativeBy}_{args.outname}', scope='last10')

        # combined probability stats for PLoss, PGain, and PPure
        combined_df = pd.concat([OddsLoss_stats, OddsGain_stats, OddsPure_stats])
        combined_df['Label'] = combined_df['Metric']
        print(f'CombinedOdds:\n{combined_df}')
        anal.Plot_OneSampleStats(combined_df, outfiletag=f'CombinedOdds_NativeBy{NativeBy}_{args.outname}', scope='last10')

        # get two sample stats for PLoss versus PGain
        anal.TwoSampleStats(OddsLoss_df, OddsGain_df, NativeBy=NativeBy, setID=setID, test='permutation', labels=['OddsLoss', 'OddsGain'], scope='last10')
      
        # Odds of Loss considering only loss, gain changes in ent
        OddsLossNoPure_df = anal.CalcChangeMetrics(NativeBy=NativeBy, Metric='OddsLossNoPure', scope='last10')
        print(f'OddsLossNoPure_df:\n{OddsLossNoPure_df}')
        OddsLossNoPure_stats = anal.OneSampleStats(OddsLossNoPure_df, setID=setID, NativeBy=NativeBy, test_mean=0.5, label='OddsLossNoPure', scope='last10')
        anal.Plot_OneSampleStats(OddsLossNoPure_stats, outfiletag=f'OddsLossNoPure_NativeBy{NativeBy}_{args.outname}', scope='last10')

        # Odds of Gain considering only loss, gain changes in ent
        OddsGainNoPure_df = anal.CalcChangeMetrics(NativeBy=NativeBy, Metric='OddsGainNoPure', scope='last10')
        print(f'OddsGainNoPure_df:\n{OddsGainNoPure_df}')
        OddsGainNoPure_stats = anal.OneSampleStats(OddsGainNoPure_df, setID=setID, NativeBy=NativeBy, test_mean=0.5, label='OddsGainNoPure', scope='last10')
        anal.Plot_OneSampleStats(OddsGainNoPure_stats, outfiletag=f'OddsGainNoPure_NativeBy{NativeBy}_{args.outname}', scope='last10')

        # combined probability stats for PLoss, PGain, and PPure
        combined_df = pd.concat([OddsLossNoPure_stats, OddsGainNoPure_stats])
        combined_df['Label'] = combined_df['Metric']
        print(f'CombinedOddsNoPure:\n{combined_df}')
        anal.Plot_OneSampleStats(combined_df, outfiletag=f'CombinedOddsNoPure_NativeBy{NativeBy}_{args.outname}', scope='last10')

        # get two sample stats for PLoss versus PGain
        anal.TwoSampleStats(OddsLossNoPure_df, OddsGainNoPure_df, NativeBy=NativeBy, setID=setID, test='permutation', labels=['OddsLossNoPure', 'OddsGainNoPure'], scope='last10')
        ###############################################################################

       
    print(f'logfile: {logfile}')

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()

print(f'NORMAL TERMINATION: {time.time() - start_time}')
logging.info(f'NORMAL TERMINATION: {time.time() - start_time}')