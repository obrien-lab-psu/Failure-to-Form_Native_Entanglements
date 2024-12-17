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
from scipy.stats import mode, permutation_test, ttest_1samp, ttest_ind
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

        ### Load candidate set 1 collected OP data
        set1_OPfile = os.path.join(self.CollectedOPpath, f'setID1/DATA/Quench_Collected_GQK.csv')
        print(f'set1_OPfile: {set1_OPfile}')
        if not os.path.exists(set1_OPfile):
            raise ValueError(f'No set 1 OP file!')
        else:
            self.set1_OP = pd.read_csv(set1_OPfile)
            self.set1_OP = self.set1_OP[self.set1_OP['Mirror'] == False]
            print(f'set1_OP:\n{self.set1_OP}')

        
        ### Load candidate set 1 collected OP data
        set2_OPfile = os.path.join(self.CollectedOPpath, f'setID2/DATA/Quench_Collected_GQK.csv')
        print(f'set2_OPfile: {set2_OPfile}')
        if not os.path.exists(set2_OPfile):
            raise ValueError(f'No set 1 OP file!')
        else:
            self.set2_OP = pd.read_csv(set2_OPfile)
            self.set2_OP = self.set2_OP[self.set2_OP['Mirror'] == False]
            print(f'set2_OP:\n{self.set2_OP}')

    #######################################################################################  

    #######################################################################################
    def CalcFractMisfolded(self, NativeBy='MSM', scope='full'):
        """
        Calculate the ratio of misfolded and native frames in each trajectory 
        Do it for the full trajectory and every 10% percentile buckets
        Apply a blanket Q threshold and only consider structures with Q > 0.6
        NativeBy: can be MSM, Ref, or None
            MSM = using the metastable state with the highest <Q> and lowest <G> as the native state. Will examine all frames outside this native meta stable state.
            Ref = using <Qref> - 3*sigma and <Gref> + 3*sigma as thresholds for selecting native frames. Will examine all frames with Q greater than the threshold and G less than. 
            None = No native defined frames. examine all frames. 
        scope: defines the scope of the analysis. full analyzes the whole GQ data and last10 only uses the last 10% of the traj
        """

        num_chunks = 10
        #print(f'num_chunks: {num_chunks}')
        #################################################
        ## Determine the scope of the analysis
        if scope in ['full', 'last10']:
            if scope == 'full':
                scope_dict = {0:'full', 1:'0-10%', 2:'10-20%', 3:'20-30%', 4:'30-40%', 5:'40-50%', 6:'50-60%', 7:'60-70%', 8:'70-80%', 9:'80-90%', 10:'90-100%'}
                fract_misfolded_df = {'tag':[], 'setID':[], 'traj':[], 'NativeBy':[], 'Metric':[],
                    'full':[], '0-10%':[], '10-20%':[], '20-30%':[], '30-40%':[], 
                    '40-50%':[], '50-60%':[], '60-70%':[], '70-80%':[], '80-90%':[], '90-100%':[]}
            else:
                scope_dict = {0:'full', 10:'90-100%'}
                fract_misfolded_df = {'tag':[], 'setID':[], 'traj':[], 'NativeBy':[], 'Metric':[],
                    'full':[], '90-100%':[]}
        else:
            raise ValueError(f'Scope can only be full or last10')


        #####################################################
        ## Get the frame chunks to split the analysis on
        FrameEnd = 26666
        print(f'FrameEnd: {FrameEnd}')
        # Split the array into 10 equal chunks
        frame_chunks = [np.arange(0, FrameEnd + 1)] + np.array_split(np.arange(0, FrameEnd + 1), 10)
        print(f'frame_chunks: {frame_chunks}')

        ####################################################
        # Calculate R for each traj in sets 1 and 2
        for gene, pdb, chain, setID, set_name in self.candidates.values:

            tag = f'{gene}_{pdb}_{chain}'
            print(gene, pdb, chain, setID, set_name)

            if setID == 1:
                loc_df = self.set1_OP[self.set1_OP['gene'] == gene]
            elif setID == 2:
                loc_df = self.set2_OP[self.set2_OP['gene'] == gene]
            else:
                continue
            #print(f'loc_df:\n{loc_df}')

            for traj, traj_df in loc_df.groupby('traj'):
                #print(traj_df)
                fract_misfolded_df['tag'] += [tag]
                fract_misfolded_df['setID'] += [setID]
                fract_misfolded_df['traj'] += [traj]
                fract_misfolded_df['NativeBy'] += [NativeBy]
                fract_misfolded_df['Metric'] += ['FracMisfolded']

                for i,label in scope_dict.items():
                    # Calculate the percentile range
                    lower_bound = frame_chunks[i][0]
                    upper_bound = frame_chunks[i][-1]
                    lower_bound, upper_bound = int(lower_bound), int(upper_bound)
                    #print(i, lower_bound, upper_bound)
                    
                    # Filter the DataFrame for rows within this percentile range
                    chunk = traj_df[(traj_df['frame'] >= lower_bound) & (traj_df['frame'] < upper_bound)]
                    #print(chunk)
                    if len(chunk) == 0:
                        R = np.nan
                    else:
                        R = 1 - np.mean(chunk[f'NativeBy{NativeBy}'].values)
                    fract_misfolded_df[label] += [R]

        fract_misfolded_df = pd.DataFrame(fract_misfolded_df)           
        logging.info(f'fract_misfolded_df:\n{fract_misfolded_df}')
        outfile = os.path.join(self.data_path, f'FracMisfolded_Scope-{scope}_NativeBy{NativeBy}.csv')
        fract_misfolded_df.to_csv(outfile, index=False)
        logging.info(f'SAVED: {outfile}')
        return fract_misfolded_df
    ########################################################################################

    #######################################################################################
    def OneSampleStats(self, df, setID=2, NativeBy='MSM', test_mean=1, label='None', scope='full'):
        """
        Calculate the OneSampleStats for any metric: Avg, Median, std, 95% ci of the mean, 1sample ttest stat and pvalue
        setID is an integer to be used to select the data
        df is a metric dataframe 
        NativeBy: can be MSM, Ref, or None
            MSM = using the metastable state with the highest <Q> and lowest <G> as the native state. Will examine all frames outside this native meta stable state.
            Ref = using <Qref> - 3*sigma and <Gref> + 3*sigma as thresholds for selecting native frames. Will examine all frames with Q greater than the threshold and G less than. 
            None = No native defined frames. examine all frames. 
        test_mean is the population mean used in the 1sample ttest
        outfiletag is just a string appended to the image file name
        scope is the level at which data was analyzed. all chunks of the trajectory were analyzed (full) or just the last 10% (last10)
        """
        logging.info(f'Performing 1 sample T-Test on {label} data for setID: {setID} and NativeBy: {NativeBy}')
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
        return df_stats   
    #######################################################################################

    #######################################################################################
    def TwoSampleStats(self, df, setIDs=(1,2), NativeBy='MSM', test='permutation', label='Fract_misfolded', scope='full'):
        """
        Calculate the TwoSampleStats for any metric between two sets of data using either permutation or ttest
        setIDs is a tuple of two setIDs to be used to select the data
        df is a metric dataframe 
        NativeBy: can be MSM, Ref, or None
            MSM = using the metastable state with the highest <Q> and lowest <G> as the native state. Will examine all frames outside this native meta stable state.
            Ref = using <Qref> - 3*sigma and <Gref> + 3*sigma as thresholds for selecting native frames. Will examine all frames with Q greater than the threshold and G less than. 
            None = No native defined frames. examine all frames. 
        test is permutation or ttest depending on which test you want to perform
        outfiletag is just a string appended to the image file name
        scope is the level at which data was analyzed. all chunks of the trajectory were analyzed (full) or just the last 10% (last10)
        """
        print(f'Performing {test} on {label} data for setID: {setIDs} and NativeBy: {NativeBy}')
        df = df[(df['setID'].isin(setIDs)) & (df['NativeBy']==NativeBy)]
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
        df_stats = {'setIDs':[], 'NativeBy':[], 'Label':[], 'test':[], 'statistic':[], 'pvalue':[]}

        for key in scope_list:

            data1 = df[df['setID'] == setIDs[0]][key].values
            data2 = df[df['setID'] == setIDs[1]][key].values
            
            if test == 'permutation':
                res = permutation_test((data1, data2), statistic, vectorized=True, alternative='greater')
            elif test == 'ttest':
                res = ttest_ind(data1, data2, alternative='greater', equal_var=False)
            statistic_val = res.statistic
            pvalue = res.pvalue

            df_stats['setIDs'] += [setIDs]
            df_stats['NativeBy'] += [NativeBy]
            df_stats['Label'] += [key]
            df_stats['test'] += [test]
            df_stats['statistic'] += [statistic_val]
            df_stats['pvalue'] += [pvalue]

        df_stats = pd.DataFrame(df_stats)
        print(df_stats)
        outfile = os.path.join(self.data_path, f'{label}_{test}_TwoSampleStats_Scope-{scope}_NativeBy{NativeBy}.csv')
        df_stats.to_csv(outfile, index=False)
        logging.info(f'SAVED: {outfile}')  
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
        print(f'Plotting:\n{df}')
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
        logging.info(f'SAVED: {outfile}')
        print(f'SAVED: {outfile}')

        df.to_csv(outfile_csv, index=False)
        print(f'SAVED: {outfile_csv}')
    #######################################################################################

    ########################################################################################
    def Plot_TwoSampleStats(self, df1, df2, test_df, outfiletag='temp', scope='full'):     
        """
        Plots the two sample stats for two sets of data.
        df1 is a OneSampleStats dataframe with the mean, median, std, and conf. intervals for dataset 1
        df2 is a OneSampleStats dataframe with the mean, median, std, and conf. intervals for dataset 2
        test_df is the TwoSampleStats dataframe with the test statistic and pvalue for the comparison of datasets 1 and 2
        outfiletag is just a string appended to the image file name
        scope is the level at which data was analyzed. all chunks of the trajectory were analyzed (full) or just the last 10% (last10)
        """       
        ########################################################################################
        ## calculate the statstics between the two setIDs 1 and 2 for Rfull and all the blocks
        # save both a image with the histograms and a csv file with the results of the permutation test
        #print(f'df1:\n{df1}')
        #print(f'df2:\n{df2}')
        #print(f'test_df:\n{test_df}')

        #print(f'Plotting:\n{df}')
        outfile = os.path.join(self.plot_path, f'{outfiletag}_CompStats_Scope-{scope}_Plot.png')
        outfile_csv = os.path.join(self.plot_path, f'{outfiletag}_CompStats_Scope-{scope}_Plot.csv')

        #################################################
        ## Determine the scope of the analysis
        if scope in ['full', 'last10']:
            if scope == 'full':
                fig, axes = plt.subplots(4, 1, figsize=(4, 6), sharex=True)
            else:
                fig, axes = plt.subplots(4, 1, figsize=(2, 6), sharex=True)
        else:
            raise ValueError(f'Scope can only be full or last10')
        
        # First plot: avg with error bars (lb and ub)
        axes[0].errorbar(df1['Label'], df1['avg'], yerr=[df1['avg'] - df1['lb'], df1['ub'] - df1['avg']], 
                        fmt='o', capsize=5, label=f'Low Mis. Prone', color='red')
        axes[0].errorbar(df2['Label'], df2['avg'], yerr=[df2['avg'] - df2['lb'], df2['ub'] - df2['avg']], 
                        fmt='o', capsize=5, label='High Mis. Prone', color='blue')
        axes[0].set_ylabel('Average (avg)')
        #axes[0].set_title('Average with Error Bars')
        
        # Second plot: median
        axes[1].plot(df1['Label'], df1['median'], 'o', color='red', label=f'Low Mis. Prone')
        axes[1].plot(df2['Label'], df2['median'], 'o', color='blue', label=f'High Mis. Prone')
        axes[1].set_ylabel('Median')
        axes[0].legend(fontsize=8)

        # Third plot: pvalue with log scale
        axes[2].plot(test_df['Label'], test_df['statistic'], 'o', color='black', label='Eff. Size')
        axes[2].set_ylabel('Eff. Size')

        # Fourth plot: pvalue with log scale
        axes[3].plot(test_df['Label'], test_df['pvalue'], 'o', color='black', label='P-Value')
        axes[3].set_yscale('log')
        axes[3].set_ylabel('P-Value (log scale)')
        axes[3].axhline(0.05)

        # Set shared x-axis labels
        plt.xticks(rotation=45)
        #axes[2].set_xlabel('Label')
        # Adjust alignment of x-tick labels to be anchored at the upper right
        for tick in axes[3].get_xticklabels():
            tick.set_ha('right')  # Set horizontal alignment to right
            tick.set_va('top')    # Set vertical alignment to top
        
        # Layout adjustments
        plt.suptitle(outfiletag)
        plt.tight_layout()
        plt.savefig(outfile)
        logging.info(f'SAVED: {outfile}')
        print(f'SAVED: {outfile}')
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

############## MAIN #################
def main():
    
    script_name = f'CompareMisfoldingPropensity'
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("--candidates", type=str, required=True, help="A file containing two columns. The candidate tag and the groupID")
    parser.add_argument("--CollectedOPpath", type=str, required=True, help="path to the CollectAndProcessesOP DATA/ folder")
    parser.add_argument("--outname", type=str, required=True, help="base name for output files")

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
        # Fraction of misfolded frames and compare the distributions from set 1 and 2
        Fract_misfolded_df = anal.CalcFractMisfolded(NativeBy=NativeBy, scope='full') 
        print(f'Fract_misfolded_df:\n{Fract_misfolded_df}')
        Fract_misfolded_set1_stats = anal.OneSampleStats(Fract_misfolded_df, setID=1, NativeBy=NativeBy, test_mean=None, label='Fract_misfolded', scope='full')
        Fract_misfolded_set2_stats = anal.OneSampleStats(Fract_misfolded_df, setID=2, NativeBy=NativeBy, test_mean=None, label='Fract_misfolded', scope='full')
        anal.Plot_OneSampleStats(Fract_misfolded_set1_stats, outfiletag=f'Fract_misfolded_set1_NativeBy{NativeBy}', scope='full')
        anal.Plot_OneSampleStats(Fract_misfolded_set2_stats, outfiletag=f'Fract_misfolded_set2_NativeBy{NativeBy}', scope='full')

        perm_df = anal.TwoSampleStats(Fract_misfolded_df, setIDs=(1,2), NativeBy=NativeBy, test='permutation', label='Fract_misfolded', scope='full')
        ttest_df = anal.TwoSampleStats(Fract_misfolded_df, setIDs=(1,2), NativeBy=NativeBy, test='ttest', label='Fract_misfolded', scope='full')
        anal.Plot_TwoSampleStats(Fract_misfolded_set1_stats, Fract_misfolded_set2_stats, perm_df, outfiletag=f'Fract_misfolded_set1V2_NativeBy{NativeBy}', scope='full')
        ###############################################################################


    print(f'logfile: {logfile}')

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()

print(f'NORMAL TERMINATION: {time.time() - start_time}')
logging.info(f'NORMAL TERMINATION: {time.time() - start_time}')