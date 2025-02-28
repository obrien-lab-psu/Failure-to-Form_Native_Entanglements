import pandas as pd
import glob
from scipy.stats import permutation_test, mode
import numpy as np
import time
import argparse
import sys, os, logging
from multiprocessing import Pool, cpu_count
import matplotlib.pyplot as plt
from scipy.ndimage import label
from types import SimpleNamespace

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
        ("--toplevel", type=str, required=True, help="file containing relative paths to either native state GQ files or the MSM file for various proteins")
        ("--outname", type=str, required=True, help="base name for output files")
        ("--setID", type=int, required=True, help="setID to analyze (2 or 3)")
        """

        # parse the parameters 
        self.setID = args.setID
        self.candidates = pd.read_csv(args.candidates)
        self.candidates = self.candidates[self.candidates['set'] == self.setID]
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

        self.toplevel = args.toplevel
        logging.info(f'toplevel: {self.toplevel}')
        print(f'toplevel: {self.toplevel}')

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
        1. Load the mirror file containing those hand annotated mirror trajectories
        2. Load the collected GQK dataframe
        3. 
        """
        logging.info('Load data')

        ###############################################################################################
        ## Load the Quench Trajectory stats
        set3_files = glob.glob(os.path.join(self.CollectedOPpath, f'setID3/DATA/*_Quench_Collected_Traj_GQK_STATS.csv'))
        #print(f'set3_files: {set3_files}')
        self.Traj_GQK_STATS = []
        for f in set3_files:
            df = pd.read_csv(f)
            self.Traj_GQK_STATS += [df]
        self.Traj_GQK_STATS = pd.concat(self.Traj_GQK_STATS)
        self.Traj_GQK_STATS = self.Traj_GQK_STATS[self.Traj_GQK_STATS['NativeFolded'] == False]
        print(f'Traj_GQK_STATS:\n{self.Traj_GQK_STATS}')

        traj_stats_outfile = os.path.join(self.data_path, f'Combined_Quench_Collected_Traj_GQK_STATS.csv')
        self.Traj_GQK_STATS.to_csv(traj_stats_outfile, index=False)
        print(f'SAVED: {traj_stats_outfile}')
        ###############################################################################################

        ###############################################################################################
        ## Load the Quench Trajectory Lifetime info and get only the misfolded traj
        set3_files = glob.glob(os.path.join(self.CollectedOPpath, f'setID3/DATA/*_ChangeType_Lifetimes_summary_setID3.csv'))
        #print(f'set3_files: {set3_files}')
        self.ChangeType_Lifetimes = []
        for f in set3_files:
            df = pd.read_csv(f)
            self.ChangeType_Lifetimes += [df]
        self.ChangeType_Lifetimes = pd.concat(self.ChangeType_Lifetimes)
        print(f'ChangeType_Lifetimes:\n{self.ChangeType_Lifetimes}')

        misfolded_dfs = []
        for gene, traj in self.Traj_GQK_STATS[['gene', 'traj']].values:
            tmp = self.ChangeType_Lifetimes[(self.ChangeType_Lifetimes['gene'] == gene) & (self.ChangeType_Lifetimes['traj'] == traj)]
            if len(tmp) == 0:
                raise ValueError('empty lifetime traj but GQK stats say it is misfolded!')
            
            misfolded_dfs += [tmp]
        self.ChangeType_Lifetimes = pd.concat(misfolded_dfs)
        print(f'ChangeType_Lifetimes:\n{self.ChangeType_Lifetimes}')

        ChangeType_Lifetimes_summary_outfile = os.path.join(self.data_path, f'Combined_ChangeType_Lifetimes_summary_setID3.csv')
        self.ChangeType_Lifetimes.to_csv(ChangeType_Lifetimes_summary_outfile, index=False)
        print(f'SAVED: {ChangeType_Lifetimes_summary_outfile}') 
        ###############################################################################################

        ###############################################################################################
        ## Load the structure summary for onlyLoss, onlyGain, and BothLoss,Gain info
        set3_files = glob.glob(os.path.join(self.CollectedOPpath, f'setID3/DATA/*_StructureCategory_summary_setID3.csv'))
        #print(f'set3_files: {set3_files}')
        self.StructureCategory = []
        for f in set3_files:
            df = pd.read_csv(f)
            self.StructureCategory += [df]
        self.StructureCategory = pd.concat(self.StructureCategory)
        print(f'StructureCategory:\n{self.StructureCategory}')

        misfolded_dfs = []
        for gene, traj in self.Traj_GQK_STATS[['gene', 'traj']].values:
            tmp = self.StructureCategory[(self.StructureCategory['gene'] == gene) & (self.StructureCategory['traj'] == traj)]
            if len(tmp) == 0:
                print(gene, traj)
                continue
            
            misfolded_dfs += [tmp]
        self.StructureCategory = pd.concat(misfolded_dfs)
        print(f'StructureCategory:\n{self.StructureCategory}')   

        StructureCategory_summary_outfile = os.path.join(self.data_path, f'Combined_StructureCategory_summary_setID3.csv')
        self.StructureCategory.to_csv(StructureCategory_summary_outfile, index=False)
        print(f'SAVED: {StructureCategory_summary_outfile}') 
        quit()
        ###############################################################################################

    #######################################################################################

    #######################################################################################  
    def plot_lifetimes_cdf(self,):
        """
        Plot the CDF for the change in entanglement life times
        """
        Lifetimes = self.ChangeType_Lifetimes[self.ChangeType_Lifetimes['ChangeType'] == 'TotalChanges']
        print(Lifetimes)

        output_file_png = os.path.join(self.plot_path, f'TotalChanges_Lifetimes_CDF.png')
        output_file_svg = os.path.join(self.plot_path, f'TotalChanges_Lifetimes_CDF.svg')
        output_file_csv = os.path.join(self.plot_path, f'TotalChanges_Lifetimes_CDF.csv')

        sorted_data = np.sort(Lifetimes['meanLifetime'].values)
        # Calculate the CDF values
        cdf_y = np.arange(1, len(sorted_data) + 1) / len(sorted_data)
        cdf_x = sorted_data
        
        # Plot the CDF
        plt.figure(figsize=(8, 6))
        plt.plot(cdf_x, cdf_y, marker='.', linestyle='none', label='CDF')
        plt.xlabel('Lifetime, (ns)')
        plt.ylabel('Cumulative Probability (CDF)')
        #plt.title('Cumulative Distribution Function (CDF)')
        plt.grid(True)
        #plt.legend()
        
        # Save the plot
        plt.savefig(output_file_png)
        plt.savefig(output_file_svg)
        plt.close()  # Close the figure to free memory
        print(f'SAVED: {output_file_png}')

        df = pd.DataFrame({'x':cdf_x, 'y':cdf_y})
        self.lifetimes_cdf = df
        df.to_csv(output_file_csv, index=False)
        print(f'SAVED: {output_file_csv}')
    ##########################################################################################

    #######################################################################################
    def calc_thresholding_metrics(self,):
        """
        get the percentiles of the lifetimes and scan 
        for each percentile calculate the following
        <Qmode>
        <Zeta>
        <FractionTotalLoss>
        <FractionTotalGain>
        and all its components
        """
        change_keys = ['TotalChanges', 'TotalLoss', 'TotalGain', 'TrueLoss', 'PartialLoss', 'TrueGain', 'PartialGain', 'TrueLossOverlap', 'PartialLossOverlap', 'TrueGainOverlap', 'PartialGainOverlap']
        overlap_keys = ['TotalChanges', 'TotalLoss', 'TotalGain', 'TrueLoss', 'PartialLoss', 'TrueGain', 'PartialGain', 'TrueLossOverlap', 'PartialLossOverlap', 'TrueGainOverlap', 'PartialGainOverlap']

        #print(f'Calculating thresholds and metrics df:\n{self.ChangeType_Lifetimes}')
        Mechanism_stats_summary_outfile = os.path.join(self.data_path, f'Mechanism_stats_summary_setID{self.setID}.csv')
        print(f'Mechanism_stats_summary_outfile: {Mechanism_stats_summary_outfile}')
        if not os.path.exists(Mechanism_stats_summary_outfile):
            column = 'meanLifetime'
            stats_df = {'lower_bound':[], 'upper_bound':[], 'type':[], 'mean':[], 'lb':[], 'ub':[], 'n':[]}
            Lifetimes = self.ChangeType_Lifetimes[self.ChangeType_Lifetimes['ChangeType'] == 'TotalChanges']
            Lifetimes[column] = Lifetimes[column].fillna(0.0)
            Lifetimes = Lifetimes[column].values

            bounds = [(0,100), (100, 200.025)]
            print(f'bounds: {bounds}')
            logging.info(f'bounds: {bounds}')

            n_resamples = 10000

            # Loop through each 10th percentile
            pvalues = []
            for i, (lower_bound, upper_bound) in enumerate(bounds):

                ################################################################################
                # Filter rows within this percentile range
                Mech_chunk = []
                QZ_chunk = []
                Struct_chunk = []
                n = 0
                for rowi, row in self.ChangeType_Lifetimes.iterrows():
                    if row['ChangeType'] == 'TotalChanges':
                        if row['meanLifetime'] > lower_bound and row['meanLifetime'] <= upper_bound:
                            gene = row['gene']
                            traj = row['traj']
                            n += 1

                            tmp = self.ChangeType_Lifetimes[(self.ChangeType_Lifetimes['gene'] == gene) & (self.ChangeType_Lifetimes['traj'] == traj)]
                            Mech_chunk += [tmp]

                            tmp = self.Traj_GQK_STATS[(self.Traj_GQK_STATS['gene'] == gene) & (self.Traj_GQK_STATS['traj'] == traj)]
                            QZ_chunk += [tmp]

                            tmp = self.StructureCategory[(self.StructureCategory['gene'] == gene) & (self.StructureCategory['traj'] == traj)]
                            Struct_chunk += [tmp]                            
             
                Mech_chunk = pd.concat(Mech_chunk)
                QZ_chunk = pd.concat(QZ_chunk)
                Struct_chunk = pd.concat(Struct_chunk)

                # Process the chunk (e.g., print or perform operations)
                print(f"\nChunk {i + 1}: {lower_bound:.6f} to {upper_bound:.6f}")
                logging.info(f"\nChunk {i + 1}: {lower_bound:.6f} to {upper_bound:.6f}")
                print(Mech_chunk)
                print(QZ_chunk)
                print(Struct_chunk)

                ################################################################################
                ## get stats for each change type fraction
                change_keys = ['TotalChanges', 'TotalLoss', 'TotalGain', 'TrueLoss', 'PartialLoss', 'TrueGain', 'PartialGain']
                for key in change_keys:
                    key_chunk = Mech_chunk[Mech_chunk['ChangeType'] == key]
                    key_chunk.loc[:, 'meanFract'] = key_chunk['meanFract'].fillna(0.0)
                    #print(key_chunk)
                    key_chunk_vals = key_chunk['meanFract'].values
                    median, mean, std, lb, ub = get_stats(key_chunk_vals, n_resamples=n_resamples)
                    logging.info(f'Fraction {key} stats: {mean} ({lb}, {ub})')
                    print(f'Fraction {key} stats: {mean} ({lb}, {ub})')
                    #stats_df = {'lower_bound':[], 'upper_bound':[], 'type':[], 'mean':[], 'lb':[], 'ub':[]}
                    stats_df['lower_bound'] += [lower_bound]
                    stats_df['upper_bound'] += [upper_bound]
                    stats_df['type'] += [key]
                    stats_df['mean'] += [mean]
                    stats_df['lb'] += [lb]
                    stats_df['ub'] += [ub]
                    stats_df['n'] += [n] 

                ## get stats for Loss overlap fraction and 1 minus that
                overlap_keys = ['TrueLossOverlap']
                for key in overlap_keys:
                    key_chunk = Mech_chunk[Mech_chunk['ChangeType'] == key]
                    key_chunk.loc[:, 'meanFract'] = key_chunk['meanFract'].fillna(0.0)
                    #print(key_chunk)
                    key_chunk_vals = key_chunk['meanFract'].values
                    median, mean, std, lb, ub = get_stats(key_chunk_vals, n_resamples=n_resamples)
                    logging.info(f'Fraction {key} stats: {mean} ({lb}, {ub})')
                    print(f'Fraction {key} stats: {mean} ({lb}, {ub})')
                    #stats_df = {'lower_bound':[], 'upper_bound':[], 'type':[], 'mean':[], 'lb':[], 'ub':[]}
                    stats_df['lower_bound'] += [lower_bound]
                    stats_df['upper_bound'] += [upper_bound]
                    stats_df['type'] += [key]
                    stats_df['mean'] += [mean]
                    stats_df['lb'] += [lb]
                    stats_df['ub'] += [ub]
                    stats_df['n'] += [n] 

                    key_chunk_vals_minus = 1 - key_chunk_vals
                    median, mean, std, lb, ub = get_stats(key_chunk_vals_minus, n_resamples=n_resamples)
                    logging.info(f'Fraction 1 - {key} stats: {mean} ({lb}, {ub})')
                    print(f'Fraction 1 - {key} stats: {mean} ({lb}, {ub})')
                    #stats_df = {'lower_bound':[], 'upper_bound':[], 'type':[], 'mean':[], 'lb':[], 'ub':[]}
                    stats_df['lower_bound'] += [lower_bound]
                    stats_df['upper_bound'] += [upper_bound]
                    stats_df['type'] += [f'{key}_minusOne']
                    stats_df['mean'] += [mean]
                    stats_df['lb'] += [lb]
                    stats_df['ub'] += [ub]
                    stats_df['n'] += [n]                 

                for key in ['meanQmode', 'meanZmode']:
                    key_chunk_vals = QZ_chunk[key].values
                    median, mean, std, lb, ub = get_stats(key_chunk_vals, n_resamples=n_resamples)
                    logging.info(f'Fraction {key} stats: {mean} ({lb}, {ub})')
                    print(f'Fraction {key} stats: {mean} ({lb}, {ub})')
                    stats_df['lower_bound'] += [lower_bound]
                    stats_df['upper_bound'] += [upper_bound]
                    stats_df['type'] += [key]
                    stats_df['mean'] += [mean]
                    stats_df['lb'] += [lb]
                    stats_df['ub'] += [ub]   
                    stats_df['n'] += [n]     

                for key in ['LossOnly', 'GainOnly', 'BothLossGain']:
                    key_chunk = Struct_chunk[Struct_chunk['ChangeType'] == key]
                    #key_chunk.loc[:, 'mean'] = key_chunk['mean'].fillna(0.0)
                    #print(key_chunk)
                    key_chunk_vals = key_chunk['mean'].values
                    median, mean, std, lb, ub = get_stats(key_chunk_vals, n_resamples=n_resamples)
                    logging.info(f'Fraction {key} stats: {mean} ({lb}, {ub})')
                    print(f'Fraction {key} stats: {mean} ({lb}, {ub})')
                    #stats_df = {'lower_bound':[], 'upper_bound':[], 'type':[], 'mean':[], 'lb':[], 'ub':[]}
                    stats_df['lower_bound'] += [lower_bound]
                    stats_df['upper_bound'] += [upper_bound]
                    stats_df['type'] += [key]
                    stats_df['mean'] += [mean]
                    stats_df['lb'] += [lb]
                    stats_df['ub'] += [ub]
                    stats_df['n'] += [n] 
        
                ## Permutations
                pvalues += [np.nan]*14
                for (x, y) in [('LossOnly', 'GainOnly')]:
                    x_vals = Struct_chunk[Struct_chunk['ChangeType'] == x]['mean'].values
                    y_vals = Struct_chunk[Struct_chunk['ChangeType'] == y]['mean'].values
                    res = permutation_testing((x_vals,y_vals))
                    print(res)
                    pvalues[-3] = res.pvalue
                    pvalues[-2] = res.pvalue

            stats_df['pvalue'] = pvalues
            stats_df = pd.DataFrame(stats_df)
            logging.info(f'stats_df:\n{stats_df.to_string()}')
            print(f'stats_df:\n{stats_df.to_string()}')
            ## save the Mechanism stats summary file
            Mechanism_stats_summary_outfile = os.path.join(self.data_path, f'Mechanism_stats_summary_setID{self.setID}.csv')
            stats_df.to_csv(Mechanism_stats_summary_outfile, index=False)
            logging.info(f'SAVED: {Mechanism_stats_summary_outfile}')
            print(f'SAVED: {Mechanism_stats_summary_outfile}')

        else:
            Mechanism_stats_summary_outfile = os.path.join(self.data_path, f'Mechanism_stats_summary_setID{self.setID}.csv')
            stats_df = pd.read_csv(Mechanism_stats_summary_outfile)
            logging.info(f'LOADED: {Mechanism_stats_summary_outfile}')
            print(f'LOADED: {Mechanism_stats_summary_outfile}')
        return stats_df
        #######################################################################################
    #######################################################################################

    #######################################################################################
    def plot_thresholding_metrics(self, df):
        """
        Create a two-row plot of 'mean' vs. 'type' with error bars using [lb, ub].
        The first row uses rows where lower_bound == 0, and the second where lower_bound == 100.

        Parameters:
            df (pd.DataFrame): The input DataFrame.
        
        Plot 1 row by 
        """
        # Filter rows based on lower_boun
        # Setup the figure and axes
        print(df)
        mapper = {0:'Short', 100:'Long'}
        df['xlables'] = [mapper[i] for i in df['lower_bound'].values]
        print(df)
        plt.rcParams.update({'font.size': 7})
        fig, axes = plt.subplots(1, 3, figsize=(7, 3))

        # Plot the CDF
        axes[0].plot(self.lifetimes_cdf['x'], self.lifetimes_cdf['y'], marker='.', linestyle='none', label='CDF')
        axes[0].set_xlabel('Lifetime, (ns)')
        axes[0].set_ylabel('Cumulative Probability (CDF)')
        axes[0].tick_params(axis='y', labelsize=6)
        axes[0].tick_params(axis='x', labelsize=6)
        axes[0].spines['right'].set_visible(False)
        axes[0].spines['top'].set_visible(False)
        axs_position = axes[0].get_position()
        bbox_in_fig_coords = axes[0].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
        fig.text(bbox_in_fig_coords.x0, 1, 'a', fontsize=8, fontweight='bold', va='top', ha='left')


        meanQmode = df[df['type'].isin(['meanQmode'])]
        meanQmode_yerr=[meanQmode['mean'] - meanQmode['lb'], meanQmode['ub'] - meanQmode['mean']]
        axes[1].errorbar(x=meanQmode['xlables'], y=meanQmode['mean'], yerr=meanQmode_yerr,fmt='o', capsize=5)
        axes[1].set_ylabel(r"$\langle Q_{\mathrm{mode}} \rangle$")
        axes[1].set_xlim(-0.5, 1.5)
        axes[1].tick_params(axis='y', labelsize=6)
        axes[1].tick_params(axis='x', labelsize=6)
        axes[1].spines['right'].set_visible(False)
        axes[1].spines['top'].set_visible(False)
        axs_position = axes[1].get_position()
        bbox_in_fig_coords = axes[1].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
        fig.text(bbox_in_fig_coords.x0, 1, 'b', fontsize=8, fontweight='bold', va='top', ha='left')

        meanZmode = df[df['type'].isin(['meanZmode'])]
        meanZmode_yerr=[meanZmode['mean'] - meanZmode['lb'], meanZmode['ub'] - meanZmode['mean']]
        axes[2].errorbar(x=meanZmode['xlables'], y=meanZmode['mean'], yerr=meanZmode_yerr,fmt='o', capsize=5)
        axes[2].set_ylabel(r"$\langle \zeta_{\mathrm{mode}} \rangle$")
        axes[2].set_xlim(-0.5, 1.5)
        axes[2].tick_params(axis='y', labelsize=6)
        axes[2].tick_params(axis='x', labelsize=6)
        axes[2].spines['right'].set_visible(False)
        axes[2].spines['top'].set_visible(False)
        axs_position = axes[2].get_position()
        bbox_in_fig_coords = axes[2].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
        fig.text(bbox_in_fig_coords.x0, 1, 'c', fontsize=8, fontweight='bold', va='top', ha='left')

        plt.tight_layout()
        Mechanism_stats_summary_outfile_png = os.path.join(self.data_path, f'Mechanism_stats_summary_setID{self.setID}.svg')
        plt.savefig(Mechanism_stats_summary_outfile_png)
        print(f'SAVED: {Mechanism_stats_summary_outfile_png}')

        Mechanism_stats_summary_outfile_png = os.path.join(self.data_path, f'Mechanism_stats_summary_setID{self.setID}.png')
        plt.savefig(Mechanism_stats_summary_outfile_png)
        print(f'SAVED: {Mechanism_stats_summary_outfile_png}')
        quit()
    #######################################################################################

#######################################################################################
def get_stats(arr, n_resamples=100, metric='mean'):
    """
    Get the <> and 95% ci for a stats array
    """
    if metric == 'mean':
        mean = np.mean(arr)
        std = np.std(arr)
        (lb, ub) = bootstrap(arr, n_resamples=100)
        median = np.median(arr)
        return (median, mean, std, lb, ub)
    
    elif metric == 'sum':
        sum = np.sum(arr)
        std = np.std(arr)
        (lb, ub) = bootstrap(arr, n_resamples=100, metric='sum')
        median = np.median(arr)
        return (median, sum, std, lb, ub)       
#######################################################################################

#######################################################################################
def bootstrap(data, n_resamples=10000, metric='mean'):
    boot_vals = []
    for b in range(n_resamples):
        boot_samp = np.random.choice(data, size=len(data))

        if metric == 'mean':
            boot_metric = np.mean(boot_samp)
        elif metric == 'sum':
            boot_metric = np.sum(boot_samp)

        #print(b, boot_metric)
        boot_vals += [boot_metric]

    lb = np.percentile(boot_vals, 2.5)
    ub = np.percentile(boot_vals, 97.5)
    return (lb, ub)
#######################################################################################

#######################################################################################
def permutation_testing(data, n_resamples=10000):
    perm_vals = []
    combined = np.hstack(data)
    d1N = len(data[0])
    d2N = len(data[1])
    GT = statistic(data[0], data[1], 0)
    #print(f'd1 shape: {d1N} | d2 shape: {d2N} | GT: {GT} | combdined: {combined.shape}')
    rng = np.random.default_rng()
    pvalue = 1
    for p in range(n_resamples):
        permute_samp = rng.permutation(combined)
        p1 = permute_samp[:d1N]
        p2 = permute_samp[d1N:]

        perm_val = statistic(p1, p2, 0)

        if GT > 0:
            if perm_val >= GT:
                pvalue += 1
        
        if GT < 0:
            if perm_val <= GT:
                pvalue += 1 

    pvalue /= n_resamples 
    
    # Create a namespace object
    results = SimpleNamespace(pvalue=pvalue, statistic=GT)
    
    return results
#######################################################################################

#######################################################################################
def statistic(x, y, axis):
    return np.mean(x, axis=axis) - np.mean(y, axis=axis)
#######################################################################################

#######################################################################################
def statistic_sum(x, y, axis):
    return np.sum(x, axis=axis) - np.sum(y, axis=axis)
#######################################################################################

############## MAIN #################
def main():
    
    script_name = f'CompareMisfoldingMechanism'
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("--candidates", type=str, required=True, help="A file containing two columns. The candidate tag and the groupID")
    parser.add_argument("--CollectedOPpath", type=str, required=True, help="path to the CollectAndProcessesOP DATA/ folder")
    parser.add_argument("--toplevel", type=str, required=True, help="file containing relative paths to either native state GQ files or the MSM file for various proteins")
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
    
    # Setup logging configuration
    logfile = os.path.join(logs, f'{args.outname}.log')
    print(f'logfile: {logfile}')
    logging.basicConfig(filename=logfile, level=logging.INFO, format='%(asctime)s %(message)s')
    logging.info(f'{"#"*50}NEW RUN {script_name}{"#"*50}')


    anal = Analysis(args)
    anal.load_OP()
    anal.plot_lifetimes_cdf()
    stats_df = anal.calc_thresholding_metrics()
    anal.plot_thresholding_metrics(stats_df)

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()

print(f'NORMAL TERMINATION: {time.time() - start_time}')
logging.info(f'NORMAL TERMINATION: {time.time() - start_time}')