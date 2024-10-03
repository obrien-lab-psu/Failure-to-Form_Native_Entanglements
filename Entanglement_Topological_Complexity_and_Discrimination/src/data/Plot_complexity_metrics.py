import sys, os, re, time, logging
from sklearn.utils import shuffle
import ast
from collections import defaultdict
import multiprocessing 
import argparse
import math
import numpy as np
import pandas as pd
from operator import itemgetter
from glob import glob
import pickle
from scipy.spatial import distance
from itertools import product, combinations
from Bio import SeqIO
from Bio import PDB
import requests
import mdtraj as md
from scipy.spatial.distance import pdist, squareform
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import permutation_test, ttest_ind, false_discovery_control, mode, bootstrap, mannwhitneyu
import matplotlib.pyplot as plt
#pd.set_option('display.max_rows', 500)

class DataAnalysis:
    """
    This class is meant to calculate the distribution of various entanglement complexity metrics between Ess and NonEssential gene list

    """

    def __init__(self, outpath, merged_files):
        """
        Initializing the DataAnalysis object and ensure atleast the top level outpath is present and makes it if not. 
        """

        # Make outpath if it doesnt exists
        self.outpath = outpath
        #self.DataAnalysisOutpath = os.path.join(self.outpath, 'DataAnalysisOutput/Plots/')
        self.DataAnalysisOutpath = outpath
        if not os.path.exists(self.DataAnalysisOutpath):
            os.makedirs(self.DataAnalysisOutpath)
            print(f'Made directory: {self.DataAnalysisOutpath}')

        self.merged_files = glob(merged_files)
        for f in self.merged_files:
            print(f)
        self.merged_dfs = {f:pd.read_csv(f, sep='|') for f in self.merged_files}

    #################################################################################################################
    def Plot(self, tag):
        
        #        C  fold  balanced_accuracy  accuracy  
        # Create subplots
        plot_df = []
        for buff_i, buff in enumerate(['C', 'CD', 'CG']):

            merged_df_f = [f for f in self.merged_files if f'_{buff}_' in f][0]
            print(buff_i, buff, merged_df_f)

            merged_df = self.merged_dfs[merged_df_f]
            merged_df['buff'] = buff
            plot_df += [merged_df]

        merged_df = pd.concat(plot_df)
        merged_df['qvalue'] = false_discovery_control(merged_df['pvalues'].values)
        print(f'merged_df:\n{merged_df}')
        outfile = f'{self.DataAnalysisOutpath}merged_plotting_df_{tag}.csv'
        merged_df.to_csv(outfile, sep='|', index=False)
        print(f'SAVED: {outfile}')

        fig, axes = plt.subplots(3, 6, figsize=(14, 6))
        # Flatten the axes array for easy iteration

        ylabels = {'Gn': r'$G_n$',
                    'N_term_thread': r'$N_{Thread}$',
                    'Gc': r'$G_c$',
                    'C_term_thread': r'$C_{Thread}$',
                    'loopsize': r'$N_{Loop}$',
                    'num_zipper_nc': r'$N_{zipper}$',
                    'perc_bb_loop': r'$f_{BB,Loop}$',
                    'num_loop_contacting_res': r'$N_{Loop-cont}$',
                    'num_cross_nearest_neighbors': r'$N_{Cross-cont}$',
                    'ent_coverage': r'$f_{entR}$',
                    'min_N_prot_depth_left': r'$d_P(N)$',
                    'min_N_thread_depth_left': r'$d_T(N)$',
                    'min_N_thread_slippage_left': r'$d_S(N)$',
                    'min_C_prot_depth_right': r'$d_P(C)$',
                    'min_C_thread_depth_right': r'$d_T(C)$',
                    'min_C_thread_slippage_right': r'$d_S(C)$',
                    'ACO': 'ACO',
                    'RCO': 'RCO'}

        axes = axes.flatten()
        for buff_i, buff in enumerate(['C', 'CD', 'CG']):
            buff_df = merged_df[merged_df['buff'] == buff]
            print(f'buff_df:\n{buff_df}')
            
            # Plot each metric
            for i, (idx, row) in enumerate(buff_df.iterrows()):
                metric = row['metric']
                print(metric)
                qvalue = row['qvalue']
                highest_ub = 0
                lowest_lb = 999
                for label in ['Ess', 'NonEss']:
                    stat = row[f'{label}_mean']
                    median = row[f'{label}_median']
                    stat_lb = row[f'{label}_mean_lb']
                    stat_ub = row[f'{label}_mean_ub']
                    if stat_ub > highest_ub:
                        highest_ub = stat_ub
                        highest_lb = stat_lb
                    if stat_lb < lowest_lb:
                        lowest_lb = stat_lb
                    
                    # Plot stat with error bars
                    if label == 'Ess':
                        color = 'blue'
                    if label == 'NonEss':
                        color = 'red'
                    axes[i].errorbar(buff_i, stat, yerr=[[stat - stat_lb], [stat_ub - stat]], fmt='o', capsize=5, color=color, label=label)
                    axes[i].plot(buff_i, median, color=color, marker='X')
                    # Set title and labels
                    #axes[i].set_title(metric)
                    axes[i].set_xlabel(f'Buffer')
                    axes[i].set_ylabel(ylabels[metric])
                    
                    # Add grid
                    #axes[i].grid(True)
                    
                    # Set x-ticks to be empty since it's just a single point
                    axes[i].set_xticks([0, 1, 2], labels=['C', 'CD', 'CG'])
                highest_ub = highest_ub + (highest_ub - highest_lb)*0.1

                axes[i].set_ylim(top=highest_ub + (highest_ub - highest_lb))
                axes[i].set_ylim(bottom=lowest_lb - (highest_ub - highest_lb))
                axes[i].set_xlim(-0.5, 2.5)

                # Add annotation based on p-value
                annotation = ''
                if qvalue < 0.001:
                    annotation = '***'
                elif qvalue < 0.01:
                    annotation = '**'
                elif qvalue < 0.05:
                    annotation = '*'

                if annotation:
                    axes[i].text(buff_i, highest_ub, annotation, ha='center', va='bottom', fontsize=12, color='red')

        axes[5].legend(loc='upper right', bbox_to_anchor=(2, 1))

        plt.tight_layout()
        #plt.show()
        outfile = f'{self.DataAnalysisOutpath}Complexity_measures_{tag}.png'

        plt.savefig(outfile)
        print(f'SAVED: {outfile}')

#################################################################################################################

def main():
    """
    Main function to control workflow. 
    """

    # Parse the user supplied arguments
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-m", "--merged_files", type=str, required=True, help=f"path to merged stats file")
    parser.add_argument("-l", "--log_file", type=str, required=True, help="Path to logging file")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="path to output directory. will be made if doesnt exist")
    #parser.add_argument("-s", "--stat_type", type=str, required=True, help="mean or median to use in plots")
    parser.add_argument("-t", "--tag", type=str, required=True, help="tag for output file")
    args = parser.parse_args()

    log_file = args.log_file
    outpath = args.outpath
    merged_files = args.merged_files
    #stat_type = args.stat_type
    tag = args.tag

    # Setup logging configuration
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s %(message)s') 
    logging.info(f'{"#"*100}\nNEW RUN')

    # Make outpath if it doesnt exists
    if not os.path.exists(outpath):
        os.makedirs(outpath)
        print(f'Made directory: {outpath}')

    # Initalize the DataAnalysis class object
    Analyzer = DataAnalysis(outpath, merged_files)
    print(f'Analyzer: {Analyzer}')

    Analyzer.Plot(tag)
if __name__ == "__main__":
    main()

print('NORMAL TERMINATION')
