import time, sys
import pickle
import multiprocessing as mp
from scipy.stats import bootstrap
import logging
import argparse
import glob
import numpy as np
import pandas as pd
from sklearn.preprocessing import OneHotEncoder, LabelEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.model_selection import train_test_split, StratifiedKFold, KFold, cross_validate, GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, average_precision_score, f1_score, recall_score, precision_score, roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from scipy.spatial.distance import euclidean
import matplotlib.pyplot as plt
import os
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import poisson, binom, fisher_exact, chi2, norm, kstest
import scipy.stats as st
from matplotlib.ticker import MultipleLocator
import math, random
import itertools
from tqdm import tqdm

pd.set_option('display.max_rows', 4000)
pd.options.mode.chained_assignment = None  # default='warn'

class DataAnalysis:
    """
    A class to handle the data analysis process including encoding, regression, and statistical tests.
    """
    ##########################################################################################################
    def __init__(self, inpdir, outpath, tag):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:

        """
        self.inpdir = inpdir
        self.outpath = outpath
        self.tag = tag

        if not os.path.exists(f'{self.outpath}'):
            os.makedirs(f'{self.outpath}')
            print(f'Made output directories {self.outpath}')
    ##########################################################################################################

    ##########################################################################################################
    def load_data(self, buff):
        """
        Loads the residue feature files and filters the data for analysis.
        """
        
        # get all files in inpdir
        files = []
        for dirpath, dirnames, filenames in os.walk(self.inpdir):
            for filename in filenames:
                if f'_{buff}_' in filename:
                    files += [os.path.join(dirpath, filename)]
        print(f'Number of files found for {buff}: {len(files)}')

        ## parse the datafiles for final gene lists
        final_genelist = {}
        final_traj = {}
        Final_step_reg = {}
        for f in files:
            traj = f.split('/')[-2]

            if 'final_genelist' in f:
                state = f.split('/')[-1].split('_')[0].replace('State','')
                if traj not in final_genelist:
                    final_genelist[traj] = {}
                final_genelist[traj][state] = np.loadtxt(f, dtype=str)

            if 'final_traj' in f:
                state = f.split('/')[-1].split('_')[0].replace('State','')
                if traj not in final_traj:
                    final_traj[traj] = {}
                df = pd.read_csv(f)
                final_traj[traj][state] = pd.read_csv(f)

            if 'Final_step_reg' in f:
                if traj not in Final_step_reg:
                    Final_step_reg[traj] = {}
                Final_step_reg[traj] = pd.read_csv(f)

        return final_genelist, final_traj, Final_step_reg
    ##########################################################################################################

    ##########################################################################################################
    def Plot_traj(self, data, title):
        ## Data is a tuple with three elements: final_genelist, final_traj, Final_step_reg
        traj_data = data[1]
        final_reg_data = data[2]

        for traj, traj_states in traj_data.items():
            fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(6, 6))
            outfile = os.path.join(self.outpath, f'Traj{traj}_{title}.png')
            for state, state_df in traj_states.items():
                #print(state, state_df)

                x = state_df['step'].values
                OR = state_df['OR'].values
                pvals = state_df['pvalue'].values
                cuts = state_df['cuts'].values
                size_means = state_df['psize_mean'].values

                # plot OR
                axes[0].plot(x, OR, label=f't{traj}s{state}', ls='-')
                axes[0].axhline(y=1, color='red', linestyle='--')
                axes[0].set_ylabel('OR')
                axes[0].set_xlabel('MC step')
                axes[0].legend(loc='upper left', bbox_to_anchor=(1,1))

                # plot pvals
                axes[1].plot(x, pvals, label=f't{traj}s{state}', ls='-')
                axes[1].axhline(y=0.05, color='red', linestyle='--')
                axes[1].set_yscale('log')
                axes[1].set_ylabel('pvals')
                axes[1].set_xlabel('MC step')

                # plot cuts
                axes[2].plot(x, cuts, label=f't{traj}s{state}', ls='-')
                axes[2].set_ylabel('cuts')
                axes[2].set_xlabel('MC step')
                #axes[2].set_ylim(50,500)
                axes[2].set_yticks(range(50, 550, 150))
                axes[2].grid(True, which='both', linestyle='--', linewidth=0.5)

                # plot size_means
                axes[3].plot(x, size_means, label=f't{traj}s{state}', ls='-')
                axes[3].set_ylabel('prot. size')
                axes[3].set_xlabel('MC step')

            plt.suptitle(title)
            plt.tight_layout()
            #plt.show()
            plt.savefig(outfile)
            plt.close()
            print(f'SAVED: {outfile}')
    ##########################################################################################################

    ##########################################################################################################
    def State_overlap(self, data, title):
        ## Data is a tuple with three elements: final_genelist, final_traj, Final_step_reg
        final_genelist_data = data[0]
        traj_data = data[1]
        final_reg_data = data[2]

        ## identify the highest OR state from each traj
        max_state_genes = {}
        total_genes = []
        for traj, traj_states in traj_data.items():
            max_state = 0
            max_OR = 0
            for state, state_df in traj_states.items():
                print(state, state_df)
                final_OR = state_df['OR'].values[-1]
                if final_OR > max_OR:
                    max_OR = final_OR
                    max_state = state
            
            print(f'traj {traj} with max_state: {max_state} with OR {max_OR}')
            max_state_genes[traj] = final_genelist_data[traj][max_state]
            print(max_state_genes[traj])
            total_genes += [max_state_genes[traj]]

        total_genes = np.hstack(total_genes)
        total_genes = set(total_genes)
        print(f'total_genes: {total_genes} {len(total_genes)}')

        ## for each gene in total_genes calculate the fraction of traj it was present in the top state
        dist = []
        for gene in total_genes:
            F = 0
            for traj, traj_genes in max_state_genes.items():
                if gene in traj_genes:
                    F += 1
            F /= 10
            print(gene, F)
            dist += [F]
        print(f'dist: {dist}')

        ## Plot the data
        # Create a figure with two subplots in the same row
        fig, ax = plt.subplots(1, 2, figsize=(12, 5))
        outfile = os.path.join(self.outpath, f'Fract_High-OR_dists_{title}.png')

        # First subplot: Normalized Histogram
        ax[0].hist(dist, bins=30, density=False, alpha=0.7, color='blue', width=0.05)
        ax[0].set_title('Normalized Histogram')
        ax[0].set_xlabel('Fraction traj in High-OR')
        ax[0].set_ylabel('Frequency')

        # Second subplot: CDF
        sorted_data = np.sort(dist)
        cdf = np.arange(1, len(sorted_data) + 1) / len(sorted_data)
        
        ax[1].bar(sorted_data, cdf, width=0.05, color='red', alpha=0.7)
        ax[1].set_title('Cumulative Distribution Function (CDF)')
        ax[1].set_xlabel('Fraction traj in High-OR')
        ax[1].set_ylabel('CDF')

        # Display the plots
        plt.suptitle(title)
        plt.tight_layout()
        #plt.show()
        plt.savefig(outfile)
        plt.close()
        print(f'SAVED: {outfile}')
    ##########################################################################################################

    ##########################################################################################################
    def Plot_ensem(self, data, title):
        ## Data is a tuple with three elements: final_genelist, final_traj, Final_step_reg
        traj_data = data[1]
        final_reg_data = data[2]

        ## determine rank orders of states
        ensemble_state_idxs = {s:[] for s in range(5)}
        print(ensemble_state_idxs)

        for traj, traj_states in traj_data.items():
            traj_ranking = []
            for state, state_df in traj_states.items():
                final_OR = state_df['OR'].values[-1]
                traj_ranking += [[state, final_OR]]
            traj_ranking = np.asarray(traj_ranking)
            traj_ranking = traj_ranking[traj_ranking[:,1].argsort()[::-1]]

            for s, (state, OR) in enumerate(traj_ranking):
                print(traj, s, state, OR)
                ensemble_state_idxs[s] += [traj_states[state]]

        # quality checck that each ensemble state has 10 dfs
        for s, s_dfs in ensemble_state_idxs.items():
            if len(s_dfs) != 10:
                raise ValueError('Ensemble state does not have 10 traj in it')
        
 
        fig, axes = plt.subplots(nrows=4, ncols=5, figsize=(12, 10), sharex=True)
        outfile = os.path.join(self.outpath, f'Ensemble_Traj_{title}.png')
        for s_idx, (s, s_dfs) in enumerate(ensemble_state_idxs.items()):
            
            # get the average across all dataframe in state
            means_df, ci_lower_df, ci_upper_df = running_average_and_ci(s_dfs, ['step', 'OR', 'cuts', 'pvalue', 'psize_mean'])
            print(means_df)
            print(ci_lower_df)
            print(ci_upper_df)

            x = means_df['step'].values

            OR = means_df['OR'].values
            OR_lb = OR - ci_lower_df['OR']
            OR_ub = ci_upper_df['OR'] - OR

            pvals = means_df['pvalue'].values
            pvals_lb = pvals - ci_lower_df['pvalue']
            pvals_ub = ci_upper_df['pvalue'] - pvals

            cuts = means_df['cuts'].values
            cuts_lb = cuts - ci_lower_df['cuts']
            cuts_ub = ci_upper_df['cuts'] - cuts

            size_means = means_df['psize_mean'].values
            size_means_lb = size_means - ci_lower_df['psize_mean']
            size_means_ub = ci_upper_df['psize_mean'] - size_means

            # plot OR
            axes[0, s_idx].errorbar(x, OR, yerr=[OR_lb, OR_ub], label=f'state-{s}', marker='o', ls='none', fillstyle='none', capsize=3)
            axes[0, s_idx].plot(x, OR, label=f'state-{s}', ls='-', color='k', zorder=10)
            axes[0, s_idx].axhline(y=1, color='red', linestyle='--')
            axes[0, s_idx].set_ylabel('OR')
            #axes[0, s_idx].set_xlabel('MC step')
            axes[0, s_idx].set_title(f'State: {s}')


            # plot pvals
            #axes[1].errorbar(x, pvals, yerr=[pvals_lb, pvals_ub], label=f'state-{s}', marker='o', ls='none', fillstyle='none', capsize=3)
            axes[1, s_idx].plot(x, pvals, label=f'state-{s}', ls='-')
            axes[1, s_idx].axhline(y=0.05, color='red', linestyle='--')
            axes[1, s_idx].set_yscale('log')
            axes[1, s_idx].set_ylabel('pvals')
            #axes[1, s_idx].set_xlabel('MC step')

            # plot cuts
            axes[2, s_idx].errorbar(x, cuts, yerr=[cuts_lb, cuts_ub], label=f'state-{s}', marker='o', ls='none', fillstyle='none', capsize=3)
            axes[2, s_idx].plot(x, cuts, label=f'state-{s}', ls='-', color='k', zorder=10)
            axes[2, s_idx].set_ylabel('cuts')
            #axes[2, s_idx].set_xlabel('MC step')
            #axes[2].set_ylim(50,500)
            axes[2, s_idx].set_yticks(range(50, 550, 150))
            axes[2, s_idx].grid(True, which='both', linestyle='--', linewidth=0.5)

            # plot size_means
            axes[3, s_idx].errorbar(x, size_means, yerr=[size_means_lb, size_means_ub], label=f'state-{s}', marker='o', ls='none', fillstyle='none', capsize=3)
            axes[3, s_idx].plot(x, size_means, label=f'state-{s}', ls='-', color='k', zorder=10)
            axes[3, s_idx].set_ylabel('prot. size')
            axes[3, s_idx].set_xlabel('MC step')

        for row in range(0,4):
            for col in range(0,5):
                axes[row, col].sharey(axes[row, 0])

        #axes[0, 4].legend(loc='upper left', bbox_to_anchor=(1,1))
        plt.suptitle(title)
        plt.tight_layout()
        #plt.show()
        plt.savefig(outfile)
        plt.close()
        print(f'SAVED: {outfile}')
    ##########################################################################################################


    ##########################################################################################################
    def run(self):
        """
        Orchestrates the workflow by loading data, performing regression, and saving results.
        """
        start_time = time.time()

        # Load data for cyto-serum only
        C_data = self.load_data('C')

        # Plot the trajectories of each state and the final state regression data
        #self.Plot_traj(C_data, f'{self.tag}_C')

        # Get the overlap of gene histograms and CDF across the highest OR state in each traj
        self.State_overlap(C_data, f'{self.tag}_C')

        # plot the enesmble of traj in a single plot
        self.Plot_ensem(C_data, f'Ensemble_{self.tag}_C')

        logging.info(f'NORMAL TERMINATION {time.time() - start_time}')
    ##########################################################################################################

def running_average_and_ci(df_list, columns):
    # Stack the dataframes vertically for easy computation across rows
    stacked_data = np.stack([df[columns].values for df in df_list], axis=2)

    # Calculate the mean across the third dimension (which contains the dataframes)
    means = np.mean(stacked_data, axis=2)

    # Calculate the standard deviation across the third dimension
    std_devs = np.std(stacked_data, axis=2, ddof=1)  # ddof=1 for sample standard deviation

    # Number of dataframes
    n = len(df_list)

    # Calculate the 95% confidence interval
    margin_of_error = 1.96 * (std_devs / np.sqrt(n))
    lower_bound = means - margin_of_error
    upper_bound = means + margin_of_error

    # Convert to DataFrames for ease of use
    means_df = pd.DataFrame(means, columns=columns)
    ci_lower_df = pd.DataFrame(lower_bound, columns=columns)
    ci_upper_df = pd.DataFrame(upper_bound, columns=columns)

    return means_df, ci_lower_df, ci_upper_df

def main():
    """
    Main function to parse arguments and run the DataAnalysis class.
    """
    script_name = f'Optimizer_SimulatedAnnealing'
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-i", "--inpdir", type=str, required=True, help="Path to top level directory containing trajectories")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("-t", "--tag", type=str, required=True, help="Tag for output filenames")
    parser.add_argument("-l", "--log", type=str, required=True, help="Path to logging file")
    args = parser.parse_args()

    # Setup logging configuration
    logging.basicConfig(filename=args.log, level=logging.INFO, format='%(asctime)s %(message)s')
    logging.info(f'{"#"*100}\nNEW RUN {script_name}')

    analysis = DataAnalysis(
        inpdir=args.inpdir,
        outpath=args.outpath,
        tag=args.tag)
    analysis.run()

if __name__ == "__main__":
    main()

