import sys, os, re, time, logging
import scipy.stats as st
from scipy.stats import bootstrap, t
from scipy.optimize import curve_fit
from sklearn.utils import shuffle
import ast
from collections import defaultdict
from multiprocessing import Pool, cpu_count
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
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import permutation_test, ttest_ind, false_discovery_control, fisher_exact, mannwhitneyu
import matplotlib.pyplot as plt
from math import comb
pd.set_option('display.max_rows', 500)

class Analyzer:
    """
Test 1: What is the fraction of proteins with atleast 1 SLFC and is that fraction higher in non-essential proteins than essential proteins
        i. for the set of essential and non-essential genes calculate the fraction of proteins with atleast 1 SLFC and then do a permuration test
        ii. get the distributions of fraction LFC for essential and non-essential and do permutation tests to determine with they are significantly different

Test 2: Is the OR (misfolding Y/N, entangled region Y/N) dependant on the number of SLFC at a given SPA threshold
        i. Make the following groups of proteins
                (1) has atleast n SLFC (n = 0, 1, 2, 3, 4, 5)
                (2) has n SLFC (n = 0, 1, 2, 3, 4, 5)
        ii. For each group calculate the OR using the contingency table and estimate bootstrapped 95% ci
        iii. calculate the delta between n=0 and the other n values for each group and determine if the slope is signifcantly different from zero
        iv. permutation test for differences between n=0 and other n valued states
        v. Chi-squared test for association (not recommended)
        vi. log-linear model "count ~ cut * region * n"
    """

    def __init__(self, outpath, tag, buff, resFeat_files, spa, LiPMScov, gene_list, ess_gene_list, noness_gene_list, loopcontact_df):
        """
        Initializing the FrequencyGenerator object and ensure atleast the top level outpath is present and makes it if not.
        """

        # Make outpath if it doesnt exists
        self.Outpath = outpath
        if not os.path.exists(self.Outpath):
            os.makedirs(self.Outpath)
            print(f'Made directory: {self.Outpath}')

        self.resFeatsOutpath = os.path.join(outpath, 'resFeats/')
        if not os.path.exists(self.resFeatsOutpath):
            os.makedirs(self.resFeatsOutpath)
            print(f'Made directory: {self.resFeatsOutpath}')

        self.resFeat_files = glob(os.path.join(resFeat_files, '*'))

        self.tag = tag
        self.buff = buff
        self.spa = spa
        self.LiPMScov = LiPMScov
        print(f'Tag: {self.tag} | Buff: {self.buff} | SPA: {self.spa} | LiPMScov: {self.LiPMScov}')

        ## Load gene lists
        self.gene_list = np.loadtxt(gene_list, dtype=str)
        print(f'Loaded gene_list: {len(self.gene_list)}')
        self.ess_gene_list = np.loadtxt(ess_gene_list, dtype=str)
        print(f'Loaded ess_gene_list: {len(self.ess_gene_list)}')
        self.noness_gene_list = np.loadtxt(noness_gene_list, dtype=str)
        print(f'Loaded noness_gene_list: {len(self.noness_gene_list)}')


        ## Load the residue feature files 
        dfs = []
        for gene in self.gene_list:
            f = [f for f in self.resFeat_files if gene in f]
            if len(f) != 1:
                raise ValueError(f"No residue feature file for gene {gene}")
            else:
                f = f[0]
            df = pd.read_csv(f, sep='|')
            dfs += [df]
        self.resFeat_df = pd.concat(dfs)
        print(f'resFeat_df:\n{self.resFeat_df}')


        ## Load loop forming contact classification df
        self.loopcontact_df = pd.read_csv(loopcontact_df, sep='|')
        print(f'loopcontact_df:\n{self.loopcontact_df}')
    ##################################################################

    ##################################################################
    def FractionGenes_w_n_LoopContacts(self,):
        """
        What is the fraction of proteins with atleast 1 SLFC and is that fraction higher in non-essential proteins than essential proteins
            i. for the set of essential and non-essential genes calculate the fraction of proteins with atleast 1 SLFC and then do a permuration test
        """
        
        ## make directory to output test 1 info to 
        self.FractionGenes_w_LoopContactsOutpath = os.path.join(self.Outpath, 'FractionGenes_w_LoopContacts/')
        if not os.path.exists(self.FractionGenes_w_LoopContactsOutpath):
            os.makedirs(self.FractionGenes_w_LoopContactsOutpath)
            print(f'Made directory: {self.FractionGenes_w_LoopContactsOutpath}')

        ## for each of the loop contacting types in the self.loopcontact_df object get the fraction of essential proteins with atleast n of those contact where n = 1 2 3
        # List of columns to exclude from the operation
        outfile = f'{self.FractionGenes_w_LoopContactsOutpath}EssVSNonEss_stats_{self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}.csv'
        if not os.path.exists(outfile):
            exclude_columns = ['gene', 'Total', 'Essential']
            outdf = []
            for n in [1, 2, 3, 4, 5]:
                n_loopcontact_df = self.loopcontact_df.copy()

                # Apply the condition to all columns except those in exclude_columns
                for col in n_loopcontact_df.columns:
                    if col not in exclude_columns:
                        n_loopcontact_df[col] = n_loopcontact_df[col] >= n
                
                Ess_loopcontact_df = n_loopcontact_df[n_loopcontact_df['gene'].isin(self.ess_gene_list)]
                Ess_loopcontact_df['Essential'] = True
                NonEss_loopcontact_df = n_loopcontact_df[n_loopcontact_df['gene'].isin(self.noness_gene_list)]
                NonEss_loopcontact_df['Essential'] = False
                n_loopcontact_df = pd.concat([Ess_loopcontact_df, NonEss_loopcontact_df], ignore_index=True)
                #print(n_loopcontact_df)

                # Get the stats results for this n threshold
                stats = []
                for col in n_loopcontact_df.columns:
                    if col not in exclude_columns:
                        stats += [self.analyze_groups(n_loopcontact_df, metric_col=col, group_col='Essential', stat='mean')]

                # merge into a single dictionary and make a df
                n_stats_df = defaultdict(list)
                for d in stats:
                    for k,v in d.items():
                        n_stats_df[k] += [v]
                n_stats_df = pd.DataFrame(n_stats_df)
                n_stats_df['n'] = n
                n_stats_df['Essential True n'] = len(Ess_loopcontact_df)
                n_stats_df['Essential False n'] = len(NonEss_loopcontact_df)
                #print(n_stats_df)
                outdf += [n_stats_df]

            outdf = pd.concat(outdf, ignore_index=True)
            outfile = f'{self.FractionGenes_w_LoopContactsOutpath}EssVSNonEss_stats_{self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}.csv'
            outdf.to_csv(outfile, sep='|', index=False)
            print(f'SAVED: {outfile}')

        else:
            outdf = pd.read_csv(outfile, sep='|')
            print(f'LOADED: {outfile}')
        print(outdf)

        ##############################################
        ## Plot each n level 
        #outdf keys: ['Metric', 'Essential True Mean', 'Essential True 95% CI low',
        #'Essential True 95% CI high', 'Essential False Mean',
        #'Essential False 95% CI low', 'Essential False 95% CI high',
        #'Permutation Test p-value', 'n', 'Essential True n',
        #'Essential False n']
        for n, n_df in outdf.groupby('n'):
            outfile = f'{self.FractionGenes_w_LoopContactsOutpath}EssVSNonEss_stats_atleast_n{n}_LFC_{self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}.png'
            outfile_csv = f'{self.FractionGenes_w_LoopContactsOutpath}EssVSNonEss_stats_atleast_n{n}_LFC_{self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}.csv'

            n1_df = n_df[['Metric', 'Essential True Mean', 'Essential True 95% CI low', 'Essential True 95% CI high', 'Permutation Test p-value', 'Essential True n']] 
            n2_df = n_df[['Metric', 'Essential False Mean', 'Essential False 95% CI low', 'Essential False 95% CI high', 'Permutation Test p-value', 'Essential False n']] 
            print(n1_df)
            print(n2_df)

            # Plot the data and the regression line
            fig,ax = plt.subplots(4,1, figsize=(6,6), sharex=True)
            
            x = n_df['Metric']
            y1 = n_df['Essential True Mean'].values
            y1err = [y1 - n_df['Essential True 95% CI low'], n_df['Essential True 95% CI high'] - y1]
            y2 = n_df['Essential False Mean'].values
            y2err = [y2 - n_df['Essential False 95% CI low'], n_df['Essential False 95% CI high'] - y2]
            
            diff_y1y2 = []
            diff_y1y2_err = []
            for row in n_df[['Essential True Mean', 'Essential True 95% CI low', 'Essential True 95% CI high', 'Essential False Mean', 'Essential False 95% CI low', 'Essential False 95% CI high']].values:
                diff, err, ci = self.difference_with_error_propagation(row[0], (row[1], row[2]), row[3], (row[4], row[5]))
                diff_y1y2 += [diff]
                diff_y1y2_err += [err]

            pvalues = n_df['Permutation Test p-value'].values
            pvalues = false_discovery_control(pvalues) 
            y1_size = n_df['Essential True n']
            y2_size = n_df['Essential False n']

            ## make the plot dataframe and save it
            plot_df = {'Metric':x, 'Essential':y1, 'Essential_lb':y1err[0], 'Essential_ub':y1err[1], 'Essential_n':y1_size, 
                        'NonEssential':y2, 'NonEssential_lb':y2err[0], 'NonEssential_ub':y2err[1], 'NonEssential_n':y2_size, 'pvalues':pvalues}
            plot_df = pd.DataFrame(plot_df) 
            print(f'plot_df:\n{plot_df}')
            plot_df.to_csv(outfile_csv)
            logging.info(f'SAVED: {outfile_csv}')

            ax[0].errorbar(x, y1, yerr=y1err, label='Essential', capsize=3, marker='o', ls='None', markerfacecolor='none')
            ax[0].errorbar(x, y2, yerr=y2err, label='NonEssential', capsize=3, marker='o', ls='None', markerfacecolor='none')
            ax[0].set_ylabel(f'Fraction of genes with {n}\nloop forming contacts', fontsize=7)
            ax[0].legend(loc='upper left', bbox_to_anchor=(1., 1))
            ax[0].tick_params(axis='both', labelsize=7)

            ax[1].errorbar(x, diff_y1y2, yerr=diff_y1y2_err, capsize=3, marker='o', ls='None', markerfacecolor='none')
            ax[1].set_ylabel(f'Diff', fontsize=7)
            ax[1].axhline(y=0, color='red', linestyle='--', linewidth=2)
            ax[1].tick_params(axis='both', labelsize=7)

            ax[2].scatter(x, pvalues)
            ax[2].set_yscale('log')
            ax[2].axhline(y=0.05, color='red', linestyle='--', linewidth=2, label='y = 0.05')
            ax[2].set_ylabel('Permute pvalue (FDR)', fontsize=7)
            ax[2].tick_params(axis='both', labelsize=7)
            
            ax[3].scatter(x, y1_size, label='Essential')
            ax[3].scatter(x, y2_size, label='NonEssential')
            ax[3].set_ylabel('Sample size', fontsize=7)
            ax[3].tick_params(axis='both', labelsize=7)
            #ax[3].set_yscale('log')

            # R3tate the x-axis tick labels by 45 degrees for all subplots
            plt.setp(ax[3].get_xticklabels(), rotation=45, ha="right", fontsize=7)            

            #plt.title(f'')
            plt.tight_layout()
            #plt.grid(True)
            
            # Show the plot
            #plt.show()
            plt.savefig(outfile)
            logging.info(f'SAVED: {outfile}')
            plt.close()
    ##################################################################


    ##############################################################################################################
    def analyze_groups(self, df, metric_col='A', group_col='B', n_resamples=10000, confidence_level=0.95, stat='mean'):
        """
        (1) Calculates the mean of A for each group defined by B.
        (2) Calculates the 95% confidence intervals for each group's mean using bootstrapping.
        (3) Performs a permutation test to evaluate the difference in means between the two groups.
        """

        # Perform permutation test for the difference in means
        def mean_diff(x, y):
            return np.mean(x) - np.mean(y)
        def median_diff(x, y):
            return np.median(x) - np.median(y)

        # Separate data into two groups based on the binary flag in column B
        group_1 = df[df[group_col] == True][metric_col].values
        group_2 = df[df[group_col] == False][metric_col].values
        
        if stat == 'mean':
            # Calculate means for each group
            mean_1 = np.mean(group_1)
            mean_2 = np.mean(group_2)
            
            # Bootstrapping to calculate 95% confidence intervals
            ci_1 = bootstrap((group_1,), np.mean, confidence_level=confidence_level, n_resamples=n_resamples).confidence_interval
            ci_2 = bootstrap((group_2,), np.mean, confidence_level=confidence_level, n_resamples=n_resamples).confidence_interval
            
            perm_test_result = permutation_test((group_1, group_2), mean_diff, alternative='two-sided', n_resamples=n_resamples)

            # Compile results
            results = {
                f'Metric': metric_col,
                f'{group_col} True Mean': mean_1,
                f'{group_col} True 95% CI low': ci_1.low,
                f'{group_col} True 95% CI high': ci_1.high,
                f'{group_col} False Mean': mean_2,
                f'{group_col} False 95% CI low': ci_2.low,
                f'{group_col} False 95% CI high': ci_2.high,
                'Permutation Test p-value': perm_test_result.pvalue}

        elif stat == 'median':
            # Calculate medians for each group
            median_1 = np.median(group_1)
            median_2 = np.median(group_2)
            
            # Bootstrapping to calculate 95% confidence intervals
            ci_1 = bootstrap((group_1,), np.median, confidence_level=confidence_level, n_resamples=n_resamples).confidence_interval
            ci_2 = bootstrap((group_2,), np.median, confidence_level=confidence_level, n_resamples=n_resamples).confidence_interval
            
            perm_test_result = permutation_test((group_1, group_2), median_diff, alternative='two-sided', n_resamples=n_resamples)

            # Compile results
            results = {
                f'Metric': metric_col,
                f'{group_col} True Median': median_1,
                f'{group_col} True 95% CI low': ci_1.low,
                f'{group_col} True 95% CI high': ci_1.high,
                f'{group_col} False Median': median_2,
                f'{group_col} False 95% CI low': ci_2.low,
                f'{group_col} False 95% CI high': ci_2.high,
                'Permutation Test p-value': perm_test_result.pvalue}

        else:
            raise ValueError("Stat option not mean or median")

        return results
    ##################################################################

    ####################################################################################################
    def LoopContactsDistsCompare(self,):

        """
        Get the distributions of fraction LFC for essential and non-essential and do permutation tests to determine with they are significantly different
        """
        
        ## make directory to output test 1 info to 
        self.LoopContactsDistsCompareOutpath = os.path.join(self.Outpath, 'LoopContactsDistsCompare/')
        if not os.path.exists(self.LoopContactsDistsCompareOutpath):
            os.makedirs(self.LoopContactsDistsCompareOutpath)
            print(f'Made directory: {self.LoopContactsDistsCompareOutpath}')

        ## for each of the loop contacting types in the self.loopcontact_df object get the fraction of essential proteins with atleast n of those contact where n = 1 2 3
        # List of columns to exclude from the operation
        outfile = f'{self.LoopContactsDistsCompareOutpath}EssVSNonEss_dist_means_stats_{self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}.csv'
        if not os.path.exists(outfile):
            exclude_columns = ['gene', 'Total', 'Essential']
            df_copy = self.loopcontact_df.copy()
                
            Ess_loopcontact_df = df_copy[df_copy['gene'].isin(self.ess_gene_list)]
            Ess_loopcontact_df['Essential'] = True
            NonEss_loopcontact_df = df_copy[df_copy['gene'].isin(self.noness_gene_list)]
            NonEss_loopcontact_df['Essential'] = False
            df_copy = pd.concat([Ess_loopcontact_df, NonEss_loopcontact_df], ignore_index=True)
            print(df_copy)


            ## Compare means
            # Get the stats results for this n threshold
            stats = []
            for col in df_copy.columns:
                if col not in exclude_columns:
                    stats += [self.analyze_groups(df_copy, metric_col=col, group_col='Essential', stat='mean')]

            # merge into a single dictionary and make a df
            stats_df = defaultdict(list)
            for d in stats:
                for k,v in d.items():
                    stats_df[k] += [v]
            stats_df = pd.DataFrame(stats_df)
            stats_df['Essential True n'] = len(Ess_loopcontact_df)
            stats_df['Essential False n'] = len(NonEss_loopcontact_df)
            print(stats_df)

            outfile = f'{self.LoopContactsDistsCompareOutpath}EssVSNonEss_dist_means_stats_{self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}.csv'
            stats_df.to_csv(outfile, sep='|', index=False)
            print(f'SAVED: {outfile}')
        else:
            stats_df = pd.read_csv(outfile, sep='|')
            print(f'LOADED: {outfile}')
            print(stats_df)


        ##############################################
        ## Plot each n level 
        #outdf keys: ['Metric', 'Essential True Mean', 'Essential True 95% CI low',
        #'Essential True 95% CI high', 'Essential False Mean',
        #'Essential False 95% CI low', 'Essential False 95% CI high',
        #'Permutation Test p-value', 'n', 'Essential True n',
        #'Essential False n']
        outfile = f'{self.LoopContactsDistsCompareOutpath}EssVSNonEss_dist_stats_LFC_{self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}.png'

        n1_df = stats_df[['Metric', 'Essential True Mean', 'Essential True 95% CI low', 'Essential True 95% CI high', 'Permutation Test p-value', 'Essential True n']] 
        n2_df = stats_df[['Metric', 'Essential False Mean', 'Essential False 95% CI low', 'Essential False 95% CI high', 'Permutation Test p-value', 'Essential False n']] 
        print(n1_df)
        print(n2_df)
        # Plot the data and the regression line
        fig,ax = plt.subplots(4,1, figsize=(6,6), sharex=True)
        
        x = stats_df['Metric']
        y1 = stats_df['Essential True Mean'].values
        y1err = [y1 - stats_df['Essential True 95% CI low'], stats_df['Essential True 95% CI high'] - y1]
        y2 = stats_df['Essential False Mean'].values
        y2err = [y2 - stats_df['Essential False 95% CI low'], stats_df['Essential False 95% CI high'] - y2]
        
        diff_y1y2 = []
        diff_y1y2_err = []
        for row in stats_df[['Essential True Mean', 'Essential True 95% CI low', 'Essential True 95% CI high', 'Essential False Mean', 'Essential False 95% CI low', 'Essential False 95% CI high']].values:
            diff, err, ci = self.difference_with_error_propagation(row[0], (row[1], row[2]), row[3], (row[4], row[5]))
            diff_y1y2 += [diff]
            diff_y1y2_err += [err]

        pvalues = stats_df['Permutation Test p-value'].values
        pvalues = false_discovery_control(pvalues) 
        y1_size = stats_df['Essential True n']
        y2_size = stats_df['Essential False n']

        ax[0].errorbar(x, y1, yerr=y1err, label='Essential', capsize=3, marker='o', ls='None', markerfacecolor='none')
        ax[0].errorbar(x, y2, yerr=y2err, label='NonEssential', capsize=3, marker='o', ls='None', markerfacecolor='none')
        ax[0].set_ylabel(f'n', fontsize=7)
        ax[0].legend(loc='upper left', bbox_to_anchor=(1., 1))
        ax[0].tick_params(axis='both', labelsize=7)

        ax[1].errorbar(x, diff_y1y2, yerr=diff_y1y2_err, capsize=3, marker='o', ls='None', markerfacecolor='none')
        ax[1].set_ylabel(f'Diff', fontsize=7)
        ax[1].axhline(y=0, color='red', linestyle='--', linewidth=2)
        ax[1].tick_params(axis='both', labelsize=7)

        ax[2].scatter(x, pvalues)
        ax[2].set_yscale('log')
        ax[2].axhline(y=0.05, color='red', linestyle='--', linewidth=2, label='y = 0.05')
        ax[2].set_ylabel('Permute pvalue (FDR)', fontsize=7)
        ax[2].tick_params(axis='both', labelsize=7)
        
        ax[3].scatter(x, y1_size, label='Essential')
        ax[3].scatter(x, y2_size, label='NonEssential')
        ax[3].set_ylabel('Sample size', fontsize=7)
        ax[3].tick_params(axis='both', labelsize=7)
        ax[3].set_ylim(10,1000)
        ax[3].set_yscale('log')

        # R3tate the x-axis tick labels by 45 degrees for all subplots
        plt.setp(ax[3].get_xticklabels(), rotation=45, ha="right", fontsize=7)            

        #plt.title(f'')
        plt.tight_layout()
        #plt.grid(True)
        
        # Show the plot
        #plt.show()
        plt.savefig(outfile)
        plt.close()

    ##################################################################

    ##################################################################
    def OR_trends(self,):

        """
        Test 2: Is the OR (misfolding Y/N, entangled region Y/N) dependant on the number of SLFC at a given SPA threshold
                i. Make the following groups of proteins
                        (1) has atleast n SLFC (n = 0, 1, 2, 3, 4, 5)
                        (2) has n SLFC (n = 0, 1, 2, 3, 4, 5)
                ii. For each group calculate the OR using the contingency table and estimate bootstrapped 95% ci
                iii. calculate the delta between n=0 and the other n values for each group and determine if the slope is signifcantly different from zero
        """

        ## make directory to output test 1 info to 
        self.OR_trendsOutpath = os.path.join(self.Outpath, 'OR_trends/')
        if not os.path.exists(self.OR_trendsOutpath):
            os.makedirs(self.OR_trendsOutpath)
            print(f'Made directory: {self.OR_trendsOutpath}')

        
        ############################################################################
        ## for each of the loop contacting types in the self.loopcontact_df object get the list of genes with n or => n loop forming contacts
        # List of columns to exclude from the operation
        exclude_columns = ['gene', 'Total', 'Essential']
        atleast_n_loopcontacts = {}
        n_loopcontacts = {}
        n_loopcontacts_gene_list_outfile = os.path.join(self.OR_trendsOutpath, f'n_loopcontacts_gene_lists_{self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}.pkl')
        atleast_n_loopcontacts_gene_list_outfile = os.path.join(self.OR_trendsOutpath, f'atleast_n_loopcontacts_gene_lists_{self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}.pkl')
        for n in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:

            n_loopcontact_df = self.loopcontact_df[self.loopcontact_df['gene'].isin(self.gene_list)].copy()
            n_loopcontacts[n] = {}

            atleast_n_loopcontact_df = self.loopcontact_df[self.loopcontact_df['gene'].isin(self.gene_list)].copy()
            atleast_n_loopcontacts[n] = {}

            # Apply the condition to all columns except those in exclude_columns
            for col in n_loopcontact_df.columns:
                if col not in exclude_columns:
                    n_loopcontact_df[col] = n_loopcontact_df[col] == n
                    n_loopcontact_df_genes = n_loopcontact_df[n_loopcontact_df[col] == True]
                    n_loopcontact_df_genes = n_loopcontact_df_genes['gene'].values
                    n_loopcontacts[n][col] = n_loopcontact_df_genes

                    atleast_n_loopcontact_df[col] = atleast_n_loopcontact_df[col] >= n
                    atleast_n_loopcontact_df_genes = atleast_n_loopcontact_df[atleast_n_loopcontact_df[col] == True]
                    atleast_n_loopcontact_df_genes = atleast_n_loopcontact_df_genes['gene'].values
                    atleast_n_loopcontacts[n][col] = atleast_n_loopcontact_df_genes

        with open(n_loopcontacts_gene_list_outfile, 'wb') as fh:
            pickle.dump(n_loopcontacts, fh)
        print(f'SAVED: {n_loopcontacts_gene_list_outfile}')

        with open(atleast_n_loopcontacts_gene_list_outfile, 'wb') as fh:
            pickle.dump(atleast_n_loopcontacts, fh)
        print(f'SAVED: {atleast_n_loopcontacts_gene_list_outfile}')


        ############################################################################
        ## do a binomial regression for each set of genes to get OR, ci, and pvalues
        keys = [f'cut_{self.buff}_Rall', 'AA', 'region']
        n_loopcontacts_reg = []
        n_loopcontacts_outfile = os.path.join(self.OR_trendsOutpath, f'n_loopcontacts_reg_{self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}.csv')
        atleast_n_loopcontacts_reg = []
        atleast_n_loopcontacts_outfile = os.path.join(self.OR_trendsOutpath, f'atleast_n_loopcontacts_reg_{self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}.csv')

        if not os.path.exists(n_loopcontacts_outfile):
            for n in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:

                for k,v in n_loopcontacts[n].items():
                    print(f'{n} loop closing contacts of {k}: {len(v)}')
                    if len(v) != 0:
                        reg_df = self.resFeat_df[self.resFeat_df['gene'].isin(v)].copy()[keys]
                        reg_df[f'cut_{self.buff}_Rall'] = reg_df[f'cut_{self.buff}_Rall'].astype(int)
                        #print(reg_df)
                        reg_res = self.binom_reg(reg_df, f"cut_{self.buff}_Rall ~ AA + region")
                        reg_res['class'] = k
                        reg_res['nContacts'] = n
                        reg_res['nContacts_cond'] = 'exact'
                        reg_res['n_prot'] = len(v)
                        reg_res = reg_res[reg_res['var'] == 'region']
                        #print(reg_res)
                        n_loopcontacts_reg += [reg_res]
                
                for k,v in atleast_n_loopcontacts[n].items():
                    print(f'{n} loop closing contacts of {k}: {len(v)}')
                    if len(v) != 0:
                        reg_df = self.resFeat_df[self.resFeat_df['gene'].isin(v)].copy()[keys]
                        reg_df[f'cut_{self.buff}_Rall'] = reg_df[f'cut_{self.buff}_Rall'].astype(int)
                        #print(reg_df)
                        reg_res = self.binom_reg(reg_df, f"cut_{self.buff}_Rall ~ AA + region")
                        reg_res['class'] = k
                        reg_res['nContacts'] = n
                        reg_res['nContacts_cond'] = 'atleast'
                        reg_res['n_prot'] = len(v)
                        reg_res = reg_res[reg_res['var'] == 'region']
                        #print(reg_res)
                        atleast_n_loopcontacts_reg += [reg_res]

            n_loopcontacts_reg = pd.concat(n_loopcontacts_reg)
            n_loopcontacts_reg['OR'] = np.exp(n_loopcontacts_reg['coef'].astype(float))
            n_loopcontacts_reg['OR_lb'] = np.exp(n_loopcontacts_reg['[0.025'].astype(float))
            n_loopcontacts_reg['OR_ub'] = np.exp(n_loopcontacts_reg['0.975]'].astype(float))
            n_loopcontacts_reg['buff'] = self.buff
            n_loopcontacts_reg['spa'] = self.spa
            n_loopcontacts_reg['LiPMScov'] = self.LiPMScov
            n_loopcontacts_reg.to_csv(n_loopcontacts_outfile, sep='|', index=False)
            print(f'SAVED: {n_loopcontacts_outfile}')

            atleast_n_loopcontacts_reg = pd.concat(atleast_n_loopcontacts_reg)
            atleast_n_loopcontacts_reg['OR'] = np.exp(atleast_n_loopcontacts_reg['coef'].astype(float))
            atleast_n_loopcontacts_reg['OR_lb'] = np.exp(atleast_n_loopcontacts_reg['[0.025'].astype(float))
            atleast_n_loopcontacts_reg['OR_ub'] = np.exp(atleast_n_loopcontacts_reg['0.975]'].astype(float))
            atleast_n_loopcontacts_reg['buff'] = self.buff
            atleast_n_loopcontacts_reg['spa'] = self.spa
            atleast_n_loopcontacts_reg['LiPMScov'] = self.LiPMScov
            atleast_n_loopcontacts_reg.to_csv(atleast_n_loopcontacts_outfile, sep='|', index=False)
            print(f'SAVED: {atleast_n_loopcontacts_outfile}')

        else:
            n_loopcontacts_reg = pd.read_csv(n_loopcontacts_outfile, sep='|')
            print(f'LOADED: {n_loopcontacts_outfile}')
            atleast_n_loopcontacts_reg = pd.read_csv(atleast_n_loopcontacts_outfile, sep='|')
            print(f'LOADED: {atleast_n_loopcontacts_outfile}')

        print(f'n_loopcontacts_reg:\n{n_loopcontacts_reg}')
        print(f'atleast_n_loopcontacts_reg:\n{atleast_n_loopcontacts_reg}')
       

        ######################################################################################################################
        ### using the set of genes with 0 loop forming contacts as a ref get the set of differences for both n and atleast n OR 
        ref = n_loopcontacts_reg[n_loopcontacts_reg['nContacts'] == 0]
        print(ref)
        outdf = {'LoopContactClass': [], 'nContacts_cond':[], 'slope':[], 'p_value':[]}
        for LoopContactClass, LoopContactClass_df in ref.groupby('class'):
            ref_OR = LoopContactClass_df['OR'].values[0]
            ref_OR_ci = (LoopContactClass_df['OR_lb'].values[0], LoopContactClass_df['OR_ub'].values[0])
            print(f'{"#"*100}\nREF: {LoopContactClass} OR {ref_OR} {ref_OR_ci}')

            # for the genes with exactly n loop forming contacts of the class get the diff and fit a line to get the pvalue
            n_loopcontacts_reg_class = n_loopcontacts_reg[n_loopcontacts_reg['class'] == LoopContactClass] 
            #n_loopcontacts_reg_class = n_loopcontacts_reg_class[n_loopcontacts_reg_class['nContacts'].isin([0, 1, 2, 3, 4, 5])]
            n_loopcontacts_reg_class = n_loopcontacts_reg_class[n_loopcontacts_reg_class['nContacts'].isin([0, 1, 2, 3])]
            print(n_loopcontacts_reg_class)
            x = []
            y = []
            yerr = []
            n_prot = []
            for row_i, row in enumerate(n_loopcontacts_reg_class[['OR', 'OR_lb', 'OR_ub']].values):
                diff, err, ci = self.difference_with_error_propagation(ref_OR, ref_OR_ci, row[0], (row[1], row[2]))
                if err != np.inf:
                    x += [n_loopcontacts_reg_class['nContacts'].values[row_i]]
                    y += [(-1)*diff]
                    yerr += [err]
                    n_prot += [n_loopcontacts_reg_class['n_prot'].values[row_i]]
                    print(diff, err, ci)
            #x = n_loopcontacts_reg_class['nContacts'].values
            #n_prot = n_loopcontacts_reg_class['n_prot'].values
            outfile = os.path.join(self.OR_trendsOutpath, f'n_loopcontacts_ORdiff_linear_reg_{LoopContactClass}_{self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}.png')
            fit_res, slope, p_value = self.linear_fit_and_plot(x, y, yerr, n_prot, f'Exactly n {LoopContactClass} LFC | {self.tag}', outfile, ylabel=f'OR(n) - OR(n=0)')
            #print(fit_res)
            
            outdf['LoopContactClass'] += [LoopContactClass]
            outdf['nContacts_cond'] += ['exact']
            outdf['slope'] += [slope]
            outdf['p_value'] += [p_value]


            # for the genes with atleast n loop forming contacts of the class get the diff and fit a line to get the pvalue
            atleast_n_loopcontacts_reg_class = atleast_n_loopcontacts_reg[atleast_n_loopcontacts_reg['class'] == LoopContactClass] 
            #atleast_n_loopcontacts_reg_class = atleast_n_loopcontacts_reg_class[atleast_n_loopcontacts_reg_class['nContacts'].isin([0, 1, 2, 3, 4, 5])]
            atleast_n_loopcontacts_reg_class = atleast_n_loopcontacts_reg_class[atleast_n_loopcontacts_reg_class['nContacts'].isin([0, 1, 2, 3])]
            print(atleast_n_loopcontacts_reg_class)
            y = []
            yerr = []
            x = []
            n_prot = []
            for row_i, row in enumerate(atleast_n_loopcontacts_reg_class[['OR', 'OR_lb', 'OR_ub']].values):
                diff, err, ci = self.difference_with_error_propagation(ref_OR, ref_OR_ci, row[0], (row[1], row[2]))
                if err != np.inf:
                    x += [atleast_n_loopcontacts_reg_class['nContacts'].values[row_i]]
                    y += [(-1)*diff]
                    yerr += [err]
                    n_prot += [atleast_n_loopcontacts_reg_class['n_prot'].values[row_i]]
                    print(diff, err, ci)

            #x = atleast_n_loopcontacts_reg_class['nContacts'].values
            #n_prot = atleast_n_loopcontacts_reg_class['n_prot'].values
            outfile = os.path.join(self.OR_trendsOutpath, f'atleast_n_loopcontacts_ORdiff_linear_reg_{LoopContactClass}_{self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}.png')
            fit_res, slope, p_value = self.linear_fit_and_plot(x, y, yerr, n_prot, f'atleast n {LoopContactClass} LFC | {self.tag}', outfile, ylabel='OR(atleast-n) - OR(n=0)')
            #print(fit_res)
            
            outdf['LoopContactClass'] += [LoopContactClass]
            outdf['nContacts_cond'] += ['atleast']
            outdf['slope'] += [slope]
            outdf['p_value'] += [p_value]
            

        outdf = pd.DataFrame(outdf)
        outdf = outdf.sort_values(by=['p_value'])
        print(outdf)
        outfile = os.path.join(self.OR_trendsOutpath, f'ORdiff_linear_reg_{self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}.csv')
        outdf.to_csv(outfile, sep='|', index=False)
        print(f'SAVED: {outfile}')
    ################################################################################################


    ######################################################################
    def linear_fit_and_plot(self, x, y, yerr, n, title, outfile, ylabel='Y'):
        # Define the linear model: y = a * x + b
        def linear_model(x, a, b):
            return a * x + b

        # calculate the weights for y which are the inverse of the variance (1/sigma^2)
        print(f'y: {y}')
        print(f'x: {x}')
        yerr = np.asarray(yerr, dtype=float)
        print(f'yerr: {yerr}')
        weights  = 1/yerr**2
        print(f'yerr: {yerr}')
        print(f'weights: {weights}')


        # Perform the weighted fit using curve_fit
        popt, pcov = curve_fit(linear_model, x, y, sigma=yerr, absolute_sigma=True)

        # Extract the parameters and their standard errors
        slope, intercept = popt
        a_err, b_err = np.sqrt(np.diag(pcov))

        # Calculate the t-statistic for the slope
        t_stat = slope / a_err

        # Degrees of freedom
        df = len(x) - 2

        # Calculate the p-value
        p_value = 2 * (1 - t.cdf(np.abs(t_stat), df))
        print(f"Slope: {slope:.4f} ±{a_err:.4f}")
        print(f"Intercept: {intercept:.4f} ± {b_err:.4f}")
        print(f"P-value for slope different from 0: {p_value:.4f}")
   
        # Generate points for the fitted line
        x_fit = np.linspace(min(x), max(x), 100)
        y_fit = linear_model(x_fit, *popt)

        ## using the stats model package if necesasry
        # Add a constant to the independent variable array (for the intercept)
        #x_with_constant = sm.add_constant(x)        
        # Fit the model
        #model = sm.OLS(y, x_with_constant, weights=weights).fit()
        # Get the model summary
        #summary = model.summary()
        
        ## save plot data
        plot_df = {'n_contacts':x, 'OddsDiff':y, 'OddsDiff_lb':yerr[0], 'OddsDiff_ub':yerr[1], 'sample_size':n}
        plot_df = pd.DataFrame(plot_df)
        print(f'plot_df:\n{plot_df}')
        outfile_csv = outfile.replace('.png', '.csv')
        plot_df.to_csv(outfile_csv)
        logging.info(f'SAVED: {outfile_csv}')

        # Plot the data and the regression line
        fig,ax = plt.subplots(2,1, figsize=(6,4), sharex=True)

        ax[0].errorbar(x, y, yerr=yerr, label='Data Points', color='k', capsize=3, marker='o', ls='None')
        ax[0].plot(x_fit, y_fit, color='red', label=f'Linear fit: y = {slope:.2f}x + {intercept:.2f}')
        #ax[0].plot(x, model.predict(x_with_constant), color='red', label='Regression Line')
        ax[0].set_ylabel(ylabel)

        # Annotate the plot with the regression equation and p-value
        #slope = model.params[1]
        #intercept = model.params[0]
        #p_value = model.pvalues[1]
        equation_text = f'y = {intercept:.2f} + {slope:.2f}x\np-value = {p_value:.3e}'
        
        ax[0].annotate(equation_text, xy=(1.05, 0.95), xycoords='axes fraction', fontsize=12,
                     ha='left', va='top', bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="white"))
        
        ax[1].plot(x, n, color='k', marker='o', ls='None')
        ax[1].set_ylabel('sample size')
        ax[1].set_yscale('log')
        ax[1].set_xlabel('n')
        ax[1].set_ylim(1,500)
        plt.suptitle(title)
        plt.tight_layout()
        #plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        #plt.grid(True)
        
        # Show the plot
        #plt.show()
        plt.savefig(outfile)
        plt.close()
        logging.info(f'SAVED: {outfile}')
        
        #return summary, slope, p_value
        return intercept, slope, p_value
    ######################################################################


    ######################################################################
    def difference_with_error_propagation(self, mean1, ci1, mean2, ci2):
        """
        Calculate the difference between two means and propagate the error.
        
        Parameters:
        mean1 (float): The first mean value.
        ci1 (tuple): Confidence interval for the first mean (lower, upper).
        mean2 (float): The second mean value.
        ci2 (tuple): Confidence interval for the second mean (lower, upper).
        
        Returns:
        diff (float): The difference between the two means.
        error (float): The propagated error for the difference.
        ci_diff (tuple): The confidence interval for the difference.
        """
        
        # Calculate the errors from the confidence intervals
        error1 = (ci1[1] - ci1[0]) / 2
        error2 = (ci2[1] - ci2[0]) / 2
        
        # Calculate the difference between the means
        diff = mean1 - mean2
        
        # Propagate the error
        propagated_error = np.sqrt(error1**2 + error2**2)
        
        # Calculate the confidence interval for the difference
        ci_diff = (diff - propagated_error, diff + propagated_error)
        
        return diff, propagated_error, ci_diff
    ######################################################################

    ######################################################################
    def binom_reg(self, df, formula):
        model = sm.GLM.from_formula(formula, family=sm.families.Binomial(), data=df)
        result = model.fit()

        ## recalculate the pvalue to add more digits as statsmodels truncates it to 0 if it is below 0.0001 for some reason. 
        table = result.summary().tables[1]
        table_df = pd.DataFrame(table.data[1:], columns=table.data[0])
        pvalues = []
        for z in table_df['z']:
            z = float(z)
            if z < 0:
                p = st.norm.cdf(z)
            else:
                p = 1 - st.norm.cdf(z)
            pvalues += [p*2]

        table_df['P>|z|'] = pvalues
        table_df = table_df.rename(columns={"": "var"})
        return table_df
    ######################################################################


    ######################################################################
    def loglinear(self,):

        """
        Test 2: Is the OR (misfolding Y/N, entangled region Y/N) dependant on the number of SLFC at a given SPA threshold
                i. Make the following groups of proteins
                        (1) has atleast n SLFC (n = 0, 1, 2, 3, 4, 5)
                        (2) has n SLFC (n = 0, 1, 2, 3, 4, 5)
                ii. For each group calculate the OR using the contingency table and estimate bootstrapped 95% ci
                iii. calculate the delta between n=0 and the other n values for each group and determine if the slope is signifcantly different from zero
        """

        ## make directory to output test 1 info to 
        self.LogLinearOutpath = os.path.join(self.Outpath, 'LogLinear/')
        if not os.path.exists(self.LogLinearOutpath):
            os.makedirs(self.LogLinearOutpath)
            print(f'Made directory: {self.LogLinearOutpath}')

        
        ############################################################################
        ## for each of the loop contacting types in the self.loopcontact_df object get the list of genes with n or => n loop forming contacts
        # List of columns to exclude from the operation
        exclude_columns = ['gene', 'Total', 'Essential']

        ############################################################################
        ## do a binomial regression for each set of genes to get OR, ci, and pvalues
        n_loopcontact_df = self.loopcontact_df[self.loopcontact_df['gene'].isin(self.gene_list)].copy()
        loglinear_dfs = []
        for col in n_loopcontact_df.columns:
            if col not in exclude_columns:
                print(f'{"#"*50}\nCreating loglinear df for LoopContactingCLass: {col}')

                keys = ['gene', f'cut_{self.buff}_Rall', 'region', 'AA']
                loglinear_df = self.resFeat_df[self.resFeat_df['gene'].isin(self.gene_list)].copy()[keys]
                dfs = []
                for gene, gene_df in loglinear_df.groupby('gene'):
                    n = n_loopcontact_df[n_loopcontact_df['gene'] == gene]
                    n = n[col].values[0]
                    gene_df['n'] = n
                    dfs += [gene_df]

                loglinear_df = pd.concat(dfs, ignore_index=True)
                loglinear_df[f'cut_{self.buff}_Rall'] = loglinear_df[f'cut_{self.buff}_Rall'].astype(int)
                reg_df = loglinear_df[[f'cut_{self.buff}_Rall', 'region', 'n', 'AA']]
                loglinear_df = loglinear_df[[f'cut_{self.buff}_Rall', 'region', 'n']]

                ## # Create a contingency table from the data
                contingency_table = pd.crosstab(index=[loglinear_df['region'], loglinear_df['n']], columns=loglinear_df[f'cut_{self.buff}_Rall'])

                # Flatten the contingency table
                contingency_table = contingency_table.stack().reset_index(name='count')
                print(contingency_table)

                # Fit the log-linear model
                # The formula now includes interactions between A, B (with 4 levels), and C
                #model = sm.GLM.from_formula(f'count ~ region + n + cut_{self.buff}_Rall', data=contingency_table, family=sm.families.Poisson()).fit()
                model = sm.GLM.from_formula(f'count ~ region * n * cut_{self.buff}_Rall', data=contingency_table, family=sm.families.Poisson()).fit()

                # Display the summary of the model
                table = model.summary().tables[1]
                print(model.summary())
                table_df = pd.DataFrame(table.data[1:], columns=table.data[0])
                print(table_df)
                quit()
                pvalues = []
                for z in table_df['z']:
                    z = float(z)
                    if z < 0:
                        p = st.norm.cdf(z)
                    else:
                        p = 1 - st.norm.cdf(z)
                    pvalues += [p*2]

                table_df['P>|z|'] = pvalues
                table_df = table_df.rename(columns={"": "var"})
                #print(model.summary())
                table_df['buff'] = self.buff
                table_df['spa'] = self.spa
                table_df['LiPMScov'] = self.LiPMScov
                table_df['class'] = col
                print(table_df)

                loglinear_dfs += [table_df]
        
        loglinear_dfs = pd.concat(loglinear_dfs, ignore_index=True)
        print(f'loglinear_dfs:\n{loglinear_dfs}')
        loglinear_dfs = loglinear_dfs.rename(columns={"P>|z|": "pvalue"})
        print(f'loglinear_dfs:\n{loglinear_dfs}')
        title = 'log(count) ~ B0 + B1*n + B2*region + B3*cut +...\n...+ B12*n*region + B23*region*cut + B123*n*region*cut'
        #title = 'log(count) ~ B0 + B1*n + B2*region + B3*cut + e'
        for var, var_df in loglinear_dfs.groupby('var'):
            self.plot_OLSorLogLinear(var_df, os.path.join(self.LogLinearOutpath, f'LogLinear_reg_var-{var}_{self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}.png'), title, ylabel=var)
        outfile = os.path.join(self.LogLinearOutpath, f'LogLinear_reg_{self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}.csv')
        loglinear_dfs.to_csv(outfile, sep='|', index=False)
        print(f'SAVED: {outfile}')
    ######################################################################

    ######################################################################
    def OLSfit(self, Contact_cond='atleast'):

        """
        Test 4: Model the ratio of cut proportions as a function of the number of loop closing contacts
        (Ncut-entR/NentR)/(Ncut-nonentR/NnonentR) ~ B0 + B1*n + e
        where the group of proteins used to calcualte the OR can have either exactly the n number of loop closing contacts specified by the user or atleast depending on the Contact_cond
        """

        ## make directory to output test 1 info to 
        self.OLSfitOutpath = os.path.join(self.Outpath, 'OLSfit/')
        if not os.path.exists(self.OLSfitOutpath):
            os.makedirs(self.OLSfitOutpath)
            print(f'Made directory: {self.OLSfitOutpath}')

        
        ############################################################################
        ## for each of the loop contacting types in the self.loopcontact_df object get the list of genes with n or => n loop forming contacts
        # List of columns to exclude from the operation
        exclude_columns = ['gene', 'Total', 'Essential']

        ############################################################################
        ## do a binomial regression for each set of genes to get OR, ci, and pvalues
        n_loopcontact_df = self.loopcontact_df[self.loopcontact_df['gene'].isin(self.gene_list)].copy()
        OLS_dfs = []
        for col in n_loopcontact_df.columns:
            if col not in exclude_columns:
                print(f'{"#"*50}\nCreating OLS df for LoopContactingCLass: {col}')

                keys = ['gene', f'cut_{self.buff}_Rall', 'region', 'AA']
                OLS_df = self.resFeat_df[self.resFeat_df['gene'].isin(self.gene_list)].copy()[keys]
                reg_df = {'OR':[], 'n':[]}
                for n in n_loopcontact_df[col].unique():
                    #print(col, n)

                    if Contact_cond == 'exact':
                        genes_with_n = n_loopcontact_df[n_loopcontact_df[col] == n]['gene'].values
                    elif Contact_cond == 'atleast':
                        genes_with_n = n_loopcontact_df[n_loopcontact_df[col] >= n]['gene'].values
                    #print(genes_with_n)

                    genes_with_n_resfeat = OLS_df[OLS_df['gene'].isin(genes_with_n)][keys]
                    #print(genes_with_n_resfeat)

                    table = pd.crosstab(genes_with_n_resfeat[f'cut_{self.buff}_Rall'], genes_with_n_resfeat['region'])
                    table += 1
                    if table.values.shape == (2,2):
                        #print(table)
                        cutY_regionY = table.values[1,1]
                        cutY_regionN = table.values[1,0]
                        cutN_regionY = table.values[0,1]
                        cutN_regionN = table.values[0,0]

                        #diff = (cutY_regionY/(cutY_regionY + cutN_regionY)) - (cutY_regionN/(cutY_regionN + cutN_regionN))
                        #reg_df['OR'] += [diff]
                        #reg_df['n'] += [n]
                        #print(col, n, diff)
                        OR, pvalue = fisher_exact(table)
                        if OR != np.inf:
                            #print(col, n, OR, pvalue, len(genes_with_n))
                            reg_df['OR'] += [OR]
                            reg_df['n'] += [n]

                reg_df = pd.DataFrame(reg_df)
                reg_df = reg_df.sort_values(by=['n'])
                #print(reg_df)

                # Fit the model
                model = smf.ols('OR ~ n', data=reg_df).fit()

                # Print the summary of the model
                #print(model.summary())
                # Display the summary of the model
                table = model.summary().tables[1]
                table_df = pd.DataFrame(table.data[1:], columns=table.data[0])
                pvalues = []
                for t in table_df['t']:
                    t = float(t)
                    if t < 0:
                        p = st.norm.cdf(t)
                    else:
                        p = 1 - st.norm.cdf(t)
                    pvalues += [p*2]

                #table_df['P>|t|'] = pvalues
                table_df = table_df.rename(columns={"": "var"})
                #print(model.summary())
                table_df['buff'] = self.buff
                table_df['spa'] = self.spa
                table_df['LiPMScov'] = self.LiPMScov
                table_df['class'] = col
                print(table_df)

                OLS_dfs += [table_df]

        OLS_dfs = pd.concat(OLS_dfs, ignore_index=True)
        OLS_dfs = OLS_dfs[OLS_dfs['var'] == 'n']
        #OLS_dfs = OLS_dfs.sort_values('P>|t|')
        OLS_dfs = OLS_dfs.rename(columns={"P>|t|": "pvalue"})
        print(f'OLS_dfs:\n{OLS_dfs}')
        title = 'OR ~ B0 + B1*n + e'
        self.plot_OLSorLogLinear(OLS_dfs, os.path.join(self.OLSfitOutpath, f'OLSfit_reg_{Contact_cond}-n_{self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}.png'), title, ylabel='n')
        outfile = os.path.join(self.OLSfitOutpath, f'OLSfit_reg_{Contact_cond}-n_{self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}.csv')
        OLS_dfs.to_csv(outfile, sep='|', index=False)
        print(f'SAVED: {outfile}')
    #################################################################################################################

    #################################################################################################################
    def plot_OLSorLogLinear(self,df, outfile, title, ylabel):

        print(df)
        # Plot the data and the regression line
        fig,ax = plt.subplots(2,1, figsize=(5,3.5), sharex=True)
        
        x = df['class']
        y1 = []
        y1err = []
        for row in df[['coef', '[0.025', '0.975]']].values.astype(float):
            print(row)
            y1 += [row[0]]
            y1err += [[row[0] - row[1], row[2] - row[0]]]
        print(y1)
        print(y1err)
        y1err = np.asarray(y1err).T

        pvalues = df['pvalue'].values.astype(float)
        pvalues = false_discovery_control(pvalues) 
        pvalues = np.where(pvalues == 0, 0.000001, pvalues)

        ax[0].errorbar(x, y1, yerr=y1err, label='Essential', capsize=3, marker='o', ls='None', markerfacecolor='none')
        ax[0].set_ylabel(ylabel, fontsize=7)
        ax[0].tick_params(axis='both', labelsize=7)
        ax[0].axhline(y=0, color='red', linestyle='--', linewidth=2)

        ax[1].scatter(x, pvalues)
        ax[1].set_yscale('log')
        ax[1].axhline(y=0.05, color='red', linestyle='--', linewidth=2)
        ax[1].set_ylabel('pvalue (FDR)', fontsize=7)
        ax[1].tick_params(axis='both', labelsize=7)
        
        # R3tate the x-axis tick labels by 45 degrees for all subplots
        plt.setp(ax[1].get_xticklabels(), rotation=45, ha="right", fontsize=7)            

        plt.suptitle(title, fontsize=8)
        plt.tight_layout()
        #plt.grid(True)
        
        # Show the plot
        #plt.show()
        
        plt.savefig(outfile)
        plt.close()
    #################################################################################################################


    #################################################################################################################
    def Binom_reg_n_dep(self,):
        """
        Do a binomial logistic regression on the odds of cutting but include the atleast n variable
        cut ~ AA + region + n
        """

        ## make directory to output test 1 info to 
        self.Binom_reg_n_depOutpath = os.path.join(self.Outpath, 'Binom_reg_n_dep/')
        if not os.path.exists(self.Binom_reg_n_depOutpath):
            os.makedirs(self.Binom_reg_n_depOutpath)
            print(f'Made directory: {self.Binom_reg_n_depOutpath}')

        
        ############################################################################
        ## for each of the loop contacting types in the self.loopcontact_df object get the list of genes with n or => n loop forming contacts
        # List of columns to exclude from the operation
        exclude_columns = ['gene', 'Total', 'Essential']

        ############################################################################
        ## do a binomial regression for each set of genes to get OR, ci, and pvalues
        n_loopcontact_df = self.loopcontact_df[self.loopcontact_df['gene'].isin(self.gene_list)].copy()
        BinomReg_dfs = []
        for col in n_loopcontact_df.columns:
            if col not in exclude_columns:
                print(f'{"#"*100}\nCreating Binomial df for LoopContactingCLass: {col}')

                keys = ['gene', f'cut_{self.buff}_Rall', 'region', 'AA']
                BinomReg_df = self.resFeat_df[self.resFeat_df['gene'].isin(self.gene_list)].copy()[keys]
                dfs = []
                for gene, gene_df in BinomReg_df.groupby('gene'):
                    n = n_loopcontact_df[n_loopcontact_df['gene'] == gene]
                    n = n[col].values[0]
                    gene_df['n'] = n
                    dfs += [gene_df]

                BinomReg_df = pd.concat(dfs, ignore_index=True)
                BinomReg_df[f'cut_{self.buff}_Rall'] = BinomReg_df[f'cut_{self.buff}_Rall'].astype(int)
                BinomReg_df = BinomReg_df[[f'cut_{self.buff}_Rall', 'region', 'n', 'AA']]
                BinomReg_df = BinomReg_df.sort_values(by=['n'])
                print(f'BinomReg_df:\n{BinomReg_df}')

                reg_res = self.binom_reg(BinomReg_df, f"cut_{self.buff}_Rall ~ AA + region + n")
                #reg_res = self.binom_reg(BinomReg_df, f"cut_{self.buff}_Rall ~ AA * region * n")
                print(reg_res)

                #print(model.summary())
                reg_res['buff'] = self.buff
                reg_res['spa'] = self.spa
                reg_res['LiPMScov'] = self.LiPMScov
                reg_res['class'] = col
                print(reg_res)

                BinomReg_dfs += [reg_res]
        
        BinomReg_dfs = pd.concat(BinomReg_dfs, ignore_index=True)
        BinomReg_dfs = BinomReg_dfs[BinomReg_dfs['var'].isin(['region', 'n'])]
        BinomReg_dfs = BinomReg_dfs.sort_values(by=['var', 'P>|z|'])
        BinomReg_dfs = BinomReg_dfs.rename(columns={"P>|z|": "pvalue"})
        print(f'BinomReg_dfs:\n{BinomReg_dfs}')
        title = 'log(odds) ~ B0 + B1*n + B2*region + B3*cut + B12*n*region + B23*region*cut + B123*n*region*cut'
        for var, var_df in BinomReg_dfs.groupby('var'):
            self.plot_OLSorLogLinear(var_df, os.path.join(self.Binom_reg_n_depOutpath, f'Binom_reg_n_dep_reg_var-{var}_{self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}.png'), title, ylabel=var)
        outfile = os.path.join(self.Binom_reg_n_depOutpath, f'Binom_reg_n_dep_reg_{self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}.csv')
        BinomReg_dfs.to_csv(outfile, sep='|', index=False)
        print(f'SAVED: {outfile}')


#################################################################################################################

def statistic(x, y, axis):
    return np.mean(x, axis=axis) - np.mean(y, axis=axis)

def main():
    """
    Main function to control workflow.
    """

    # Parse the user supplied arguments
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-g", "--gene_list", type=str, required=True, help=f"path to all ent gene list use")
    parser.add_argument("-Eg", "--ess_gene_list", type=str, required=True, help=f"path to essentail ent gene list use")
    parser.add_argument("-NEg", "--noness_gene_list", type=str, required=True, help=f"path to nonessentail ent gene list use")
    parser.add_argument("-r", "--resFeat_files", type=str, required=True, help="path to residue Feature files")
    parser.add_argument("-c", "--loopcontact_df", type=str, required=True, help="path to dataframe containing loop contact classifications")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="path to output directory. will be made if doesnt exist")
    parser.add_argument("-t", "--tag", type=str, required=True, help="tag for final output image")
    parser.add_argument("-b", "--buff", type=str, required=True, help="buffer used C CD CG")
    parser.add_argument("-s", "--spa", type=str, required=True, help="spa used 0, 10, 20, 30, 40, 50, 60, 70, 80, 90")
    parser.add_argument("--LiPMScov", type=str, required=True, help="LiPMS coverage used 0, 10, 20, 30, 40, 50, 60, 70, 80, 90")
    parser.add_argument("-l", "--log_file", type=str, required=True, help="Path to logging file")
    args = parser.parse_args()

    resFeat_files = args.resFeat_files
    log_file = args.log_file
    gene_list = args.gene_list
    ess_gene_list = args.ess_gene_list
    noness_gene_list = args.noness_gene_list
    outpath = args.outpath
    tag = args.tag
    buff = args.buff
    spa = args.spa
    LiPMScov = args.LiPMScov
    loopcontact_df = args.loopcontact_df

    # Setup logging configuration
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s %(message)s')
    logging.info(f'{"#"*100}\nNEW RUN')

    # Make outpath if it doesnt exists
    if not os.path.exists(outpath):
        os.makedirs(outpath)
        print(f'Made directory: {outpath}')

    # Initalize the FrequencyGenerator class object
    analysis = Analyzer(outpath, tag, buff, resFeat_files, spa, LiPMScov, gene_list, ess_gene_list, noness_gene_list, loopcontact_df)
    print(analysis)

    # Do Analysis 1: fraction of genes with atleast one of the significant loop forming contacts
    analysis.FractionGenes_w_n_LoopContacts()

    # Do analysis 2: comparing distributions of loop closing contacts for each class
    #analysis.LoopContactsDistsCompare()

    # Do analysis 3: trends in OR 
    analysis.OR_trends()
    
    # Do analysis 3: log linear model
    #analysis.loglinear()
    

    # Do analysis 4: Ordinary least squares fit to model ratio of cut proportions in each bucket of n
    #analysis.OLSfit(Contact_cond='exact')
    #analysis.OLSfit(Contact_cond='atleast')

    # Do analysis 5: Binomial regression with dependance on n -> cut ~ AA + region + n
    #analysis.Binom_reg_n_dep()

    print(f'Logfile: {log_file}')


if __name__ == "__main__":
    main()

print('NORMAL TERMINATION')
logging.info('NORMAL TERMINATION')
