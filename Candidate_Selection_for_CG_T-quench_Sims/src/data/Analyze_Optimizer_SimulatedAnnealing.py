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
    def __init__(self, inpdir, outpath, tag, resFeat_files, reg_formula):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:

        """
        self.inpdir = inpdir
        self.outpath = outpath
        self.tag = tag
        self.resFeat_files = resFeat_files
        self.reg_formula = reg_formula

        if not os.path.exists(f'{self.outpath}'):
            os.makedirs(f'{self.outpath}')
            print(f'Made output directories {self.outpath}')
    ##########################################################################################################

    ##########################################################################################################

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
                final_traj[traj][state] = df

            if 'Final_step_reg' in f:
                if traj not in Final_step_reg:
                    Final_step_reg[traj] = {}
                Final_step_reg[traj] = pd.read_csv(f)


        #load res features files
        res_files = glob.glob(os.path.join(self.resFeat_files, '*'))
        print(f'Number of res_files: {len(res_files)}')
        res_feat_df = pd.concat([pd.read_csv(f, sep='|') for f in res_files])
        print(f'res_feat_df:\n{res_feat_df}')
                # Quality check to ensure all columns are present in the df
        reg_vars = [v for v in self.reg_formula.split(' ') if v not in ['~', '+']]
        keys = []
        cut_key = ''
        for reg_var in reg_vars:
            if '*' in reg_var:
                reg_var = reg_var.split('*')
            else:
                reg_var = [reg_var]
            for v in reg_var:

                if 'cut_' in v:
                    cut_key = v

                if v not in res_feat_df:
                    logging.info(f"Regression variable {reg_var} not in the DataFrame")
                    sys.exit()
                else:
                    keys += [v]

        logging.info(f'cut_key: {cut_key}')
        print(f'cut_key: {cut_key}')
        self.cut_key = cut_key
        keys = list(set(keys))
        bin_keys = keys.copy()
        if 'AA' in keys:
            bin_keys.remove('AA')
        logging.info(f'bin_keys: {bin_keys}')

        # Encode boolean columns
        encoded_df = self.encode_boolean_columns(res_feat_df, bin_keys)
        encoded_df = encoded_df[['gene']+reg_vars]
        print(f'encoded_df:\n{encoded_df}')

        # Get the per gene number of residues
        prot_size = {'gene':[], 'prot_size':[]}
        for gene, gene_df in encoded_df.groupby('gene'):
            prot_size['gene'] += [gene]
            prot_size['prot_size'] += [len(gene_df)]
        prot_size = pd.DataFrame(prot_size)
        print(f'prot_size:\n{prot_size}')
        self.prot_size = prot_size

        return final_genelist, final_traj, Final_step_reg, encoded_df
    ##########################################################################################################

    ##########################################################################################################
    def encode_boolean_columns(self, df: pd.DataFrame, boolean_columns: list) -> pd.DataFrame:
        """
        Encodes boolean-like columns in a DataFrame to binary 0 and 1.

        Parameters:
        - df (pd.DataFrame): The input DataFrame.
        - boolean_columns (list): A list of column names to be encoded.

        Returns:
        - pd.DataFrame: The DataFrame with encoded columns.
        """
        label_encoder = LabelEncoder()

        for column in boolean_columns:
            if column in df.columns:
                df[column] = label_encoder.fit_transform(df[column])
            else:
                print(f"Column '{column}' does not exist in the DataFrame.")

        return df
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
                #print(state, state_df)
                final_OR = float(state_df['OR'].values[-1])
                if final_OR > max_OR:
                    max_OR = final_OR
                    max_state = state
            
            print(f'traj {traj} with max_state: {max_state} with OR {max_OR}')
            max_state_genes[traj] = final_genelist_data[traj][max_state]
            #print(max_state_genes[traj])
            total_genes += [max_state_genes[traj]]
      
        total_traj = len(traj_data)
        print(f'Total number of traj: {total_traj}')
        total_genes = np.hstack(total_genes)
        total_genes = set(total_genes)
        total_genes = list(total_genes)
        print(f'total_genes: {total_genes} {len(total_genes)}')

        ## for each gene in total_genes calculate the fraction of traj it was present in the top state
        dist_df = {'gene':[], 'F':[]}
        for gene in total_genes:
            F = 0
            for traj, traj_genes in max_state_genes.items():
                if gene in traj_genes:
                    F += 1
            F /= total_traj
            #print(gene, F)
            dist_df['gene'] += [gene]
            dist_df['F'] += [F]
        dist_df = pd.DataFrame(dist_df)
        #print(f'dist: {dist}')

        ## Plot the data
        # Create a figure with two subplots in the same row
        fig, ax = plt.subplots(1, 2, figsize=(12, 5))
        plot_outfile = os.path.join(self.outpath, f'Fract_High-OR_dists_{title}.png')
        hist_outfile = os.path.join(self.outpath, f'Fract_High-OR_Hdists_{title}.csv')
        cdf_outfile = os.path.join(self.outpath, f'Fract_High-OR_Cdists_{title}.csv')

        # First subplot: Normalized Histogram
        H, bins, _ = ax[0].hist(dist_df['F'].values, bins=20, density=False, alpha=0.7, color='blue', width=0.01)
        ax[0].set_title('Normalized Histogram')
        ax[0].set_xlabel('Fraction traj in High-OR')
        ax[0].set_ylabel('Frequency')
        ax[0].set_xlim(0,1)
        bins = bins[:-1]
       
        outdf = pd.DataFrame({'Hist':H, 'BinEdges':bins})
        outdf.to_csv(hist_outfile, index=False)
        print(f'SAVED: {hist_outfile}')


        # Second subplot: CDF
        #sorted_data = np.sort(dist)
        sorted_data = dist_df.sort_values(by=['F'])
        sorted_data = sorted_data.reset_index(drop=True)
        #print(sorted_data)
        sorted_F = sorted_data['F'].values
        cdf = np.arange(1, len(sorted_F) + 1) / len(sorted_F)
        
        ax[1].bar(sorted_F, cdf, width=0.01, color='red', alpha=0.7)
        ax[1].set_title('Cumulative Distribution Function (CDF)')
        ax[1].set_xlabel('Fraction traj in High-OR')
        ax[1].set_ylabel('CDF')
        ax[1].set_xlim(0,1)

        sorted_data['CDF'] = cdf
        print(sorted_data)
        #outdf = pd.DataFrame({'CDF':cdf, 'sorted_F':sorted_data})
        sorted_data.to_csv(cdf_outfile, index=False)
        print(f'SAVED: {cdf_outfile}')

        # Display the plots
        plt.suptitle(title)
        plt.tight_layout()
        #plt.show()
        plt.savefig(plot_outfile)
        plt.close()
        print(f'SAVED: {plot_outfile}')
    ##########################################################################################################

    ##########################################################################################################
    def Plot_ensem(self, data, title, num_states):
        ## Data is a tuple with three elements: final_genelist, final_traj, Final_step_reg
        final_gene_lists = data[0]
        traj_data = data[1]
        final_reg_data = data[2]
        res_feat_df = data[3]

        ## determine if any of the states switch order within the last 10 frames
        for traj, traj_states in traj_data.items():
            last10_order = {}
            # loop over all the final frames of each state of each traj
            for iframe in range(-10, 0):
                traj_ranking = []
                for state, state_df in traj_states.items():
                    
                    #final_OR = state_df['OR'].values[-1]
                    final_OR = state_df['OR'].values[iframe]
                    traj_ranking += [[int(state), float(final_OR)]]

                traj_ranking = np.asarray(traj_ranking)
                traj_ranking = traj_ranking[traj_ranking[:,1].argsort()[::-1]]
                print(f'traj_ranking:\n{traj_ranking}')

                print('iframe, traj, s, state, OR')
                for s, (state, OR) in enumerate(traj_ranking):
                    if s not in last10_order:
                        last10_order[s] = [state]
                    else:
                        last10_order[s] += [state]
                    print(iframe, traj, s, state, OR)

            for ranked_s, states in last10_order.items():
                print(ranked_s, states)
                if len(np.unique(states)) != 1:
                    logging.info(f'NOTICE: one of the states here switched order in traj: {traj} for the ranked ordered state: {ranked_s}')

        logging.info(f'Final set of ranked states')
        ensemble_state_idxs = {s:[] for s in range(num_states)}
        ranking_summary_df = {'traj':[], 'ranked_stateID':[], 'old_stateID':[], 'OR':[]}
        ## get final set of ensemble states
        for traj, traj_states in traj_data.items():
                       
                traj_ranking = []
                for state, state_df in traj_states.items():
                    
                    #final_OR = state_df['OR'].values[-1]
                    final_OR = np.mean(state_df['OR'].values[-10:])
                    traj_ranking += [[int(state), float(final_OR)]]
                traj_ranking = np.asarray(traj_ranking)
                traj_ranking = traj_ranking[traj_ranking[:,1].argsort()[::-1]]

                logging.info('traj, ranked-s, state, OR')
                for s, (state, OR) in enumerate(traj_ranking):
                    print(f'{traj}, {s}, {state}, {OR}')
                    logging.info(f'{traj}, {s}, {state}, {OR}')
                    ranking_summary_df['traj'] += [traj]
                    ranking_summary_df['ranked_stateID'] += [s]
                    ranking_summary_df['old_stateID'] += [str(int(state))]
                    ranking_summary_df['OR'] += [OR]
                    df = traj_states[str(int(state))]
                    df['traj'] = traj
                    ensemble_state_idxs[s] += [df]
    
        ## save ths ranking summary dataframe so you can track back the original states if necessary
        ranking_summary_df = pd.DataFrame(ranking_summary_df)
        print(ranking_summary_df)
        ranking_summary_outfile = os.path.join(self.outpath, f'ranking_summary_{title}.csv')

        ## get the contingency tables for each ranked state and save it
        ranking_summary_df = self.get_ranked_state_ctables(ranking_summary_df, res_feat_df, final_gene_lists, ranking_summary_outfile)
        print(f'ranking_summary_df:\n{ranking_summary_df}')
 
        fig, axes = plt.subplots(nrows=4, ncols=num_states, figsize=(12, 10), sharex=True)
        outfile = os.path.join(self.outpath, f'Ensemble_Traj_{title}.png')
        data_outfile = os.path.join(self.outpath, f'Ensemble_Traj_{title}.csv')
        outdfs = []
        for s_idx, (s, s_dfs) in enumerate(ensemble_state_idxs.items()):
            
            ## get protein sizes for this state and calcualte the per reisdue number of cuts to plot since it is intensive
            for df in s_dfs:
                
                old_s = df['state'].values[-1]
                traj = df['traj'].values[-1]
                total = ranking_summary_df[(ranking_summary_df['ranked_stateID'] == s) & (ranking_summary_df['old_stateID'] == str(old_s)) & (ranking_summary_df['traj'] == str(traj))]
                total = total['total'].values[0]
                #print(f'new_state_ID: {s} | old_state_ID: {old_s} | traj: {traj} | total res: {total}')
                df['pcuts'] = df['cuts']/total

            # get the average across all dataframe in state
            means_df, ci_lower_df, ci_upper_df = running_average_and_ci(s_dfs, ['step', 'OR', 'pcuts', 'cuts', 'pvalue', 'psize_mean'])

            x = means_df['step'].values

            OR = means_df['OR'].values
            OR_lb = OR - ci_lower_df['OR']
            OR_ub = ci_upper_df['OR'] - OR

            pvals = means_df['pvalue'].values
            pvals_lb = pvals - ci_lower_df['pvalue']
            pvals_ub = ci_upper_df['pvalue'] - pvals

            pcuts = means_df['pcuts'].values
            pcuts_lb = pcuts - ci_lower_df['pcuts']
            pcuts_ub = ci_upper_df['pcuts'] - pcuts

            cuts = means_df['cuts'].values
            cuts_lb = cuts - ci_lower_df['cuts']
            cuts_ub = ci_upper_df['cuts'] - cuts

            size_means = means_df['psize_mean'].values
            size_means_lb = size_means - ci_lower_df['psize_mean']
            size_means_ub = ci_upper_df['psize_mean'] - size_means

            # make df for output
            outdf = {'MCsteps':x, 'OR':OR, 'OR_lb':OR_lb, 'OR_ub':OR_ub, 
                    'pvals':pvals, 'pvals_lb':pvals_lb, 'pvals_ub':pvals_ub,
                    'pcuts':pcuts, 'pcuts_lb':pcuts_lb, 'pcuts_ub':pcuts_ub,
                    'cuts':cuts, 'cuts_lb':cuts_lb, 'cuts_ub':cuts_ub,
                    'size_means':size_means, 'size_means_lb':size_means_lb, 'size_means_ub':size_means_ub, }
            outdf = pd.DataFrame(outdf)
            outdf['ranked_state'] = s
            outdfs += [outdf]
        
            # plot OR
            axes[0, s_idx].errorbar(x, OR, yerr=[OR_lb, OR_ub], label=f'state-{s}', marker='o', ls='none', fillstyle='none', capsize=3)
            axes[0, s_idx].plot(x, OR, label=f'state-{s}', ls='-', color='k', zorder=10)
            axes[0, s_idx].axhline(y=1, color='red', linestyle='--')
            axes[0, s_idx].set_ylabel('OR')
            #axes[0, s_idx].set_xlabel('MC step')
            axes[0, s_idx].set_title(f'State: {s}')


            # plot pvals
            axes[1, s_idx].errorbar(x, pvals, yerr=[pvals_lb, pvals_ub], label=f'state-{s}', marker='o', ls='none', fillstyle='none', capsize=3)
            axes[1, s_idx].plot(x, pvals, label=f'state-{s}', ls='-', color='k', zorder=10)
            axes[1, s_idx].axhline(y=0.05, color='red', linestyle='--', zorder=10)
            axes[1, s_idx].set_yscale('log')
            axes[1, s_idx].set_ylabel('pvals')
            #axes[1, s_idx].set_xlabel('MC step')

            # plot cuts
            axes[2, s_idx].errorbar(x, pcuts, yerr=[pcuts_lb, pcuts_ub], label=f'state-{s}', marker='o', ls='none', fillstyle='none', capsize=3)
            axes[2, s_idx].plot(x, pcuts, label=f'state-{s}', ls='-', color='k', zorder=10)
            axes[2, s_idx].set_ylabel('cuts/res')
            #axes[2, s_idx].set_xlabel('MC step')
            #axes[2].set_ylim(50,500)
            #axes[2, s_idx].set_yticks(range(50, 550, 150))
            axes[2, s_idx].grid(True, which='both', linestyle='--', linewidth=0.5)

            # plot size_means
            axes[3, s_idx].errorbar(x, size_means, yerr=[size_means_lb, size_means_ub], label=f'state-{s}', marker='o', ls='none', fillstyle='none', capsize=3)
            axes[3, s_idx].plot(x, size_means, label=f'state-{s}', ls='-', color='k', zorder=10)
            axes[3, s_idx].set_ylabel('prot. size')
            axes[3, s_idx].set_xlabel('MC step')

        for row in range(0,4):
            for col in range(0,num_states):
                if row != 2:
                    axes[row, col].sharey(axes[row, 0])

        #axes[0, 4].legend(loc='upper left', bbox_to_anchor=(1,1))
        plt.suptitle(title)
        plt.tight_layout()
        #plt.show()
        plt.savefig(outfile)
        plt.close()
        print(f'SAVED: {outfile}')

        ## save plot data
        outdfs = pd.concat(outdfs)
        outdfs.to_csv(data_outfile, index=False)
        print(f'SAVED: {data_outfile}')

    ##########################################################################################################

    ##########################################################################################################
    def plot_phaseT(self, data, title):
        ## Data is a tuple with three elements: final_genelist, final_traj, Final_step_reg
        traj_data = data[1]
        dfs = []
        for traj, traj_states in traj_data.items():
            for state, state_df in traj_states.items():
                #print(state_df)
                cols = state_df.columns
                if 'beta' not in cols or 'E' not in cols:
                    print(f'Dataframe for traj {traj} state {state} does not have beta or E columns')
                    continue
                dfs += [state_df]
        if len(dfs) == 0:
            raise ValueError('No Data to plot!')

        total_df = pd.concat(dfs)
        #print(total_df)
        total_df['T'] = 1/total_df['beta']

        # get the psuedo (Cv)p = dH/dT ~ dE/dT
        ## get the plot data
        plot_df = {'T':[], 'mean_E':[], 'E_lb_delta':[], 'E_ub_delta':[]}
        for T, Tdf in total_df.groupby('T'):
            #print(Tdf)

            E = Tdf['E'].values
            #print(Cvs)
            mean_E = np.mean(E)
            std = np.std(E)
            n = len(E)
            margin_of_error = 1.96 * (std / np.sqrt(n))

            #print(T, mean_Cv, margin_of_error)
            E_lb_delta = abs(mean_E - margin_of_error)
            E_ub_delta = abs(margin_of_error - mean_E)
            plot_df['T'] += [T]
            plot_df['mean_E'] += [mean_E]
            plot_df['E_lb_delta'] += [E_lb_delta]
            plot_df['E_ub_delta'] += [E_ub_delta]
        plot_df = pd.DataFrame(plot_df)
        #print(plot_df)
       # Extract columns
        T = plot_df['T']
        E = plot_df['mean_E']
        
        # Calculate dE/dT using numpy's gradient function
        dE_dT = np.gradient(E, T)
        
        # Create a figure with two subplots side by side
        fig, ax = plt.subplots(1, 2, figsize=(12, 5))
        
        # First subplot: E(T)
        ax[0].plot(T, E, color='blue', marker='o', linestyle='-')
        ax[0].set_title('E(T)')
        ax[0].set_xlabel('T')
        ax[0].set_ylabel('E')
        ax[0].set_xscale('log')
        
        # Second subplot: dE/dT as a function of T
        ax[1].plot(T, dE_dT, color='red', marker='x', linestyle='-')
        ax[1].set_title('dE/dT as a function of T')
        ax[1].set_xlabel('T')
        ax[1].set_ylabel('dE/dT')
        ax[1].set_xscale('log')
        ax[1].set_ylim(0,100)
        
        ## save plot data
        outdf = pd.DataFrame({'T':T, 'E':E, 'dE/dT':dE_dT})
        #print(outdf)
        outfile = os.path.join(self.outpath, f'PhaseT_{title}.csv')
        outdf.to_csv(outfile, index=False)
        print(f'SAVED: {outfile}')

        # Adjust layout for better spacing
        plt.tight_layout()
        #plt.show()
        outfile = os.path.join(self.outpath, f'PhaseT_{title}.png')
        plt.savefig(outfile)
        print(f'SAVED: {outfile}')
        
    ##########################################################################################################
    def get_ranked_state_ctables(self, ranking_summary_df, res_feat_df, final_gene_lists, ranking_summary_outfile):
        """
        For each ranked state in each traj get the contingency table of the last frame used to select genes and plot it. OR maybe different from ranking OR as those are an average over the last 10 steps. 
        """
        print(f'Getting contingency tables for ranked states')
        noReg_noCut = []
        yesReg_noCut = []
        noReg_yesCut = []
        yesReg_yesCut = []
        lastOR = []
        lastPvalue = []
        total = []
        for rowi, row in ranking_summary_df.iterrows():
            #print(row)
            gene_list = final_gene_lists[row['traj']][row['old_stateID']]
            #print(gene_list)

            loc_res_feat_df = res_feat_df[res_feat_df['gene'].isin(gene_list)]
            #print(loc_res_feat_df)

            ctable = pd.crosstab(loc_res_feat_df[self.cut_key], loc_res_feat_df['region'])
            #print(f'ctable:\n{ctable}')
            res = fisher_exact(ctable)
            OR, pvalue = res.statistic, res.pvalue
            #print(f'OR: {OR} with p-value {pvalue}')

            noReg_noCut += [ctable.values[0,0]]
            yesReg_noCut += [ctable.values[0,1]]
            noReg_yesCut += [ctable.values[1,0]]
            yesReg_yesCut += [ctable.values[1,1]]
            lastOR += [OR]
            lastPvalue += [pvalue]
            total += [np.sum(ctable.values)]
            
        res = ranking_summary_df.copy()
        res['noReg_noCut'] = noReg_noCut
        res['yesReg_noCut'] = yesReg_noCut
        res['noReg_yesCut'] = noReg_yesCut
        res['yesReg_yesCut'] = yesReg_yesCut
        res['total'] = total
        res['lastOR'] = lastOR
        res['lastPvalue'] = lastPvalue
        res.to_csv(ranking_summary_outfile, index=False)
        print(f'SAVED: {ranking_summary_outfile}')
        return res
        
    ##########################################################################################################

    ##########################################################################################################
    def run(self):
        """
        Orchestrates the workflow by loading data, performing regression, and saving results.
        """
        start_time = time.time()

        # Load data for cyto-serum only
        C_data = self.load_data('C')
        num_states = len(C_data[0]['0'])
        print(f'Number of states: {num_states}')
     
        # Plot the trajectories of each state and the final state regression data
        #self.Plot_traj(C_data, f'{self.tag}_C')

        # plot the enesmble of traj in a single plot
        #self.Plot_ensem(C_data, f'{self.tag}_C', num_states)
       
        # Get the overlap of gene histograms and CDF across the highest OR state in each traj
        self.State_overlap(C_data, f'{self.tag}_C')

        # plot the phase transition
        #self.plot_phaseT(C_data, f'{self.tag}_C')

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
    parser.add_argument("-f", "--resFeat_files", type=str, required=True, help="Path to residue feature files")
    parser.add_argument("-r", "--reg_formula", type=str, required=True, help="Regression formula")
    args = parser.parse_args()

    # Setup logging configuration
    logging.basicConfig(filename=args.log, level=logging.INFO, format='%(asctime)s %(message)s')
    logging.info(f'{"#"*100}\nNEW RUN {script_name}')

    analysis = DataAnalysis(
        inpdir=args.inpdir,
        outpath=args.outpath,
        tag=args.tag, 
        resFeat_files=args.resFeat_files,
        reg_formula=args.reg_formula)
    analysis.run()

if __name__ == "__main__":
    main()

