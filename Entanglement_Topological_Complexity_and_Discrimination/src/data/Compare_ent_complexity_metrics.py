import sys, os, re, time, logging
from sklearn.utils import shuffle
import ast
from collections import defaultdict
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LassoCV, Lasso
from sklearn.datasets import make_regression
from sklearn.model_selection import StratifiedKFold, KFold, cross_validate, GridSearchCV
from sklearn.metrics import mean_squared_error
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score, balanced_accuracy_score, average_precision_score,f1_score,recall_score,precision_score,roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn import metrics, preprocessing
from sklearn.preprocessing import LabelEncoder
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

    def __init__(self, outpath, uent_files, buff, spa):
        """
        Initializing the DataAnalysis object and ensure atleast the top level outpath is present and makes it if not. 
        """

        # Make outpath if it doesnt exists
        self.outpath = outpath
        self.DataAnalysisOutpath = os.path.join(self.outpath, 'DataAnalysisOutput/')
        if not os.path.exists(self.DataAnalysisOutpath):
            os.makedirs(self.DataAnalysisOutpath)
            print(f'Made directory: {self.DataAnalysisOutpath}')

        self.DataAnalysisDistPlotsOutpath = os.path.join(self.outpath, 'DataAnalysisOutput/Dist_plots/')
        if not os.path.exists(self.DataAnalysisDistPlotsOutpath):
            os.makedirs(self.DataAnalysisDistPlotsOutpath)
            print(f'Made directory: {self.DataAnalysisDistPlotsOutpath}')

        self.three_to_one_letter = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'MSE': 'M', 'PHE': 'F', 
        'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 
        'VAL': 'V'}
        
        self.uent_files = glob(os.path.join(uent_files, '*'))

        self.keys = ['Gn', 'N_term_thread', 'Gc', 'C_term_thread',
        'loopsize', 'num_zipper_nc', 'perc_bb_loop',
        'num_loop_contacting_res', 'num_cross_nearest_neighbors',
        'ent_coverage', 'min_N_prot_depth_left', 'min_N_thread_depth_left',
        'min_N_thread_slippage_left', 'min_C_prot_depth_right',
        'min_C_thread_depth_right', 'min_C_thread_slippage_right', 'ACO', 'RCO']

        self.buff = buff
        self.spa = spa

    #################################################################################################################
    def load_files(self, mask):
        """
        Method to load unique entanglement files and keep only certain columns for analysis
        'Gn', 
        'N_term_thread', 
        'Gc', 
        'C_term_thread',
        'loopsize', 
        'num_zipper_nc', 
        'perc_bb_loop',
        'num_loop_contacting_res', 
        'num_cross_nearest_neighbors',
        'ent_coverage', 
        'min_N_prot_depth_left', 
        'min_N_thread_depth_left',
        'min_N_thread_slippage_left', 
        'min_C_prot_depth_right',
        'min_C_thread_depth_right', 
        'min_C_thread_slippage_right', 
        'prot_size'

        """
        dfs = []
        for gene in mask:
            print(gene)
            gene_uent_file = [f for f in self.uent_files if gene in f][0]

            gene_uent = pd.read_csv(gene_uent_file, sep='|')
            #print(gene, gene_uent_file)
            #print(gene_uent)
            num_uent = len(gene_uent)

            if num_uent != 0:
                gene_uent['gene'] = gene
                dfs += [gene_uent]

        uent_df = pd.concat(dfs) 
        return uent_df


    def DistStats(self, df, keys, dist_tag, n_resamples=1000, alpha=0.05):
        """
        Calculate various parameters of the data distributions
        (1) mean, median, mode
        (2) confidence intervals for mean, median, and mode
        * change the self.keys parameters to suit your needs
        """
        df_copy = df.copy()
        df_copy = df_copy[keys]
        print(f'df_copy:\n{df_copy}')

        results = {'metric':[], 'mean':[], 'mean_lb':[], 'mean_ub':[], 'median':[], 'median_lb':[], 'median_ub':[]}

        def calculate_bootstrap_ci(data, statistic_func):
            if statistic_func == 'mean':
                res = bootstrap((data,), np.mean, confidence_level=1-alpha, n_resamples=n_resamples)
                return res.confidence_interval.low, res.confidence_interval.high

            elif statistic_func == 'median':
                medians = []
                for b in range(n_resamples):
                    boot = np.random.choice(data, replace=True)
                    medians += [np.median(boot)]
                lb = np.percentile(medians, 2.5)
                ub = np.percentile(medians, 97.5)
                return lb, ub


        for column in df_copy.columns:
            print(f'column: {column}')
            col_data = df_copy[column].dropna().values

            ## drop 0s if column == C_term_thread or N_term_thread
            if column in ['N_term_thread', 'C_term_thread']:
                col_data = col_data[np.where(col_data != 0)]
            print(col_data)

            mean_val = np.mean(col_data)
            median_val = np.median(col_data)

            mean_ci = calculate_bootstrap_ci(col_data, 'mean')
            print('Mean', mean_val, mean_ci)
            median_ci = calculate_bootstrap_ci(col_data, 'median')
            print('Median', median_val, median_ci)

            results['metric'] += [column]
            results['mean'] += [mean_val]
            results['mean_lb'] += [mean_ci[0]]
            results['mean_ub'] += [mean_ci[1]]
            results['median'] += [median_val]
            results['median_lb'] += [median_ci[0]]
            results['median_ub'] += [median_ci[1]]

            # plot histogram
            plot_filename = f'{self.DataAnalysisDistPlotsOutpath}{dist_tag}_{column}_{self.buff}_{self.spa}.png'
            plt.hist(col_data, bins=100, color='blue', edgecolor='black', density=True)  # 100 bins
            plt.xlabel(column)
            plt.ylabel('PDF')
            plt.title(dist_tag)
            plt.savefig(plot_filename)
            print(f'SAVED: {plot_filename}')
            plt.close()

        return pd.DataFrame(results)


    def Permutation(self, df1_full, df2_full, keys, n_resamples=10000):
        """
        For the columns in the two dataframes calculate the pairwise pvalue by permutation. 
        """

        df1 = df1_full.copy()
        df2 = df2_full.copy()
        df1 = df1[keys]
        df2 = df2[keys]

        if not all(df1.columns == df2.columns):
            raise ValueError("Both DataFrames must have identical column names")

        results = []

        for column in df1.columns:
            data1 = df1[column].dropna().values
            data2 = df2[column].dropna().values

            # Define the statistic function for the permutation test
            def statistic(x, y, axis):
                return np.mean(x, axis=axis) - np.mean(y, axis=axis)

            # Perform the permutation test
            #res = permutation_test((data1, data2), statistic, vectorized=True, n_resamples=n_resamples)
            stat = statistic(data1, data2, 0)
            if stat > 0:
                res = mannwhitneyu(data1, data2, alternative='greater')
            elif stat < 0:
                res = mannwhitneyu(data1, data2, alternative='less')
            else:
                res = mannwhitneyu(data1, data2)

            #results[column] = res.pvalue
            results += [res.pvalue]

            print(f'column: {res.pvalue}')

        return results

        
    def LassoRegression(self, df, X_keys, Y_key):

        ### Get input data 
        X = df[X_keys]
        X = X.fillna(0)
        # Initialize the LabelEncoder and make sure y categories are integer incoded
        y = df[Y_key]
        label_encoder = LabelEncoder()
        y_encoded = label_encoder.fit_transform(y)

        scaler = preprocessing.StandardScaler()
        std_scaled_df = scaler.fit_transform(X)
        std_scaled_df = pd.DataFrame(std_scaled_df, columns=X_keys)
        X = std_scaled_df
        X = X.values
        y = y.values
        #print(f'X:\n{X}')
        #print(f'y:\n{y}')

        logistic_regression =  LogisticRegression(penalty='l1', solver='liblinear')
        Cs = np.linspace(0.00001, 10, 1000)
        #Cs = np.linspace(0.00001, 10, 10)

        fit_data = {'C':[], 'fold':[], 'balanced_accuracy':[], 'accuracy':[]}
        for C in Cs:

            print(f'{"#"*100}\nTESTING C: {C}')

            for col in X_keys:
                if col not in fit_data:
                    fit_data[col] = []

            ## make folds and fit model
            skf = StratifiedKFold(n_splits=5, shuffle=True)
            #skf = StratifiedKFold(n_splits=5)

            for i, (train_index, test_index) in enumerate(skf.split(X, y)):
                #print(f"Fold {i}:")
                #print(f"  Train: index={train_index} {len(train_index)}")
                #print(f"  Test:  index={test_index} {len(test_index)}")

                X_train = X[train_index]
                y_train = y[train_index]
                X_test = X[test_index]
                y_test = y[test_index]

                # Get features for optimal regularization ceof
                logistic_regression =  LogisticRegression(penalty='l1', solver='liblinear', C=C)
                logistic_regression.fit(X_train, y_train)

                coefs = logistic_regression.coef_[0].tolist()
                #print(f'coefs: {coefs}')

                # Predict on the testing data
                y_pred = logistic_regression.predict(X_test)
                #print(f'y_pred: {y_pred}')
                # Calculate balanced accuracy
                balanced_accuracy = balanced_accuracy_score(y_test, y_pred)
                accuracy = accuracy_score(y_test, y_pred)

                fit_data['C'] += [C]
                fit_data['fold'] += [i]
                fit_data['balanced_accuracy'] += [balanced_accuracy]
                fit_data['accuracy'] += [accuracy]
                for col_i, col in enumerate(X_keys):
                    fit_data[col] += [coefs[col_i]]
                    
        fit_data = pd.DataFrame(fit_data)
        print(f'fit_data:\n{fit_data}')
        return fit_data

    def Plot_Lasso(self, df, keys, outfile):
        
        #        C  fold  balanced_accuracy  accuracy  
        # Create subplots
        fig, axes = plt.subplots(1, 3, figsize=(12, 4))
        X = []
        num_nonzero = []
        num_robust_nonzero = []
        BA = []
        for C, C_df in df.groupby('C'):
            #print(C_df)

            num_nonzero_coef = []
            num_robust_nonzero_coef = []
            for key in keys:
                coefs = C_df[key].values
                all_nonzero = np.all(coefs != 0)
                same_sign = np.all(coefs >= 0) or np.all(coefs <= 0)
                #print(C, key, coefs, all_nonzero, same_sign)
                num_nonzero_coef += [all_nonzero]
                num_robust_nonzero_coef += [same_sign]

            X += [C]
            num_nonzero += [np.sum(num_nonzero_coef)]
            num_robust_nonzero += [np.sum(num_robust_nonzero_coef)]
            BA += [np.mean(C_df['balanced_accuracy'].values)]

        axes[0].plot(X, num_nonzero)
        axes[0].set_ylabel('# non-zero ceof.')
        axes[0].set_xlabel('inverse regularization strength')
        axes[0].set_ylim(0, 20)
        axes[1].plot(X, num_robust_nonzero)
        axes[1].set_ylabel('# non-zero & robust ceof.')
        axes[1].set_xlabel('inverse regularization strength')
        axes[1].set_ylim(0, 20)
        axes[2].plot(X, BA)
        axes[2].set_ylabel('Balanced Accuracy')
        axes[2].set_xlabel('inverse regularization strength')
        axes[2].set_ylim(0, 1)

        plt.tight_layout()
        plt.savefig(outfile)
        print(f'SAVED: {outfile}')

#################################################################################################################

def main():
    """
    Main function to control workflow. 
    """

    # Parse the user supplied arguments
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-Eg", "--Ess_gene_list", type=str, required=True, help=f"path to Essential gene list used for mask in Fc_ab and FcG_ab calcs")
    parser.add_argument("-NEg", "--NonEss_gene_list", type=str, required=True, help=f"path to Non-Essential gene list used for mask in Fc_ab and FcG_ab calcs")
    parser.add_argument("-l", "--log_file", type=str, required=True, help="Path to logging file")
    parser.add_argument("-e", "--uent_files", type=str, required=True, help="path to unique entanglement files")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="path to output directory. will be made if doesnt exist")
    parser.add_argument("-p", "--num_permute", type=int, required=True, help="Number of permutations")
    parser.add_argument("-b", "--buffer", type=str, required=True, help="Buffer system to use")
    parser.add_argument("-s", "--spa", type=str, required=True, help="SPA threshold")
    args = parser.parse_args()

    uent_files = args.uent_files
    log_file = args.log_file
    outpath = args.outpath
    Ess_gene_list = args.Ess_gene_list
    NonEss_gene_list = args.NonEss_gene_list
    num_permute = args.num_permute
    buff = args.buffer
    spa = args.spa

    # Setup logging configuration
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s %(message)s') 
    logging.info(f'{"#"*100}\nNEW RUN')

    # Make outpath if it doesnt exists
    if not os.path.exists(outpath):
        os.makedirs(outpath)
        print(f'Made directory: {outpath}')

    # Initalize the DataAnalysis class object
    Analyzer = DataAnalysis(outpath, uent_files, buff, spa)
    print(f'Analyzer: {Analyzer}')

    ## Get the Essential gene data and stats
    Ess_gene_mask = np.loadtxt(Ess_gene_list, dtype=str)
    print(f'Ess_gene_mask: {Ess_gene_mask} {len(Ess_gene_mask)}')

    Ess_ent_data = Analyzer.load_files(Ess_gene_mask)
    Ess_ent_data['Gn'] = Ess_ent_data['Gn'].abs()
    Ess_ent_data['Gn'] = Ess_ent_data['Gn'].apply(lambda x: np.nan if x < 0.6 else x)
    Ess_ent_data['Gc'] = Ess_ent_data['Gc'].abs()
    Ess_ent_data['Gc'] = Ess_ent_data['Gc'].apply(lambda x: np.nan if x < 0.6 else x)
    Ess_ent_data['Essential'] = 1
    print(f'Ess_ent_data: {Ess_ent_data}')
    Ess_ent_data_outfile = f'{Analyzer.outpath}Ess_uent_data_{buff}_{spa}.csv'
    Ess_ent_data.to_csv(Ess_ent_data_outfile, sep='|', index=False)
    print(f'SAVED: {Ess_ent_data_outfile}')

    Ess_stats = Analyzer.DistStats(Ess_ent_data, Analyzer.keys, 'Ess', n_resamples=num_permute)
    print(f'Ess_stats:\n{Ess_stats}')
    Ess_stats_data_outfile = f'{Analyzer.outpath}Ess_stats_uent_data_{buff}_{spa}.csv'
    Ess_stats.to_csv(Ess_stats_data_outfile, sep='|', index=False)
    print(f'SAVED: {Ess_stats_data_outfile}')


    ## Get the Non-essential gene data and stats
    NonEss_gene_mask = np.loadtxt(NonEss_gene_list, dtype=str)
    print(f'NonEss_gene_mask: {NonEss_gene_mask} {len(NonEss_gene_mask)}')

    NonEss_ent_data = Analyzer.load_files(NonEss_gene_mask)
    NonEss_ent_data['Gn'] = NonEss_ent_data['Gn'].abs()
    NonEss_ent_data['Gn'] = NonEss_ent_data['Gn'].apply(lambda x: np.nan if x < 0.6 else x)
    NonEss_ent_data['Gc'] = NonEss_ent_data['Gc'].abs()
    NonEss_ent_data['Gc'] = NonEss_ent_data['Gc'].apply(lambda x: np.nan if x < 0.6 else x)
    NonEss_ent_data['Essential'] = 0
    print(f'NonEss_ent_data: {NonEss_ent_data}')
    NonEss_ent_data_outfile = f'{Analyzer.outpath}NonEss_uent_data_{buff}_{spa}.csv'
    NonEss_ent_data.to_csv(NonEss_ent_data_outfile, sep='|', index=False)
    print(f'SAVED: {NonEss_ent_data_outfile}')

    NonEss_stats = Analyzer.DistStats(NonEss_ent_data, Analyzer.keys, 'NonEss', n_resamples=num_permute)
    print(f'NonEss_stats:\n{NonEss_stats}')
    NonEss_stats_data_outfile = f'{Analyzer.outpath}NonEss_stats_uent_data_{buff}_{spa}.csv'
    NonEss_stats.to_csv(NonEss_stats_data_outfile, sep='|', index=False)
    print(f'SAVED: {NonEss_stats_data_outfile}')


    ## Compare the entanglement complexity between Essential and Nonessential data
    EssVNonEss_pvalues = Analyzer.Permutation(Ess_ent_data, NonEss_ent_data, Analyzer.keys, n_resamples=num_permute)
    print(f'EssVNonEss_pvalues:\n{EssVNonEss_pvalues}')


    # Merge the DataFrames along axis=1
    metrics = Ess_stats['metric']
    Ess_stats.drop(columns=['metric'])
    NonEss_stats.drop(columns=['metric'])

    Ess_stats = Ess_stats.add_prefix('Ess_')
    NonEss_stats = NonEss_stats.add_prefix('NonEss_')

    merged_df = pd.concat([metrics, Ess_stats, NonEss_stats], axis=1)
    merged_df['pvalues'] = EssVNonEss_pvalues
    print(f'merged_df:\n{merged_df}')
    merged_data_outfile = f'{Analyzer.outpath}merged_stats_uent_data_{buff}_{spa}.csv'
    merged_df.to_csv(merged_data_outfile, sep='|', index=False)
    print(f'SAVED: {merged_data_outfile}')

    # make the combined raw ent dataframes
    combined_df = pd.concat([Ess_ent_data, NonEss_ent_data])
    print(f'combined_df:\n{combined_df}')
    combined_data_outfile = f'{Analyzer.outpath}combined_uent_data_{buff}_{spa}.csv'
    combined_df.to_csv(combined_data_outfile, sep='|', index=False)
    print(f'SAVED: {combined_data_outfile}')


    Lasso_results = Analyzer.LassoRegression(combined_df, Analyzer.keys, 'Essential')
    print(f'Lasso_results:\n{Lasso_results}')
    Lasso_results_outfile = f'{Analyzer.outpath}Lasso_results_{buff}_{spa}.csv'
    Lasso_results.to_csv(Lasso_results_outfile, sep='|', index=False)
    print(f'SAVED: {Lasso_results_outfile}')

    Lasso_results_plot_outfile = f'{Analyzer.outpath}Lasso_results_{buff}_{spa}.png'
    Analyzer.Plot_Lasso(Lasso_results, Analyzer.keys, Lasso_results_plot_outfile)

if __name__ == "__main__":
    main()

print('NORMAL TERMINATION')
