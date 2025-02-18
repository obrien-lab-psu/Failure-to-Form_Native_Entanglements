import time, sys
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
import itertools
import os
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import poisson, binom, fisher_exact, chi2, norm
import scipy.stats as st
from matplotlib.ticker import MultipleLocator
import random
import warnings
warnings.simplefilter(action='ignore', category=Warning)

#pd.set_option('display.max_rows', 4000)

class DataAnalysis:
    """
    A class to handle the data analysis process including encoding, regression, and statistical tests.
    """

    def __init__(self, resFeat_files, outpath, gene_list, sub_gene_list, tag, load_style, buffer, spa, cov, reg_formula, hold_var):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
        - resFeat_files (str): Path to residue feature files.
        - outpath (str): Path to the output directory.
        - gene_lists (str): Path to gene lists to use.
        - tag (str): Tag for output filenames.
        - load_style (str): Load style (True: load by gene, False: load a single file with all genes present).
        - buffer (str): Buffer system to use.
        - spa (str): SPA threshold.
        - cov (str): LiPMS cov threshold.
        - reg_formula (str): Regression formula.
        """
        self.resFeat_files = resFeat_files
        self.outpath = outpath
        self.gene_list = gene_list
        self.tag = tag
        self.load_style = load_style
        self.buffer = buffer
        self.spa = spa
        self.cov = cov
        self.reg_formula = reg_formula
        self.data = {}
        #self.logger = self.setup_logging()
        #self.gene_list_files = glob.glob(self.gene_list)
        self.hold_var = hold_var

        if not os.path.exists(f'{self.outpath}'):
            os.makedirs(f'{self.outpath}')
            print(f'Made output directories {self.outpath}')


    def holdouts(self, df):
        """
        Performs quasi-binomial regression analysis on the provided DataFrame.

        Parameters:
        - df (pd.DataFrame): DataFrame containing the data for regression.
        - formula (str): The formula specifying the regression model.

        Returns:
        - table_1_df (pd.DataFrame): DataFrame containing the regression results with p-values.
        """
        print(df)
        genes = df['gene'].unique()
        print(f'# genes: {len(genes)}')

        combinations_of_5 = list(itertools.combinations(genes, 5))
        print(f"Total unique combinations: {len(combinations_of_5)}")
 
        outstrs = []
        for i, combo in enumerate(combinations_of_5):
            outstrs += [','.join([g for g in genes if g not in combo])]
            if i == 0:
                continue

            if i % 270000 == 0:
                print(f'Combo: {i} {combo}')
                outstrs = np.asarray(outstrs, dtype=str)
                np.savetxt(f'{self.outpath}/holdouts_{self.tag}_{i}.txt', outstrs, fmt='%s')
                print(f'Saved holdouts_{self.tag}_{i}.txt')
                outstrs = []
                
 

    def load_data(self):
        """
        Loads the residue feature files and filters the data for analysis.
        """
        if self.load_style == 'True':
            res_files = glob.glob(self.resFeat_files)
            print(f'Number of res_files: {len(res_files)}')
            for i, gene in enumerate(self.reg_genes):
                gene_resFeat = [f for f in res_files if gene in f]
                if len(gene_resFeat) == 0:
                    print(f"No residue feature file found for gene {gene}")
                    continue
                elif len(gene_resFeat) > 1:
                    print(f"More than 1 residue feature file found for gene {gene}")
                    continue
                gene_resFeat_file = gene_resFeat[0]
                #print(f'gene_resFeat_file: {gene_resFeat_file} {i}')
                if len(self.data) == 0:
                    self.data = pd.read_csv(gene_resFeat_file, sep='|')
                else:
                    self.data = pd.concat((self.data, pd.read_csv(gene_resFeat_file, sep='|')))
        else:
            self.data = pd.read_csv(self.resFeat_files, sep='|')

        self.data = self.data[self.data['gene'].isin(self.reg_genes)]
        self.data = self.data[self.data['AA'] != 'NC']
        self.data = self.data[self.data['mapped_resid'].notna()]
        self.data = self.data[self.data['AA'].notna()]
        self.data = self.data.reset_index()
        print(f"Data loaded and filtered. Number of unique genes: {len(self.data['gene'].unique())}")


    def run(self):
        """
        Orchestrates the workflow by loading data, performing regression, and saving results.
        """
        start_time = time.time()

        # Get list of genes to select for in regression
        print(f'reg_gene_list: {self.gene_list}')
        reg_genes = np.loadtxt(self.gene_list, dtype=str)
        self.reg_genes = reg_genes
        try:
            num_genes = len(reg_genes)
        except:
            print(f'Failed to get number of genes in {self.gene_list}')
            quit()
        print(f'Number of genes in list: {len(reg_genes)}')

        if len(reg_genes) == 0:
            print(f"Gene list is empty: {self.gene_list} {len(reg_genes)}")
            quit()

        # Load data
        self.load_data()
        
        # Get data for gene list requested by user
        self.data = self.data[self.data['gene'].isin(self.reg_genes)]

        # Quality check to ensure all columns are present in the df
        reg_vars = [v for v in self.reg_formula.split(' ') if v not in ['~', '+']]
        print(reg_vars)
        keys = []
        for reg_var in reg_vars:
            if '*' in reg_var:
                reg_var = reg_var.split('*')
            else:
                reg_var = [reg_var]
            for v in reg_var:
                keys += [v]

                if 'cut_' in v:
                    cut_key = v

        keys = list(set(keys))
        bin_keys = keys.copy()
        self.cut_key = cut_key
        if 'AA' in keys:
            bin_keys.remove('AA')
        print(f'bin_keys: {bin_keys}')
        self.bin_keys = bin_keys
        self.reg_vars = reg_vars


        # Perform regression
        self.holdouts(self.data)

        print(f'NORMAL TERMINATION: {time.time() - start_time} seconds')



def main():
    """
    Main function to parse arguments and run the DataAnalysis class.
    """
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-f", "--resFeat_files", type=str, required=True, help="Path to residue feature files")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("-g", "--gene_list", type=str, required=True, help="Path to gene list to use")
    parser.add_argument("-sub_g", "--sub_gene_list", type=str, required=True, help="Path to holdout gene list to use")
    parser.add_argument("-t", "--tag", type=str, required=True, help="Tag for output filenames")
    parser.add_argument("-l", "--load_style", type=str, required=True, help="Load style (True: load by gene, False: load a single file with all genes present)")
    parser.add_argument("-b", "--buffer", type=str, required=True, help="Buffer system to use")
    parser.add_argument("-s", "--spa", type=str, required=True, help="SPA threshold")
    parser.add_argument("-c", "--cov", type=str, required=True, help="LiPMS cov threshold")
    parser.add_argument("-r", "--reg_formula", type=str, required=True, help="Regression formula")
    parser.add_argument("-v", "--hold_var", type=str, required=True, help="Variable to hold constant while calculating odds")

    args = parser.parse_args()

    analysis = DataAnalysis(
        resFeat_files=args.resFeat_files,
        outpath=args.outpath,
        gene_list=args.gene_list,
        sub_gene_list=args.sub_gene_list,
        tag=args.tag,
        load_style=args.load_style,
        buffer=args.buffer,
        spa=args.spa,
        cov=args.cov,
        reg_formula=args.reg_formula,
        hold_var=args.hold_var)
    analysis.run()

if __name__ == "__main__":
    main()

