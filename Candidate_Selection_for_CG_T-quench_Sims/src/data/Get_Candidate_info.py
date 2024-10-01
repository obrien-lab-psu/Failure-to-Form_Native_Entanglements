import time, sys
import pickle
import multiprocessing as mp
from scipy.stats import bootstrap
import logging
import argparse
import glob
import numpy as np
import pandas as pd
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
import requests
from tqdm import tqdm
from scipy.spatial import cKDTree

pd.set_option('display.max_rows', 4000)
pd.options.mode.chained_assignment = None  # default='warn'

class DataAnalysis:
    """
    A class to handle the data analysis process including encoding, regression, and statistical tests.
    """
    ##########################################################################################################
    def __init__(self, refolded, outpath, misfolded, ent_features):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:

        """
        self.refolded = refolded
        self.outpath = outpath
        self.misfolded = misfolded
        self.ent_features = ent_features

        if not os.path.exists(f'{self.outpath}'):
            os.makedirs(f'{self.outpath}')
            print(f'Made output directories {self.outpath}')
    ##########################################################################################################

    ##########################################################################################################
    def load_data(self,):
        print(f'Loading data')

        ## load refolded genes df
        refolded_genes = pd.read_csv(self.refolded)
        refolded_genes = refolded_genes['gene'].values
        print(f'refolded_genes: {refolded_genes} {len(refolded_genes)}')

        ## load misfolded genes df and apply threshold
        misfolded_genes = pd.read_csv(self.misfolded)
        misfolded_genes = misfolded_genes[misfolded_genes['F'] >= 0.7]
        misfolded_genes = misfolded_genes['gene'].values
        print(f'misfolded_genes: {misfolded_genes} {len(misfolded_genes)}')

        ## Quality check that there is no overlap between the set of refolded and misfolded genes
        overlap = set(refolded_genes).intersection(set(misfolded_genes))
        print(f'overlap: {overlap}')
        if len(overlap) != 0:
            raise ValueError(f'There is overlap between the set of misfolded and refolded genes: {overlap}')

        ## Loading set of unique entanglement features for misfoldd proteins
        feature_files = glob.glob(os.path.join(self.ent_features, '*'))
        misfolded_features = []
        for g in misfolded_genes:
            f = [f for f in feature_files if g in f]
            df = pd.read_csv(f[0], sep='|')
            df.insert(0, 'gene', [g] * len(df))
            misfolded_features += [df]
        misfolded_features = pd.concat(misfolded_features)
        print(f'misfolded_features:\n{misfolded_features}')
        self.misfolded_features = misfolded_features


        ## Loading set of unique entanglement features for refoldable proteins if they exist
        refolded_features = []
        for g in refolded_genes:
            f = [f for f in feature_files if g in f]
            df = pd.read_csv(f[0], sep='|')
            if len(df) == 0:
                df.loc[0] = [np.nan] * len(df.columns)
                size = int(get_uniprot_sequence_length(g))
                df['prot_size'] = size
            df.insert(0, 'gene', [g] * len(df))
            refolded_features += [df]
        refolded_features = pd.concat(refolded_features)
        print(f'refolded_features:\n{refolded_features}')
        self.refolded_features =refolded_features
        
    ##########################################################################################################
    def match(self,):
        print(f'Matching 10 refolded to misfolded ')
        misfold_sizes = {'gene':[], 'prot_size':[]}
        for gene, gene_df in self.misfolded_features.groupby('gene'):
            misfold_sizes['gene'] += [gene]
            misfold_sizes['prot_size'] += [gene_df['prot_size'].values[0]]
        misfold_sizes = pd.DataFrame(misfold_sizes)
        print(f'misfold_sizes:\n{misfold_sizes}')

        refolded_sizes = {'gene':[], 'prot_size':[]}
        for gene, gene_df in self.refolded_features.groupby('gene'):
            refolded_sizes['gene'] += [gene]
            if np.isnan(gene_df['prot_size'].values[0]):
                size = int(get_uniprot_sequence_length(gene))
                refolded_sizes['prot_size'] += [size]
            else:
                refolded_sizes['prot_size'] += [gene_df['prot_size'].values[0]]
        refolded_sizes = pd.DataFrame(refolded_sizes)
        refolded_sizes['prot_size'] = refolded_sizes['prot_size'].astype(int)
        print(f'refolded_sizes:\n{refolded_sizes}')

        ## get 10 random indexes and those will be the 10 random rows to pick from the matched dataset
        if len(refolded_sizes) > 10:
            refolded_sizes = refolded_sizes.sample(10)
        print(f'refolded_sizes:\n{refolded_sizes}')

        ## get the closest matches
        matched_refolded, matched_misfolded = find_closest_matches(refolded_sizes, misfold_sizes)
        print(f'matched_refolded:\n{matched_refolded}')
        print(f'matched_misfolded:\n{matched_misfolded}')

        self.matched_refolded = self.refolded_features[self.refolded_features['gene'].isin(matched_refolded['gene'].values)]
        self.matched_misfolded = self.misfolded_features[self.misfolded_features['gene'].isin(matched_misfolded['gene'].values)]
        print(f'matched_refolded:\n{self.matched_refolded}')
        print(f'matched_misfolded:\n{self.matched_misfolded}')  
    ##########################################################################################################

        
    ##########################################################################################################

    ##########################################################################################################
    def run(self):
        """
        Orchestrates the workflow by loading data, performing regression, and saving results.
        """
        start_time = time.time()

        self.load_data()
        self.match()

        logging.info(f'NORMAL TERMINATION {time.time() - start_time}')
    ##########################################################################################################

def get_uniprot_sequence_length(accession_id):
    """
    Fetches the canonical protein sequence length for a given UniProt accession ID.
    
    Args:
    accession_id (str): The UniProt accession ID.
    
    Returns:
    int: The length of the canonical protein sequence, or None if not found.
    """
    # UniProt REST API URL
    url = f"https://rest.uniprot.org/uniprotkb/{accession_id}.fasta"
    
    try:
        # Make the GET request to fetch the FASTA format data
        response = requests.get(url)
        response.raise_for_status()  # Raise an error for a failed request

        # Extract the protein sequence from the response (ignoring the header lines)
        fasta_data = response.text.split('\n')
        sequence = ''.join(line for line in fasta_data if not line.startswith('>'))

        # Return the length of the sequence
        return len(sequence)
    
    except requests.exceptions.HTTPError as http_err:
        print(f"HTTP error occurred: {http_err}")
        return None
    except Exception as err:
        print(f"An error occurred: {err}")
        return None

def find_closest_matches(df1, df2, num_matches=10):
    # Make a copy of df2 so we can modify it without affecting the original
    df1 = df1.sort_values(by=['prot_size'])
    df1 = df1.reset_index(drop=True)
    df2 = df2.sort_values(by=['prot_size'])
    df2 = df2.reset_index(drop=True)
    df2_copy = df2.copy()
    #print(f'df1:\n{df1}')
    #print(f'df2_copy:\n{df2_copy}')

    matched_df2 = {'gene':[], 'prot_size':[]}
    for size in df1['prot_size'].values:
        df2_copy['temp'] = df2_copy['prot_size'] - size
        df2_copy['temp'] = np.abs(df2_copy['temp'])
        minidx = df2_copy['prot_size'].idxmin()
        match = df2_copy.loc[minidx]

        matched_df2['gene'] += [match['gene']]
        matched_df2['prot_size'] += [match['prot_size']]
        df2_copy = df2_copy.drop(minidx)
    matched_df2 = pd.DataFrame(matched_df2)
    #print(f'matched_df2:\n{matched_df2}')
    return df1, matched_df2



def main():
    """
    Main function to parse arguments and run the DataAnalysis class.
    """
    script_name = f'Optimizer_SimulatedAnnealing'
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-r", "--refolded", type=str, required=True, help="Path to top refolded gene file")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("-m", "--misfolded", type=str, required=True, help="path to misfolded genes cdf file")
    parser.add_argument("-l", "--log", type=str, required=True, help="Path to logging file")
    parser.add_argument("-e", "--ent_features", type=str, required=True, help="path to unique ent feature file")
    args = parser.parse_args()

    # Setup logging configuration
    logging.basicConfig(filename=args.log, level=logging.INFO, format='%(asctime)s %(message)s')
    logging.info(f'{"#"*100}\nNEW RUN {script_name}')

    analysis = DataAnalysis(
        refolded=args.refolded,
        outpath=args.outpath,
        misfolded=args.misfolded,
        ent_features=args.ent_features)
    analysis.run()

if __name__ == "__main__":
    main()

