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
    def __init__(self, refolded, outpath, misfolded, ent_features, res_features):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:

        """
        self.refolded = refolded
        self.outpath = outpath
        self.misfolded = misfolded
        self.ent_features = ent_features
        self.res_features = res_features

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
            print(f'WARNING: There is overlap between the set of misfolded and refolded genes: {overlap}. Check to ensure these are only misfolded at the 1min timepoint')

        ## Loading set of unique entanglement features for misfoldd proteins
        feature_files = glob.glob(os.path.join(self.ent_features, '*'))
        res_feature_files = glob.glob(os.path.join(self.res_features, '*'))
        misfolded_features = []
        for g in misfolded_genes:
            f = [f for f in feature_files if g in f]
            rf = [f for f in res_feature_files if g in f]
            #print(g, f, rf)
            df = pd.read_csv(f[0], sep='|')
            df.insert(0, 'gene', [g] * len(df))
            
            res_feat_df = pd.read_csv(rf[0], sep='|')
            res_feat_df = res_feat_df[~res_feat_df['mapped_resid'].isna()]
            
            num_mapped = len(res_feat_df)
            num_unmapped = len(pd.read_csv(rf[0], sep='|')) - num_mapped
            #print(num_mapped, num_unmapped)
            df['#-mapped'] = num_mapped
            df['#-unmapped'] = num_unmapped
            df['#-missing'] = df['prot_size'] - num_mapped

            misfolded_features += [df]
        misfolded_features = pd.concat(misfolded_features)
        print(f'misfolded_features:\n{misfolded_features}')
        self.misfolded_features = misfolded_features

       
        ## Loading set of unique entanglement features for refoldable proteins if they exist
        #feature_files = glob.glob(os.path.join(self.ent_features.replace('EXP', 'AF'), '*'))
        refolded_features = []
        for g in refolded_genes:
            f = [f for f in feature_files if g in f]
            rf = [f for f in res_feature_files if g in f]
            print(g,f)
            df = pd.read_csv(f[0], sep='|')
            res_feat_df = pd.read_csv(rf[0], sep='|')
            if len(df) == 0:
                df.loc[0] = [np.nan] * len(df.columns)
                size = int(get_uniprot_sequence_length(g))
                df['prot_size'] = size
                pdb = f[0].split('/')[-1].split('_')[1]
                chain = f[0].split('/')[-1].split('_')[2]
                df['PDB'] = pdb
                df['chain'] = chain
            df.insert(0, 'gene', [g] * len(df))

            res_feat_df = pd.read_csv(rf[0], sep='|')
            res_feat_df = res_feat_df[~res_feat_df['mapped_resid'].isna()]
            
            num_mapped = len(res_feat_df)
            num_unmapped = len(pd.read_csv(rf[0], sep='|')) - num_mapped
            #print(num_mapped, num_unmapped)
            df['#-mapped'] = num_mapped
            df['#-unmapped'] = num_unmapped
            df['#-missing'] = df['prot_size'] - num_mapped

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
        print(misfold_sizes['prot_size'].describe())

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
        print(refolded_sizes['prot_size'].describe())

        ## get 10 random indexes and those will be the 10 random rows to pick from the matched dataset
        if len(refolded_sizes) > 10:
            refolded_sizes = refolded_sizes.sample(10)
        refolded_sizes = refolded_sizes.reset_index(drop=True)
        print(f'refolded_sizes:\n{refolded_sizes}')

        ## get 10 random sample from the misfolded dataset to be closer to representing the proteome
        unmatched_misfolded_sizes = misfold_sizes.sample(10)
        print(f'unmatched_misfolded_sizes\n:{unmatched_misfolded_sizes}')

        ## get the closest matches
        matched_refolded, matched_misfolded = find_closest_matches(refolded_sizes, misfold_sizes)
        print(f'matched_refolded:\n{matched_refolded}')
        print(f'matched_misfolded:\n{matched_misfolded}')
    
        self.matched_refolded = self.refolded_features[self.refolded_features['gene'].isin(matched_refolded['gene'].values)]
        self.matched_misfolded = self.misfolded_features[self.misfolded_features['gene'].isin(matched_misfolded['gene'].values)]
        self.unmatched_misfolded = self.misfolded_features[self.misfolded_features['gene'].isin(unmatched_misfolded_sizes['gene'].values)]

        ## replace any ',' with ', ' in the unmapped-NC and unmapped-crossings columns
        for col in ['unmapped-NC', 'unmapped-crossings', 'mapped-NC', 'mapped-crossings']:
            self.matched_refolded[col] = self.matched_refolded[col].astype(str)
            self.matched_misfolded[col] = self.matched_misfolded[col].astype(str)
            self.unmatched_misfolded[col] = self.unmatched_misfolded[col].astype(str)
            self.matched_refolded[col] = [s.replace(',', ', ') for s in self.matched_refolded[col].values]
            self.matched_misfolded[col] = [s.replace(',', ', ') for s in self.matched_misfolded[col].values]
            self.unmatched_misfolded[col] = [s.replace(',', ', ') for s in self.unmatched_misfolded[col].values]

        ## drop the unmapped-NC_wbuff and unmapped-crossings_wbuff columns
        cols = ['unmapped-NC_wbuff', 'unmapped-crossings_wbuff']
        self.matched_refolded = self.matched_refolded.drop(cols, axis=1)
        self.matched_misfolded = self.matched_misfolded.drop(cols, axis=1)
        self.unmatched_misfolded = self.unmatched_misfolded.drop(cols, axis=1)

        ## move prot_size to 3 position
        self.matched_refolded.insert(3, 'prot_size', self.matched_refolded.pop('prot_size'))
        self.matched_misfolded.insert(3, 'prot_size', self.matched_misfolded.pop('prot_size'))
        self.unmatched_misfolded.insert(3, 'prot_size', self.unmatched_misfolded.pop('prot_size'))

        ## sort by prot_size
        self.matched_refolded = self.matched_refolded.sort_values(by=['prot_size'])
        self.matched_misfolded = self.matched_misfolded.sort_values(by=['prot_size'])
        self.unmatched_misfolded = self.unmatched_misfolded.sort_values(by=['prot_size'])

        ## rename columns
        rename_map = {'prot_size':'length', 'N_term_thread':'#N-cross', 'C_term_thread':'#C-cross', 'unmapped-NC':'NC', 'unmapped-crossings':'cross', 'num_zipper_nc':'#Zip-NC', 'perc_bb_loop':'%BB-loop',
                    'num_loop_contacting_res':'#Loop-Cont.', 'num_cross_nearest_neighbors':'#Cross-Cont.', 'ent_coverage':'EntCov', 'min_N_prot_depth_left':'NprotDepth', 'min_C_prot_depth_right':'CprotDepth',
                    'min_N_thread_depth_left':'NthreadDepth', 'min_N_thread_slippage_left':'NthreadSlip', 'min_C_thread_depth_right':'CthreadDepth', 'min_C_thread_slippage_right':'CthreadSlip'}
        self.matched_refolded = self.matched_refolded.rename(rename_map, axis='columns')
        self.matched_misfolded = self.matched_misfolded.rename(rename_map, axis='columns')
        self.unmatched_misfolded = self.unmatched_misfolded.rename(rename_map, axis='columns')
        
        ## save matched refolded df        
        print(f'matched_refolded:\n{self.matched_refolded}')
        print(self.matched_refolded.describe())
        print(f'matched_misfolded:\n{self.matched_misfolded}')  
        print(self.matched_misfolded.describe())
        print(f'unmatched_misfolded:\n{self.unmatched_misfolded}')
        print(self.unmatched_misfolded.describe())
        matched_refolded_outpath = os.path.join(self.outpath, 'Matched_refolded.csv')
        matched_misfolded_outpath = os.path.join(self.outpath, 'Matched_misfolded.csv')
        unmatched_misfolded_outpath = os.path.join(self.outpath, 'UnMatched_misfolded.csv')
        self.matched_refolded.to_csv(matched_refolded_outpath, index=False)
        self.matched_misfolded.to_csv(matched_misfolded_outpath, index=False)
        self.unmatched_misfolded.to_csv(unmatched_misfolded_outpath, index=False)
        print(f'SAVED: {matched_refolded_outpath}')
        print(f'SAVED: {matched_misfolded_outpath}')
        print(f'SAVED: {unmatched_misfolded_outpath}')

        ## save matched refolded df stats
        matched_refolded_outpath = os.path.join(self.outpath, 'Matched_refolded_stats.csv')
        matched_misfolded_outpath = os.path.join(self.outpath, 'Matched_misfolded_stats.csv')    
        unmatched_misfolded_outpath = os.path.join(self.outpath, 'UnMatched_misfolded_stats.csv')   
        self.matched_refolded.describe().to_csv(matched_refolded_outpath)
        self.matched_misfolded.describe().to_csv(matched_misfolded_outpath)
        self.unmatched_misfolded.describe().to_csv(unmatched_misfolded_outpath)
        print(f'SAVED: {matched_refolded_outpath}')
        print(f'SAVED: {matched_misfolded_outpath}')
        print(f'SAVED: {unmatched_misfolded_outpath}')
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
        df2_copy['diff'] = df2_copy['prot_size'] - size
        df2_copy['diff'] = np.abs(df2_copy['diff'])
        #print(df2_copy)
        minidx = df2_copy['diff'].idxmin()
        match = df2_copy.loc[minidx]
        #print(size, match)

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
    parser.add_argument("-f", "--res_features", type=str, required=True, help="path to residue feature files")
    args = parser.parse_args()

    # Setup logging configuration
    logging.basicConfig(filename=args.log, level=logging.INFO, format='%(asctime)s %(message)s')
    logging.info(f'{"#"*100}\nNEW RUN {script_name}')

    analysis = DataAnalysis(
        refolded=args.refolded,
        outpath=args.outpath,
        misfolded=args.misfolded,
        ent_features=args.ent_features,
        res_features=args.res_features)
    analysis.run()

if __name__ == "__main__":
    main()

