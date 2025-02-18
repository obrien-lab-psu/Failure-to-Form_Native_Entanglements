import time, sys
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
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import poisson, binom, fisher_exact, chi2, norm
import scipy.stats as st
from matplotlib.ticker import MultipleLocator

#pd.set_option('display.max_rows', 4000)

class DataAnalysis:
    """
    A class to handle the data analysis process including encoding, regression, and statistical tests.
    """

    def __init__(self, outpath):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
        - outpath (str): Path to the output directory.
        """
        self.outpath = outpath
        self.data = {}

        if not os.path.exists(f'{self.outpath}'):
            os.makedirs(f'{self.outpath}')
            print(f'Made output directories {self.outpath}')

    def setup_logging(self):
        """
        Sets up the logging configuration.

        Returns:
        - logger (logging.Logger): Configured logger.
        """
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)
        return logger


    def load_data(self, dataset):
        """
        Loads the residue feature files and filters the data for analysis.
        """
        res_files = glob.glob(f'../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gen_proteome_features_{dataset}/res_features_lib/*.csv')
        print(f'Number of {dataset} res_files: {len(res_files)}')
        for i, gene_resFeat_file in enumerate(res_files):
            #gene_resFeat = [f for f in res_files if gene in f]

            #gene_resFeat_file = gene_resFeat[0]
            #print(f'gene_resFeat_file: {gene_resFeat_file} {i}')
            if len(self.data) == 0:
                self.data = pd.read_csv(gene_resFeat_file, sep='|')
            else:
                self.data = pd.concat((self.data, pd.read_csv(gene_resFeat_file, sep='|')))

        #self.data = self.data[self.data['gene'].isin(self.reg_genes)]
        self.data = self.data[self.data['AA'] != 'NC']
        self.data = self.data[self.data['mapped_resid'].notna()]
        self.data = self.data[self.data['AA'].notna()]
        self.data = self.data.reset_index()
        print(f"Data loaded and filtered. Number of unique genes: {len(self.data['gene'].unique())}")

#########################################################################################################
    def Population_anal(self):
        """

        """
        start_time = time.time()

        ### Calculate the fraction of proteins with atleast 1 cut in the entangled region

        tags = ['all_genes',
                'essential_genes',
                'nonessential_genes']

        stats_df = {'buff':[], 'gene':[], 'NativeEnt':[], 'NonRefoldable':[]} 

        for dataset in ['EXP', 'AF']:
            print(f'{"#"*100}\n{dataset}\n{"#"*100}')
            outfile_csv = os.path.join(self.outpath, f'NativeEntanglements_and_SigCuts_{dataset}.csv')
            if not os.path.exists(outfile_csv):
                self.load_data(dataset)

                for tag in tags:

                    for buff in ['C', 'CD', 'CG']:
                        loc_df = self.data[['gene', 'region', f'cut_{buff}_Rall']]

                        for spa in [50]:
                                                    
                            gene_list = f'../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/{dataset}/{dataset}_0.6g_{buff}_Rall_spa{spa}_LiPMScov50_{tag}.txt'
                            genes = np.loadtxt(gene_list, dtype=str)
                            #print(f'gene_list: {gene_list} {len(genes)}')
                            # get the binary array of genes that have cuts in their entangled regions

                            for gene, gene_df in loc_df.groupby('gene'):
                                if gene in genes:
                                    
                                    entR = gene_df[gene_df['region'] == 1]  
                                    entR_n = len(entR)  
                                    entR_cut = np.sum(entR[f'cut_{buff}_Rall'].values)

                                    nonentR = gene_df[gene_df['region'] == 0]    
                                    nonentR_n = len(nonentR)  
                                    nonentR_cut = np.sum(nonentR[f'cut_{buff}_Rall'].values)

                                    if entR_n != 0:
                                        NativeEnt = True
                                    else:
                                        NativeEnt = False
                                    
                                    if entR_cut > 0 or nonentR_cut > 0:
                                        NonRefoldable = True
                                    else:
                                        NonRefoldable = False

                                    stats_df['buff'] += [buff]
                                    stats_df['gene'] += [gene]
                                    stats_df['NativeEnt'] += [NativeEnt]                            
                                    stats_df['NonRefoldable'] += [NonRefoldable]
        
                stats_df = pd.DataFrame(stats_df)
                #print(stats_df)
                stats_df.to_csv(outfile_csv, index=False)
                print(f'SAVED: {outfile_csv}')
            else:
                stats_df = pd.read_csv(outfile_csv)
                print(f'LOADED: {outfile_csv}')
                #print(stats_df)


            ## Quick contingency table analysis 
            for buff, buff_df in stats_df.groupby('buff'):
                table, OR, pvalue = calculate_fisher_exact(buff_df)
                print(f'{"#"*100}\n{dataset} Buffer: {buff}\n{table}\nOdds Ratio = {OR} with p-value = {pvalue}')


#######################################################################################
#######################################################################################
def calculate_fisher_exact(df):
    """
    Calculates the contingency table, odds ratio, and p-value using Fisher's exact test.

    Args:
        df (pd.DataFrame): A dataframe with two columns: 'NativeEnt' and 'NonRefoldable', each containing True or False.

    Returns:
        tuple: Contingency table (2x2 numpy array), odds ratio (float), and p-value (float).
    """
    # Creating the contingency table
    contingency_table = pd.crosstab(df['NativeEnt'], df['NonRefoldable'])
    
    # Extracting the table as a 2x2 numpy array
    table = contingency_table.values
    print(table)
    
    # Calculating odds ratio and p-value using Fisher's exact test
    odds_ratio, p_value = fisher_exact(table)
    
    return contingency_table, odds_ratio, p_value
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

def main():
    """
    Main function to parse arguments and run the DataAnalysis class.
    """
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory")


    args = parser.parse_args()

    analysis = DataAnalysis(outpath=args.outpath)
    analysis.Population_anal()


if __name__ == "__main__":
    main()

