import sys, os, re, time, logging
from scipy.stats import bootstrap
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
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import permutation_test, ttest_ind, false_discovery_control, fisher_exact, mannwhitneyu
import matplotlib.pyplot as plt
from math import comb
#pd.set_option('display.max_rows', 500)

class Analyzer:
    """
    This class is meant to calculate the Fraction of hydrophobic residues in a set of the protein feature files. Cannot use just PDBs as we need to know which residues actually mapped

    """

    def __init__(self, outpath, feature_files, Ess_gene_list, NonEss_gene_list, tag):
        """
        Initializing the FrequencyGenerator object and ensure atleast the top level outpath is present and makes it if not. 
        """

        # Make outpath if it doesnt exists
        self.Outpath = outpath
        if not os.path.exists(self.Outpath):
            os.makedirs(self.Outpath)
            print(f'Made directory: {self.Outpath}')

        # Make DATA outpath if it doesnt exists
        self.Data_Outpath = os.path.join(self.Outpath, 'DATA/')
        if not os.path.exists(self.Data_Outpath):
            os.makedirs(self.Data_Outpath)
            print(f'Made directory: {self.Data_Outpath}')

        # Make PLOT outpath if it doesnt exists
        self.Plot_Outpath = os.path.join(self.Outpath, 'PLOT/')
        if not os.path.exists(self.Plot_Outpath):
            os.makedirs(self.Plot_Outpath)
            print(f'Made directory: {self.Plot_Outpath}')

        # Define the indexing for amino acids in the 2D matrix to keep things consistent
        #self.AAs2idx = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19}
        self.AAs2idx = {'W': 0, 'F': 1, 'Y': 2, 'C': 3, 'I': 4, 'V': 5, 'L': 6, 'M': 7, 'H': 8, 'T': 9, 'R': 10, 'P': 11, 'Q': 12, 'N': 13, 'A': 14, 'S': 15, 'K': 16, 'G': 17, 'E': 18, 'D': 19}

        self.three_to_one_letter = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'MSE': 'M', 'PHE': 'F', 
        'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 
        'VAL': 'V'}
        
        self.feature_files = glob(os.path.join(feature_files, '*'))
        print(f'number of feature files: {len(self.feature_files)}')
        
        self.Ess_gene_mask = np.loadtxt(Ess_gene_list, dtype=str)
        self.NonEss_gene_mask = np.loadtxt(NonEss_gene_list, dtype=str)
        print(f'Ess_gene_mask: {len(self.Ess_gene_mask)}')
        print(f'NonEss_gene_mask: {len(self.NonEss_gene_mask)}')

        self.tag = tag
    #################################################################################################################
 
    #################################################################################################################
    def Get_FractionHydrophobic(self, AAtype='strong_hydrophobic'):
        print(f'Getting Fraction of Hydrophobic residues for each protein in Ess and NonEss groups')
        logging.info(f'Getting Fraction of Hydrophobic residues for each protein in Ess and NonEss groups')
        """
        Getting Fraction of Hydrophobic residues for each protein in Ess and NonEss groups
        """

        AAtype_df = {
        'Strong_hydrophobic':['I', 'L', 'V', 'F', 'C', 'A'],
        'Weak_hydrophobic':['G', 'T', 'S', 'W', 'Y', 'P'],
        'Total_hydrophobic':['I', 'L', 'V', 'F', 'C', 'A', 'G', 'T', 'S', 'W', 'Y', 'P'],
        'General_hydrophobic':['A', 'V', 'L', 'I', 'M', 'F', 'W', 'P', 'Y'],
        'Hydrophilic':['H', 'Q', 'E', 'N', 'D', 'K', 'R'],
        'Aromatic':['H', 'W', 'Y', 'F']}

        outdf = {'gene':[], f'Frac{AAtype}':[], 'Essential':[]}
        outfile_csv = os.path.join(self.Plot_Outpath, f'{self.tag}_perProt_Fraction_{AAtype}.csv')
        stats_outfile_csv = os.path.join(self.Plot_Outpath, f'{self.tag}_perProt_Fraction_{AAtype}_STATS.csv')
        outfile = os.path.join(self.Plot_Outpath, f'{self.tag}_perProt_Fraction_{AAtype}.png')
        for essentialTAG, gene_list in [('Essential', self.Ess_gene_mask), ('NonEssential', self.NonEss_gene_mask)]:            
            ## load each feature file in the gene_list
            for gene in gene_list:
                for f in self.feature_files:
                    if gene in f:
                        #print(gene, f)

                        df = pd.read_csv(f, sep='|')
                        df  = df[~df['mapped_resid'].isnull()]
                        
                        #Hydrophobic_df = df[df['AA'].isin(total_hydrophobic)]
                        Hydrophobic_df = df[df['AA'].isin(AAtype_df[AAtype])]
                        num_Hydrophobic = len(Hydrophobic_df)
                        num_NonHydrophobic = len(df) - num_Hydrophobic
                        Frac = num_Hydrophobic / len(df)
                        #print(f'num_Hydrophobic: {num_Hydrophobic}')
                        #print(f'num_NonHydrophobic: {num_NonHydrophobic}')
                        #print(f'Frac: {Frac}')

                        outdf['gene'] += [gene]
                        outdf[f'Frac{AAtype}'] += [Frac]
                        outdf['Essential'] += [essentialTAG]

        outdf = pd.DataFrame(outdf)   
        print(f'outdf:\n{outdf}')
        outdf.to_csv(outfile_csv, sep='|', index=False)
        print(f'SAVED: {outfile_csv}')

        ## Plot the data
        analyze_and_plot(outdf, outfile, stats_outfile_csv, self.tag, f'Frac{AAtype}')


# Function to calculate confidence intervals using bootstrapping
def bootstrap_ci(data, n_resamples=10000, ci=95):
    bootstrapped_means = [
        np.mean(np.random.choice(data, size=len(data), replace=True))
        for _ in range(n_resamples)
    ]
    lower_bound = np.percentile(bootstrapped_means, (100 - ci) / 2)
    upper_bound = np.percentile(bootstrapped_means, 100 - (100 - ci) / 2)
    return np.mean(data), lower_bound, upper_bound

# Function to perform permutation test
#def permutation_test(data1, data2, n_permutations=10000):
#    observed_stat = np.abs(np.mean(data1) - np.mean(data2))
#    combined_data = np.concatenate([data1, data2])
#    permuted_stats = []
#    for _ in range(n_permutations):
#        np.random.shuffle(combined_data)
#        perm_data1 = combined_data[: len(data1)]
#        perm_data2 = combined_data[len(data1):]
#        permuted_stats.append(np.abs(np.mean(perm_data1) - np.mean(perm_data2)))
#    p_value = np.mean(np.array(permuted_stats) >= observed_stat)
#    return observed_stat, p_value
def statistic(x, y, axis):
    return np.mean(x, axis=axis) - np.mean(y, axis=axis)

# Function to analyze the data and generate the plot
def analyze_and_plot(dataframe, outfile, stats_outfile, tag, varKey):
    groups = dataframe.groupby('Essential')
    stats = {}
    data_ess = []
    stats_df = {'Essential':[], 'Mean':[], 'lb':[], 'ub':[], 'Eff-Size':[], 'pvalue':[]}
    for name, group in groups:
        mean, ci_lower, ci_upper = bootstrap_ci(group[varKey])
        stats[name] = {
            'mean': mean,
            'ci_lower': ci_lower,
            'ci_upper': ci_upper
        }
        data_ess.append((name, group[varKey].values))
        stats_df['Essential'] += [name]
        stats_df['Mean'] += [mean]
        stats_df['lb'] += [ci_lower]
        stats_df['ub'] += [ci_upper]
    
    # Perform permutation test
    (name1, data1), (name2, data2) = data_ess
    res = permutation_test((data1, data2), statistic)
    test_stat, p_value = res.statistic, res.pvalue
    #test_stat, p_value = mannwhitneyu(data1, data2)
    #test_stat, p_value = ttest_ind(data1, data2)

    ## save the stats file
    stats_df['Eff-Size'] += [test_stat]*len(groups)
    stats_df['pvalue'] += [p_value]*len(groups)
    stats_df = pd.DataFrame(stats_df)
    print(f'stats_df:\n{stats_df}')
    stats_df.to_csv(stats_outfile)
    print(f'SAVED: {stats_outfile}')

    # Plotting
    plt.figure(figsize=(10, 6))
    bins=30
    x1, x2 = np.linspace(0, 1, 1000), np.linspace(0, 1, 1000)
    plt.hist(data1, bins=bins, range=(min([min(data1), min(data2)]), max([max(data1), max(data2)])), alpha=0.5, label=f"{name1} (Mean: {stats[name1]['mean']:.3f}, CI: [{stats[name1]['ci_lower']:.3f}, {stats[name1]['ci_upper']:.3f}])", density=True)
    plt.hist(data2, bins=bins, range=(min([min(data1), min(data2)]), max([max(data1), max(data2)])), alpha=0.5, label=f"{name2} (Mean: {stats[name2]['mean']:.3f}, CI: [{stats[name2]['ci_lower']:.3f}, {stats[name2]['ci_upper']:.3f}])", density=True)
    # Adjust legend to the upper right corner
    plt.legend(loc='upper right', bbox_to_anchor=(1.0, 1.0))
    plt.title(f"PDF of {varKey} for Essential and NonEssential Genes\n{tag}")
    plt.xlabel(f"{varKey}")
    plt.ylabel("Density")
    
    # Place annotation in the upper left corner
    annotation = f"Permutation Test\nStatistic: {test_stat:.3f}\nP-value: {p_value:.3e}"
    plt.text(0.02, 0.98, annotation, fontsize=12, transform=plt.gca().transAxes,
             bbox=dict(facecolor='white', alpha=0.5), verticalalignment='top')
    
    plt.savefig(outfile)
    print(f'SAVED: {outfile}')



#################################################################################################################

def statistic(x, y, axis):
    return np.mean(x, axis=axis) - np.mean(y, axis=axis)

def main():
    """
    Main function to control workflow. 
    (1) parse user arguments 
    (2) making logging file
    (3) attempt to make F_ab, Fc_ab, FcG_ab(Ec1,...Ec4)
    """

    # Parse the user supplied arguments
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-l", "--log_file", type=str, required=True, help="Path to logging file")
    parser.add_argument("-f", "--feature_files", type=str, required=True, help="path to native contact files")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="path to output directory. will be made if doesnt exist")
    parser.add_argument("-Es", "--Ess_gene_list", type=str, required=True, help="Path to gene list for essential proteins")
    parser.add_argument("-NEs", "--NonEss_gene_list", type=str, required=True, help="Path to gene list for Nonessential proteins")
    parser.add_argument("-t", "--tag", type=str, required=True, help="tag for output files")

    args = parser.parse_args()

    feature_files = args.feature_files

    log_file = args.log_file
    outpath = args.outpath
    Ess_gene_list = args.Ess_gene_list
    NonEss_gene_list = args.NonEss_gene_list
    tag = args.tag

    # Setup logging configuration
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s %(message)s') 
    logging.info(f'{"#"*100}\nNEW RUN')

    # Make outpath if it doesnt exists
    if not os.path.exists(outpath):
        os.makedirs(outpath)
        print(f'Made directory: {outpath}')

    # Initalize the FrequencyGenerator class object
    analysis = Analyzer(outpath, feature_files, Ess_gene_list, NonEss_gene_list, tag)

    # Get fraction AAtypes in gene lists and calc some stats 
    analysis.Get_FractionHydrophobic(AAtype='Strong_hydrophobic')
    analysis.Get_FractionHydrophobic(AAtype='Weak_hydrophobic')
    analysis.Get_FractionHydrophobic(AAtype='Total_hydrophobic')
    analysis.Get_FractionHydrophobic(AAtype='General_hydrophobic')
    analysis.Get_FractionHydrophobic(AAtype='Hydrophilic')
    analysis.Get_FractionHydrophobic(AAtype='Aromatic')


if __name__ == "__main__":
    main()

print('NORMAL TERMINATION')
logging.info('NORMAL TERMINATION')
