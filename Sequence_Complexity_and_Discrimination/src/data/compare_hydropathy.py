import sys, os, re, time, logging
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
from scipy.stats import permutation_test, ttest_ind
import matplotlib.pyplot as plt
from math import comb
#pd.set_option('display.max_rows', 5000)

################################################################################################################################
class HydropathyAnalyzer:
    """
    This class is meant to analyze the hydropathic contents of key contact involved in native entanglements
    (1) given a set of genes (essential, non_essential)
        (2) get all loop forming contacts for the unique entanglements in each dataset
            (3) for each class of hydropathy (strong hydrophobic, weak hydrophobic, and hydrophilic)
                (4) determine the fraction of loop forming native contacts that involve the class
    bootstrap set 2 to get confidence intervals

    """

    def __init__(self, outpath, contact_files, uent_files):
        """
        Initializing the FrequencyGenerator object and ensure atleast the top level outpath is present and makes it if not. 
        """

        # Make outpath if it doesnt exists
        self.outpath = outpath
        if not os.path.exists(self.outpath):
            os.makedirs(self.outpath)
            print(f'Made directory: {self.outpath}')

        self.three_to_one_letter = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'MSE': 'M', 'PHE': 'F', 
        'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 
        'VAL': 'V'}

        self.strong_hydrophobic = ['I', 'L', 'V', 'F', 'C', 'A']
        self.weak_hydrophobic = ['G', 'T', 'S', 'W', 'Y', 'P']
        self.total_hydrophobic = ['I', 'L', 'V', 'F', 'C', 'A', 'G', 'T', 'S', 'W', 'Y', 'P']
        self.hydrophilic = ['H', 'Q', 'E', 'N', 'D', 'K', 'R']
        self.aromatic = ['H', 'W', 'Y', 'F']

        # Load only those contact files in the gene_mask
        self.contact_files = [f for f in glob(os.path.join(contact_files, '*'))]
        print(f'Number of contact files: {len(self.contact_files)}')
    
        # Load unique ent files
        self.uent_files = [f for f in glob(os.path.join(uent_files, '*'))]
        print(f'Number of uent files: {len(self.uent_files)}')

    def hydropathy(self, gene_mask):

        # Load only those contact files in the gene_mask
        mask_contact_files = [f for f in self.contact_files if f.split('/')[-1].split('_')[0] in gene_mask]
        print(f'Number of contact files in mask: {len(mask_contact_files)}')
    
        # Load unique ent files
        mask_uent_files = [f for f in self.uent_files if f.split('/')[-1].split('_')[0] in gene_mask]
        print(f'Number of uent files in mask: {len(mask_uent_files)}')

        contact_df = {'SHphob':[], 'WHphob':[], 'Hphil':[], 'Arom':[], 'pureHphob':[], 'pureHphil':[], 'pureArom':[]}
        for fi, f in enumerate(mask_uent_files):
            gene = f.split('/')[-1].split('_')[0] 
            print(fi, f, gene)

            # check if gene is in the mask and if not move on
            if gene not in gene_mask:
                continue

            #uent_df = pd.read_csv(f, sep='|')['contacts'].values
            uent_df = pd.read_csv(f, sep='|')
            uent_df = uent_df[uent_df['CCBond'] == False]
            print(uent_df)
            uent_df = uent_df['contacts'].values
            uent_df = ';'.join(uent_df)
            print(uent_df)
            
            loop_contacts = [[x.split('-')[0], x.split('-')[1]] for x in uent_df.split(';') if len(x.split('-')) == 2]
            #print(loop_contacts)

            #get contact file for this gene
            cf = [cf for cf in self.contact_files if gene in cf][0]
            #print(cf)

            contacts = pd.read_csv(cf, sep='|')
            #print(contacts)

            for loop in loop_contacts:
                i = int(loop[0])
                j = int(loop[1])
                #print(loop, i, j)

                loop_df = contacts[contacts['pdb_resid_i'].isin([i, j])] 
                loop_df = loop_df[loop_df['pdb_resid_j'].isin([i, j])]
                #print(loop_df)

                # get the 3 letter resname for the loop forming contact residues
                # convert to 1 letter
                iAA = loop_df.iloc[0]['pdb_resname_i']
                jAA = loop_df.iloc[0]['pdb_resname_j']
                if iAA in self.three_to_one_letter:
                    iAA = self.three_to_one_letter[iAA]
                else:
                    #print(f'{iAA} is not a recognized amino acid')
                    continue

                if jAA in self.three_to_one_letter:
                    jAA = self.three_to_one_letter[jAA]
                else:
                    #print(f'{jAA} is not a recognized amino acid')
                    continue
                #print(iAA, jAA)

                # check if either residue in the loop closing contact is in the strong hydrophobic list
                # self.strong_hydrophobic = ['I', 'L', 'V', 'F', 'C', 'A']
                if iAA in self.strong_hydrophobic:
                    contact_df['SHphob'] += [1]
                else:
                    contact_df['SHphob'] += [0]
                
                # check if either residue in the loop closing contact is in the weak hydrophobic list
                #self.weak_hydrophobic = ['G', 'T', 'S', 'W', 'Y', 'P']
                if iAA in self.weak_hydrophobic:
                    contact_df['WHphob'] += [1]
                else:
                    contact_df['WHphob'] += [0]

                # check if either residue in the loop closing contact is in the hydrophilic list
                #self.hydrophilic = ['H', 'Q', 'E', 'N', 'D', 'K', 'R']
                if iAA in self.hydrophilic:
                    contact_df['Hphil'] += [1]
                else:
                    contact_df['Hphil'] += [0]
                    
                # check if either residue in the loop closing contact is in the aromatic list
                #self.aromatic = ['H', 'W', 'Y', 'F']
                if iAA in self.aromatic:
                    contact_df['Arom'] += [1]
                else:
                    contact_df['Arom'] += [0]
                
                # check if either residue in the loop closing contact is in the strong hydrophobic list
                # self.strong_hydrophobic = ['I', 'L', 'V', 'F', 'C', 'A']
                if jAA in self.strong_hydrophobic:
                    contact_df['SHphob'] += [1]
                else:
                    contact_df['SHphob'] += [0]
                
                # check if either residue in the loop closing contact is in the weak hydrophobic list
                #self.weak_hydrophobic = ['G', 'T', 'S', 'W', 'Y', 'P']
                if jAA in self.weak_hydrophobic:
                    contact_df['WHphob'] += [1]
                else:
                    contact_df['WHphob'] += [0]

                # check if either residue in the loop closing contact is in the hydrophilic list
                #self.hydrophilic = ['H', 'Q', 'E', 'N', 'D', 'K', 'R']
                if jAA in self.hydrophilic:
                    contact_df['Hphil'] += [1]
                else:
                    contact_df['Hphil'] += [0]
                    
                # check if either residue in the loop closing contact is in the aromatic list
                #self.aromatic = ['H', 'W', 'Y', 'F']
                if jAA in self.aromatic:
                    contact_df['Arom'] += [1]
                else:
                    contact_df['Arom'] += [0]

                # check if both residues are hydrophobic
                if iAA in self.total_hydrophobic and jAA in self.total_hydrophobic:
                    contact_df['pureHphob'] += [1, 1]
                else:
                    contact_df['pureHphob'] += [0, 0]

                # check if both residues are hydrophilic
                if iAA in self.hydrophilic and jAA in self.hydrophilic:
                    contact_df['pureHphil'] += [1, 1]
                else:
                    contact_df['pureHphil'] += [0, 0]

                # check if both residues are aromatic
                if iAA in self.aromatic and jAA in self.aromatic:
                    contact_df['pureArom'] += [1, 1]
                else:
                    contact_df['pureArom'] += [0, 0]

        return pd.DataFrame(contact_df)

    def boot_stats(self, df, tag):
        """
        Calculate means and confidence intervals for all columns in a DataFrame.
        
        Parameters:
            df (DataFrame): The DataFrame to process.
            
        Returns:
            DataFrame: Results including mean, confidence intervals, and class names.
        """
        results = []
        for column in df.columns:
            mean, lower, upper = self.bootstrap_ci(df[column])
            #results.append({f'class_{tag}': column, f'mean_{tag}': mean, f'lower_ci_{tag}': lower, f'upper_ci_{tag}': upper})
            results.append({'class': column, f'mean_{tag}': mean, f'lower_ci_{tag}': lower, f'upper_ci_{tag}': upper})
        
        return pd.DataFrame(results)


    def bootstrap_ci(self, data, n_bootstrap_samples=1000, ci=95):
        """
        Bootstrap the mean and confidence interval for given data.
        
        Parameters:
            data (array-like): The data to bootstrap.
            n_bootstrap_samples (int): Number of bootstrap samples.
            ci (float): Confidence interval (0-100).
            
        Returns:
            tuple: Mean, lower confidence interval, upper confidence interval.
        """
        bootstrap_samples = np.random.choice(data, size=(n_bootstrap_samples, len(data)), replace=True)
        sample_means = np.mean(bootstrap_samples, axis=1)
        mean = np.mean(sample_means)
        lower_ci = np.percentile(sample_means, (100 - ci) / 2)
        upper_ci = np.percentile(sample_means, 100 - (100 - ci) / 2)
        return mean, lower_ci, upper_ci


    ############################################################################
    def permutation_test(self, data1, data2, num_permutations=10000):
        """
        Perform a permutation test to calculate the p-value for the hypothesis
        that the samples in data1 and data2 have the same mean.

        Parameters:
            data1 (array-like): Data from the first condition.
            data2 (array-like): Data from the second condition.
            num_permutations (int): Number of permutations to perform.

        Returns:
            float: The p-value from the permutation test.
        """
        observed_diff = np.mean(data1) - np.mean(data2)
        combined_data = np.concatenate([data1, data2])

        more_extreme = 0
        for _ in range(num_permutations):
            np.random.shuffle(combined_data)
            new_data1 = combined_data[:len(data1)]
            new_data2 = combined_data[len(data1):]
            new_diff = np.mean(new_data1) - np.mean(new_data2)

            if abs(new_diff) >= abs(observed_diff):
                more_extreme += 1

        p_value = more_extreme / num_permutations
        return p_value

    def compare_columns(self, df1, df2):
        """
        Compare each column of df1 and df2 using a permutation test and return p-values.

        Parameters:
            df1 (DataFrame): First DataFrame containing the data.
            df2 (DataFrame): Second DataFrame containing the data.

        Returns:
            Series: P-values for the comparison of each column.
        """
        p_values = {'class':[], 'pvalue':[]}
        for column in df1.columns:
            pval = self.permutation_test(df1[column].values, df2[column].values)
            p_values['class'] += [column]
            p_values['pvalue'] += [pval]
            #p_values[column] = self.permutation_test(df1[column].values, df2[column].values)
        return pd.DataFrame(p_values)

#################################################################################################################


def main():
    """
    Main function to control workflow. 
    (1) parse user arguments 
    (2) making logging file
    (3) attempt to make F_ab, Fc_ab, FcG_ab(Ec1,...Ec4)
    """

    # Parse the user supplied arguments
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-eg", "--Ess_gene_list", type=str, required=True, help=f"path to Ess gene list")
    parser.add_argument("-neg", "--NonEss_gene_list", type=str, required=True, help=f"path to Ess gene list")
    parser.add_argument("-l", "--log_file", type=str, required=True, help="Path to logging file")
    parser.add_argument("-c", "--contact_files", type=str, required=True, help="path to native contact files")
    parser.add_argument("-e", "--uent_files", type=str, required=True, help="path to unique entanglement files")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="path to output directory. will be made if doesnt exist")
    parser.add_argument("-t", "--tag",  type=str, required=True, help="tag for output file")
    args = parser.parse_args()

    contact_files = args.contact_files
    uent_files = args.uent_files
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

    # Load gene_mask
    Ess_gene_mask = np.loadtxt(Ess_gene_list, dtype=str)
    print(f'Ess_gene_mask: {Ess_gene_mask} {len(Ess_gene_mask)}')

    NonEss_gene_mask = np.loadtxt(NonEss_gene_list, dtype=str)
    print(f'NonEss_gene_mask: {NonEss_gene_mask} {len(NonEss_gene_mask)}')


    # Get hydropathy estimates
    Hanal = HydropathyAnalyzer(outpath, contact_files, uent_files)

    Ess_hydropathy_df = Hanal.hydropathy(Ess_gene_mask)
    print(f'Ess_hydropathy_df:\n{Ess_hydropathy_df}')
    print(Ess_hydropathy_df.mean())

    NonEss_hydropathy_df = Hanal.hydropathy(NonEss_gene_mask)
    print(f'NonEss_hydropathy_df:\n{NonEss_hydropathy_df}')
    print(NonEss_hydropathy_df.mean())

    # Get boot statistics
    Ess_boot_stats_df = Hanal.boot_stats(Ess_hydropathy_df, 'Ess')
    print(f'Ess_boot_stats_df:\n{Ess_boot_stats_df}')

    NonEss_boot_stats_df = Hanal.boot_stats(NonEss_hydropathy_df, 'NonEss')
    print(f'NonEss_boot_stats_df:\n{NonEss_boot_stats_df}')

    # Get permutation
    pvals = Hanal.compare_columns(Ess_hydropathy_df, NonEss_hydropathy_df)
    print(f'pvals:\n{pvals}')

    # Make final df for saving
    df = pd.concat((Ess_boot_stats_df, NonEss_boot_stats_df.loc[:, NonEss_boot_stats_df.columns != 'class']), axis=1)
    df['pvalue'] = pvals['pvalue'].values
    print(f'df:\n{df}')
    df_outfile = f'{Hanal.outpath}{tag}_stats.csv'
    df.to_csv(df_outfile, sep='|', index=False)
    print(f'SAVED: {df_outfile}')

if __name__ == "__main__":
    main()

print('NORMAL TERMINATION')
