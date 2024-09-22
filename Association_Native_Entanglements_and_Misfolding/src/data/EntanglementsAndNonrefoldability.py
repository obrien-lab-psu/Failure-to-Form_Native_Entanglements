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
    This class is meant to calculate the effect size of the signigicant differences obsersed in the compare_energetics_with_permutation_analysis

    """

    def __init__(self, outpath, all_gene_list, ent_gene_list, nonrefold_gene_list, tag, buff, spa, LiPMScov):
        """
        Initializing the FrequencyGenerator object and ensure atleast the top level outpath is present and makes it if not.
        """

        # Make outpath if it doesnt exists
        self.Outpath = outpath
        if not os.path.exists(self.Outpath):
            os.makedirs(self.Outpath)
            print(f'Made directory: {self.Outpath}')

        self.tag = tag
        self.buff = buff
        self.spa = spa
        self.LiPMScov = LiPMScov
        print(f'Tag: {self.tag} | buff: {buff} | spa: {spa} | LiPMScov: {LiPMScov}')
        

        ## load gene lists 
        all_genes = np.loadtxt(all_gene_list, dtype=str)
        print(f'Loaded all_genes: {len(all_genes)}')
        ent_genes = np.loadtxt(ent_gene_list, dtype=str)
        print(f'Loaded ent_genes: {len(ent_genes)}')
        nonrefold_genes = np.loadtxt(nonrefold_gene_list, dtype=str)
        print(f'Loaded nonrefold_genes: {len(nonrefold_genes)}')


        ## make initial dataframe for all genes 
        df = {'gene':[], 'entangled':[], 'nonrefolded':[]}
        for gene in all_genes:

            if gene in ent_genes:
                ent = True
            else:
                ent = False

            if gene in nonrefold_genes:
                nonrefolded = True
            else:
                nonrefolded = False

            df['gene'] += [gene]
            df['entangled'] += [ent]
            df['nonrefolded'] += [nonrefolded]

        df = pd.DataFrame(df)
        print(df)
        self.df = df
        


    ##################################################################
    def Fisher(self,):
        print(f'Calculating the OR and pvalue for the association of entangled genes with being nonrefoldable')
        outfile = f'{self.Outpath}Fisher_stats_{self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}.csv'
        print(f'outfile: {outfile}')

        reps = 10000
        data = {'OR':[], 'OR_lb':[], 'OR_ub':[], 'pvalue':[], 'tag':[self.tag], 'buff':[self.buff], 'spa':[self.spa], 'LiPMScov':[self.LiPMScov], 'entN_nonrefoldedN':[], 'entY_nonrefoldedN':[], 'entN_nonrefoldedY':[], 'entY_nonrefoldedY':[], 'total_res':[], 'n':[len(self.df)]}

        table = pd.crosstab(self.df['nonrefolded'], self.df['entangled'])
        print(f'ctable:\n{table}')
        if table.values.shape != (2,2):
            print(f'ERROR: Tabel is not 2x2. Sample size maybe to small for conditions {self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}_NLC{self.threshold}')
            quit()
        entN_nonrefoldedN = table.values[0,0]
        entY_nonrefoldedN = table.values[0,1]
        entN_nonrefoldedY = table.values[1,0]
        entY_nonrefoldedY = table.values[1,1]
        data['entN_nonrefoldedN'] += [entN_nonrefoldedN]
        data['entY_nonrefoldedN'] += [entY_nonrefoldedN]
        data['entN_nonrefoldedY'] += [entN_nonrefoldedY]
        data['entY_nonrefoldedY'] += [entY_nonrefoldedY]
        data['total_res'] += [entN_nonrefoldedN + entY_nonrefoldedN + entN_nonrefoldedY + entY_nonrefoldedY]

        OR, pvalue = fisher_exact(table)

        ## get the OR confidence intervals
        boot_OR = []
        for b in range(reps):
            boot_df = self.df.sample(n=len(self.df), replace=True)
            b_table = pd.crosstab(boot_df['nonrefolded'], boot_df['entangled'])
            if b_table.values.shape != (2,2):
                print(f'ERROR: Tabel is not 2x2. Sample size maybe to small for conditions {self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}')
                continue
            b_OR, b_pvalue = fisher_exact(b_table)
            boot_OR += [b_OR]
        boot_OR = np.asarray(boot_OR)
        finite_mask = np.isfinite(boot_OR)
        OR_lb = np.percentile(boot_OR[finite_mask], 2.5)
        OR_ub = np.percentile(boot_OR[finite_mask], 97.5)

        print(f'OR: {OR} ({OR_lb}, {OR_ub}) with pvalue {pvalue}')
        #data['boot_OR'] += [boot_OR]
        data['OR'] += [OR]
        data['OR_lb'] += [OR_lb]
        data['OR_ub'] += [OR_ub]
        data['pvalue'] += [pvalue]

        outdf = pd.DataFrame(data)
        print(outdf)

        outdf.to_csv(outfile, sep="|", index=False)
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
    parser.add_argument("-Ag", "--all_gene_list", type=str, required=True, help=f"path to all gene list to use")
    parser.add_argument("-Eg", "--ent_gene_list", type=str, required=True, help=f"path to ent gene list to use")
    parser.add_argument("-NRg", "--nonrefold_gene_list", type=str, required=True, help=f"path to nonrefolded gene list to use")
    parser.add_argument("-l", "--log_file", type=str, required=True, help="Path to logging file")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="path to output directory. will be made if doesnt exist")
    parser.add_argument("-t", "--tag", type=str, required=True, help="tag for final output image")
    parser.add_argument("-b", "--buff", type=str, required=True, help="buffer used C CD CG")
    parser.add_argument("-s", "--spa", type=str, required=True, help="spa used 0, 10, 20, 30, 40, 50, 60, 70, 80, 90")
    parser.add_argument("--LiPMScov", type=str, required=True, help="LiPMS coverage used 0, 10, 20, 30, 40, 50, 60, 70, 80, 90")
    args = parser.parse_args()

    log_file = args.log_file
    all_gene_list = args.all_gene_list
    ent_gene_list = args.ent_gene_list
    nonrefold_gene_list = args.nonrefold_gene_list
    outpath = args.outpath
    tag = args.tag
    buff = args.buff
    spa = args.spa
    LiPMScov = args.LiPMScov

    # Setup logging configuration
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s %(message)s')
    logging.info(f'{"#"*100}\nNEW RUN')

    # Initalize the FrequencyGenerator class object
    analysis = Analyzer(outpath, all_gene_list, ent_gene_list, nonrefold_gene_list, tag, buff, spa, LiPMScov)

    analysis.Fisher()


if __name__ == "__main__":
    main()

print('NORMAL TERMINATION')
logging.info('NORMAL TERMINATION')
