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

    def __init__(self, outpath, gene_list, ess_gene_list, ent_gene_list, client_gene_list, tag, buff, spa, LiPMScov):
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
        genes = np.loadtxt(gene_list, dtype=str)
        print(f'Loaded genes: {len(genes)}')

        ## load essential genes
        ess_genes = np.loadtxt(ess_gene_list, dtype=str)
        print(f'Loaded essential genes: {len(ess_genes)}')

        ## load essential ent genes
        ent_genes = np.loadtxt(ent_gene_list, dtype=str)
        print(f'Loaded ent genes: {len(ent_genes)}')


        ## load clients and only those that have some interaction with either DnaK or GroEL
        clients = pd.read_csv(client_gene_list, sep='|')
        print(clients)
        DnaK_clients = clients[clients['ChapClass'].isin(['D'])]
        GroEL_clients = clients[clients['ChapClass'].isin(['G'])]
        DnaK_clients = np.unique(DnaK_clients['gene'].values)
        GroEL_clients = np.unique(GroEL_clients['gene'].values)
        print(f'DnaK_clients: {len(DnaK_clients)}')
        print(f'GroEL_clients: {len(GroEL_clients)}')


        ## make initial dataframe for all genes 
        df = {'gene':[], 'essential':[], 'entangled':[], 'ess_and_ent':[],  'DnaK_client':[], 'GroEL_client':[], 'DnaKorGroEL_client':[]}
        for gene in genes:

            if gene in ess_genes:
                ess = True
            else:
                ess = False

            if gene in ent_genes:
                ent = True
            else:
                ent = False

            if ess and ent:
                ess_ent = True
            else:
                ess_ent = False

            if gene in DnaK_clients:
                DnaK_client = True
            else:
                DnaK_client = False

            if gene in GroEL_clients:
                GroEL_client = True
            else:
                GroEL_client = False

            if DnaK_client or GroEL_client:
                DnaKorGroEL_client = True
            else:
                DnaKorGroEL_client = False

            df['gene'] += [gene]
            df['essential'] += [ess]
            df['entangled'] += [ent]
            df['ess_and_ent'] += [ess_ent]
            df['DnaK_client'] += [DnaK_client]
            df['GroEL_client'] += [GroEL_client]
            df['DnaKorGroEL_client'] += [DnaKorGroEL_client]

        df = pd.DataFrame(df)
        print(df)
        self.df = df
        outfile = os.path.join(self.Outpath, f'Inpdata_{self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}.csv')
        df.to_csv(outfile, sep='|')
        print(f'SAVED: {outfile}')
        


    ##################################################################
    def Fisher(self, df, client_type):
        print(f'Calculating the OR and pvalue for the association of essential and a client of {client_type}')
        outfile = f'{self.Outpath}Fisher_stats_{self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}_clientType-{client_type}.csv'
        print(f'outfile: {outfile}')

        reps = 10000

        df_copy = df.copy()
        data = {'OR':[], 
                'OR_lb':[], 
                'OR_ub':[], 
                'pvalue':[], 
                'tag':[self.tag], 
                'buff':[self.buff], 
                'spa':[self.spa], 
                'LiPMScov':[self.LiPMScov], 
                'clientN_essentialN':[], 'clientY_essentialN':[], 'clientN_essentialY':[], 'clientY_essentialY':[], 
                'total_res':[], 'n':[]}

        df_copy = df_copy[df_copy['entangled'] == True]
        table = pd.crosstab(df_copy['essential'], df_copy[client_type])
        print(f'ctable:\n{table}')
        if table.values.shape != (2,2):
            print(f'ERROR: Tabel is not 2x2. Sample size maybe to small for conditions {self.tag}_{self.buff}_spa{self.spa}_LiPMScov{self.LiPMScov}')
            return
        clientN_essentialN = table.values[0,0]
        clientY_essentialN = table.values[0,1]
        clientN_essentialY = table.values[1,0]
        clientY_essentialY = table.values[1,1]
        data['n'] += [len(df_copy)]
        data['clientN_essentialN'] += [clientN_essentialN]
        data['clientY_essentialN'] += [clientY_essentialN]
        data['clientN_essentialY'] += [clientN_essentialY]
        data['clientY_essentialY'] += [clientY_essentialY]
        data['total_res'] += [clientN_essentialN + clientY_essentialN + clientN_essentialY + clientY_essentialY]

        OR, pvalue = fisher_exact(table)

        ## get the OR confidence intervals
        boot_OR = []
        for b in range(reps):
            boot_df = df_copy.sample(n=len(df_copy), replace=True)
            b_table = pd.crosstab(boot_df['essential'], boot_df[client_type])

            #print(b_table)
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
        logging.info(f'OR: {OR} ({OR_lb}, {OR_ub}) with pvalue {pvalue}')
        #data['boot_OR'] += [boot_OR]
        data['OR'] += [OR]
        data['OR_lb'] += [OR_lb]
        data['OR_ub'] += [OR_ub]
        data['pvalue'] += [pvalue]

        outdf = pd.DataFrame(data)
        print(outdf.to_string())
        logging.info(f'\n{outdf.to_string()}')

        outdf.to_csv(outfile, sep="|", index=False)
        print(f'SAVED: {outfile}')
        logging.info(f'SAVED: {outfile}')


#################################################################################################################

def statistic(x, y, axis):
    return np.mean(x, axis=axis) - np.mean(y, axis=axis)

def main():
    """
    Main function to control workflow.
    """

    # Parse the user supplied arguments
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-Ag", "--gene_list", type=str, required=True, help=f"path to all gene list to use")
    parser.add_argument("-Eg", "--ess_gene_list", type=str, required=True, help=f"path to essential gene list to use")
    parser.add_argument("-Entg", "--ent_gene_list", type=str, required=True, help=f"path to ent gene list to use")
    parser.add_argument("-c", "--client_gene_list", type=str, required=True, help=f"path to list of genes that are known clients")
    parser.add_argument("-l", "--log_file", type=str, required=True, help="Path to logging file")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="path to output directory. will be made if doesnt exist")
    parser.add_argument("-t", "--tag", type=str, required=True, help="tag for final output image")
    parser.add_argument("-b", "--buff", type=str, required=True, help="buffer used C CD CG")
    parser.add_argument("-s", "--spa", type=str, required=True, help="spa used 0, 10, 20, 30, 40, 50, 60, 70, 80, 90")
    parser.add_argument("--LiPMScov", type=str, required=True, help="LiPMS coverage used 0, 10, 20, 30, 40, 50, 60, 70, 80, 90")
    args = parser.parse_args()

    log_file = args.log_file
    gene_list = args.gene_list
    ess_gene_list = args.ess_gene_list
    ent_gene_list = args.ent_gene_list
    client_gene_list = args.client_gene_list
    outpath = args.outpath
    tag = args.tag
    buff = args.buff
    spa = args.spa
    LiPMScov = args.LiPMScov

    # Setup logging configuration
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s %(message)s')
    logging.info(f'{"#"*100}\nNEW RUN')

    # Initalize the FrequencyGenerator class object
    analysis = Analyzer(outpath, gene_list, ess_gene_list, ent_gene_list, client_gene_list, tag, buff, spa, LiPMScov)

    analysis.Fisher(analysis.df[['essential', 'entangled', 'DnaKorGroEL_client']], 'DnaKorGroEL_client')
    analysis.Fisher(analysis.df[['essential', 'entangled', 'DnaK_client']], 'DnaK_client')
    analysis.Fisher(analysis.df[['essential', 'entangled', 'GroEL_client']], 'GroEL_client')


if __name__ == "__main__":
    main()

print('NORMAL TERMINATION')
logging.info('NORMAL TERMINATION')
