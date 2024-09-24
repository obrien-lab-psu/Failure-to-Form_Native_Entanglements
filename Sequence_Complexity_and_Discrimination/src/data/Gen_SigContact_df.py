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
from scipy.stats import permutation_test, ttest_ind, false_discovery_control, fisher_exact
import matplotlib.pyplot as plt
from math import comb
#pd.set_option('display.max_rows', 500)

class Analyzer:
    """
    This class is meant to calculate the effect size of the signigicant differences obsersed in the compare_energetics_with_permutation_analysis

    """

    def __init__(self, outpath, contact_files, uent_files):
        """
        Initializing the FrequencyGenerator object and ensure atleast the top level outpath is present and makes it if not. 
        """

        # Make outpath if it doesnt exists
        self.Outpath = outpath
        if not os.path.exists(self.Outpath):
            os.makedirs(self.Outpath)
            print(f'Made directory: {self.Outpath}')

        self.ContactsOutpath = os.path.join(outpath, 'Contacts/')
        if not os.path.exists(self.ContactsOutpath):
            os.makedirs(self.ContactsOutpath)
            print(f'Made directory: {self.ContactsOutpath}')

        # Define the indexing for amino acids in the 2D matrix to keep things consistent
        #self.AAs2idx = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19}
        self.AAs2idx = {'W': 0, 'F': 1, 'Y': 2, 'C': 3, 'I': 4, 'V': 5, 'L': 6, 'M': 7, 'H': 8, 'T': 9, 'R': 10, 'P': 11, 'Q': 12, 'N': 13, 'A': 14, 'S': 15, 'K': 16, 'G': 17, 'E': 18, 'D': 19}

        self.three_to_one_letter = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'MSE': 'M', 'PHE': 'F', 
        'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 
        'VAL': 'V'}
        
        self.contact_files = glob(os.path.join(contact_files, '*'))
        self.uent_files = glob(os.path.join(uent_files, '*'))
        print(f'number of entanglement files: {len(self.uent_files)}')
        
        genes = [f.split('/')[-1].split('_')[0] for f in self.uent_files]
        self.genes = genes
        print(f'genes: {len(genes)}')


    #################################################################################################################
    #################################################################################################################
    def get_contact_lib(self, tag):
        print(f'Getting contact library')
        logging.info(f'Getting contact library')
        contact_outfile = f'{self.ContactsOutpath}contacts_{tag}.csv'
        
        # Load only those contact files in the self.genes
        print(f'# genes in mask: {len(self.genes)}')
        contact_files = []
        for gene in self.genes:
            gene_contact_file = [f for f in self.contact_files if gene in f]
            #print(gene, gene_contact_file)
            if len(gene_contact_file) != 1:
                raise ValueError(f"The number of found contact files for gene {gene} != 1\n{gene_contact_file}")
            else:
                contact_files += [gene_contact_file[0]]
        print(f'# contact_files: {len(contact_files)}')

        ## QC that the number of genes in mask == number of contact files
        if len(self.genes) != len(contact_files):
            raise ValueError("The number of contact files found {len(contact_files)} != {len(self.genes)} ")

        contact_df = {'gene':[], 'contact_resids':[], 'contact_resnames':[], 'contact_resnames_single':[], 'contact_mapping2matrix':[]}
        if not os.path.exists(contact_outfile):
            for fi, f in enumerate(contact_files):
                gene = f.split('/')[-1].split('_')[0]

                contacts = pd.read_csv(f, sep='|')
                contacts = contacts[['pdb_resid_i', 'pdb_resname_i', 'pdb_resid_j', 'pdb_resname_j']].values
                pool = []
                for i, iAA, j, jAA in contacts:

                    if (i, j) in pool or (j,i) in pool:
                        continue

                    if iAA in self.three_to_one_letter and jAA in self.three_to_one_letter:
                        x = self.AAs2idx[self.three_to_one_letter[iAA]]
                        y = self.AAs2idx[self.three_to_one_letter[jAA]]
                        pool += [(i, j)]
                        #0 ../proteome_features/Gen_proteome_features/logRN_1_FDR-False_0.01/contact_type2_lib/P00452_5CNS_B_native_contacts.csv 4 ASN 51 ASP 13 19
                        contact_df['gene'] += [gene]
                        contact_df['contact_resids'] += [(i, j)]
                        contact_df['contact_resnames'] += [(iAA, jAA)]
                        contact_df['contact_resnames_single'] += [(self.three_to_one_letter[iAA], self.three_to_one_letter[jAA])]
                        contact_df['contact_mapping2matrix'] += [(x, y)]
                        #print(fi, f, i, iAA, j, jAA, x, y)
            

            contact_df = pd.DataFrame(contact_df)
            contact_df.to_csv(contact_outfile, sep='|', index=False)
            print(f'SAVED: {contact_outfile}')

        print(f'LOADING: {contact_outfile}')
        contact_df = pd.read_csv(contact_outfile, sep='|')

        ## QC that the number of genes in mask == number of contact files
        return contact_df
        
    #################################################################################################################
    def get_GE_contact_lib(self, contact_df, tag):
        print(f'Getting GE contact library')
        logging.info(f'Getting GE contact library')
        GE_contact_outfile = f'{self.ContactsOutpath}GE_contacts_{tag}.csv'

        if not os.path.exists(GE_contact_outfile):
            # Load unique ent files
            #print(f'Number of uent files: {len(self.uent_files)}')

            GE_contact_df = []
            # Check if the file already exists. If so tell user and skip.
            for fi, f in enumerate(self.uent_files):
                gene = f.split('/')[-1].split('_')[0]
                logging.info(f'{gene} {f}')
                gene_contacts = contact_df[contact_df['gene'] == gene]
                logging.info(gene_contacts)

                uent_df = pd.read_csv(f, sep='|')
                uent_df = uent_df[uent_df['CCBond'] == False]
                print(uent_df)
                quit()
                uent_df = uent_df['contacts'].values
                uent_df = ';'.join(uent_df)
                logging.info(uent_df)
                
                loop_contacts = [[x.split('-')[0], x.split('-')[1]] for x in uent_df.split(';') if len(x.split('-')) == 2]
                logging.info(loop_contacts)
                
                loop_contact_strs = []
                for loop_contact in loop_contacts:
                    loop_contact_str1 = f'({loop_contact[0]}, {loop_contact[1]})'
                    loop_contact_str2 = f'({loop_contact[1]}, {loop_contact[0]})'
                    loop_contact_strs += [loop_contact_str1, loop_contact_str2]
                    #logging.info(loop_contact, loop_contact_str1, loop_contact_str2)

                loop_contact_info = gene_contacts[gene_contacts['contact_resids'].isin(loop_contact_strs)]
                logging.info(loop_contact_info)
                GE_contact_df += [loop_contact_info]

            GE_contact_df = pd.concat(GE_contact_df)
            GE_contact_df.to_csv(GE_contact_outfile, sep='|', index=False)
            logging.info(f'SAVED: {GE_contact_outfile}')

        logging.info(f'LOADING: {GE_contact_outfile}')
        print(f'LOADING: {GE_contact_outfile}')
        GE_contact_df = pd.read_csv(GE_contact_outfile, sep='|')

        return GE_contact_df

    ##################################################################
    def get_protein_groups(self, GE_contacts):
        print(f'Getting protein groups')
        logging.info(f'Getting protein groups')
        """
        Get protein groups defined by whther they have atleast n number of loop closing contacts.
        where n is the threshold speicifed by the user.

        Group 1 proteins: atleast n hydrophobic hydrophobic loop forming contacts
        Group -1 proteins: no hydrophobic hydrophobic loop forming contacts

        Group 2 proteins: atleast n strong_hydrophobic strong_hydrophobic loop forming contacts
        Group -2 proteins: no strong_hydrophobic strong_hydrophobic loop forming contacts

        Group 3 proteins: atleast n aromatic aromatic loop forming contacts
        Group -3 proteins: no aromatic aromatic loop forming contacts

        Group 4 proteins: atleast n hydrophilic hydrophilic loop forming contacts
        Group -4 proteins: no hydrophilic hydrophilic loop forming contacts

        Group 5 proteins: atleast n hydrophobic charged loop forming contacts
        Group -5 proteins: no hydrophobic charged loop forming contacts

        Group 6 proteins: atleast n F-F, F-Y, T-Y loop forming contacts
        Group -6 proteins: no F-F, F-Y, T-Y loop forming contacts

        Group 7 proteins: atleast n M-R loop forming contacts
        Group -7 proteins: no M-R loop forming contacts

        Group 8 proteins: atleast n F-F, F-Y, T-Y, S-Y loop forming contacts
        Group -8 proteins: no F-F, F-Y, T-Y, S-Y loop forming contacts

        Group 9 proteins: atleast n M-R, R-V loop forming contacts
        Group -9 proteins: no M-R, R-V loop forming contacts
        """

        strong_hydrophobic = ['I', 'L', 'V', 'F', 'C', 'A']
        weak_hydrophobic = ['G', 'T', 'S', 'W', 'Y', 'P']
        total_hydrophobic = ['I', 'L', 'V', 'F', 'C', 'A', 'G', 'T', 'S', 'W', 'Y', 'P']
        hydrophilic = ['H', 'Q', 'E', 'N', 'D', 'K', 'R']
        aromatic = ['H', 'W', 'Y', 'F']

        outdf = {'gene':[]}
        sig_contacts = {}
        ## get group 1 significant contact defs (any hydrophobic hydrophobic contacts)
        sig_contacts['Hydrophob-Hydrophob'] = list(product(total_hydrophobic, repeat=2))
        sig_contacts['Hydrophob-Hydrophob'] = [str(s) for s in sig_contacts['Hydrophob-Hydrophob']]
        outdf['Hydrophob-Hydrophob'] = []

        ## get group 2 significant contact defs (strong hydrophobic strong hydrophobic contacts)
        sig_contacts['StrHydrophob-StrHydrophob'] = list(product(strong_hydrophobic, repeat=2))
        sig_contacts['StrHydrophob-StrHydrophob'] = [str(s) for s in sig_contacts['StrHydrophob-StrHydrophob']]
        outdf['StrHydrophob-StrHydrophob'] = []

        ## get group 3 significant contact defs (any aromatic aromatic contacts)
        sig_contacts['Aromatic-Aromatic'] = list(product(aromatic, repeat=2))
        sig_contacts['Aromatic-Aromatic'] = [str(s) for s in sig_contacts['Aromatic-Aromatic']]
        outdf['Aromatic-Aromatic'] = []

        ## get group 4 significant contact defs (any hydrophilic hydrophilic contacts)
        sig_contacts['Hydrophil-Hydrophil'] = list(product(hydrophilic, repeat=2))
        sig_contacts['Hydrophil-Hydrophil'] = [str(s) for s in sig_contacts['Hydrophil-Hydrophil']]
        outdf['Hydrophil-Hydrophil'] = []

        ## get group 5 significant contact defs (any aromatic hydrophobic contacts)
        sig_contacts['Hydrophob-Aromatic'] = list(product(total_hydrophobic, aromatic)) + list(product(aromatic, total_hydrophobic))
        sig_contacts['Hydrophob-Aromatic'] = [str(s) for s in sig_contacts['Hydrophob-Aromatic']]
        outdf['Hydrophob-Aromatic'] = []

        ## get group 5 significant contact defs (any charge hydrophobic contacts)
        sig_contacts['Hydrophob-Charged'] = list(product(total_hydrophobic, ['R', 'K', 'E', 'D', 'H'])) + list(product(['R', 'K', 'E', 'D', 'H'], total_hydrophobic))
        sig_contacts['Hydrophob-Charged'] = [str(s) for s in sig_contacts['Hydrophob-Charged']]
        outdf['Hydrophob-Charged'] = []

        ## get group 5 significant contact defs (any charge hydrophobic contacts)
        sig_contacts['Hydrophob-NegCharged'] = list(product(total_hydrophobic, ['E', 'D'])) + list(product(['E', 'D'], total_hydrophobic))
        sig_contacts['Hydrophob-NegCharged'] = [str(s) for s in sig_contacts['Hydrophob-NegCharged']]
        outdf['Hydrophob-NegCharged'] = []

        ## get group 5 significant contact defs (any charge hydrophobic contacts)
        sig_contacts['Hydrophob-PosCharged'] = list(product(total_hydrophobic, ['R', 'K', 'H'])) + list(product(['R', 'K', 'H'], total_hydrophobic))
        sig_contacts['Hydrophob-PosCharged'] = [str(s) for s in sig_contacts['Hydrophob-PosCharged']]
        outdf['Hydrophob-PosCharged'] = []

        ## get group of +charged +charged
        sig_contacts['PosCharged-PosCharged'] = list(product(['H', 'R', 'K'], repeat=2))
        sig_contacts['PosCharged-PosCharged'] = [str(s) for s in sig_contacts['PosCharged-PosCharged']]
        outdf['PosCharged-PosCharged'] = []

        ## get group of -charged -charged
        sig_contacts['NegCharged-NegCharged'] = list(product(['E', 'D'], repeat=2))
        sig_contacts['NegCharged-NegCharged'] = [str(s) for s in sig_contacts['NegCharged-NegCharged']]
        outdf['NegCharged-NegCharged'] = []

        ## get group 6 significant contact defs (any F-F, F-Y, T-Y, S-Y contacts)
        sig_contacts['FF-FY-SY-TY'] = [('F', 'F'), ('F', 'Y'), ('Y', 'F'), ('S', 'Y'), ('Y', 'S'), ('T', 'Y'), ('Y', 'T')]
        sig_contacts['FF-FY-SY-TY'] = [str(s) for s in sig_contacts['FF-FY-SY-TY']]
        outdf['FF-FY-SY-TY'] = []

        ## get group 6 significant contact defs (any F-F contacts)
        sig_contacts['FF'] = [('F', 'F')]
        sig_contacts['FF'] = [str(s) for s in sig_contacts['FF']]
        outdf['FF'] = []

        ## get group 6 significant contact defs (any F-Y contacts)
        sig_contacts['FY'] = [('F', 'Y'), ('Y', 'F')]
        sig_contacts['FY'] = [str(s) for s in sig_contacts['FY']]
        outdf['FY'] = []

        ## get group 6 significant contact defs (any T-Y contacts)
        sig_contacts['TY'] = [('T', 'Y'), ('Y', 'T')]
        sig_contacts['TY'] = [str(s) for s in sig_contacts['TY']]
        outdf['TY'] = []

        ## get group 6 significant contact defs (any S-Y contacts)
        sig_contacts['SY'] = [('S', 'Y'), ('Y', 'S')]
        sig_contacts['SY'] = [str(s) for s in sig_contacts['SY']]
        outdf['SY'] = []

        ## get group 7 significant contact defs (any charge hydrophobic contacts)
        sig_contacts['MR'] = [('M', 'R'), ('R', 'M')]
        sig_contacts['MR'] = [str(s) for s in sig_contacts['MR']]
        outdf['MR'] = []

        outdf['Total'] = []
        ### Loop through genes and determine if that gene belongs in any of the buckets
        num_genes_w_ent = len(GE_contacts['gene'].unique())
        print(f'num_genes_w_ent: {num_genes_w_ent}')
        for gene in self.genes:
            print(f'{"#"*50}\nAnalyzing {gene}')
            gene_GE_contacts = GE_contacts[GE_contacts['gene'] == gene]
            print(gene_GE_contacts)
            outdf['gene'] += [gene]
            outdf['Total'] += [len(gene_GE_contacts)]
            for group_label, group in sig_contacts.items():
                #group = [str(s) for s in group]
                #print(group_label,  group)

                gene_GE_group_contacts = gene_GE_contacts[gene_GE_contacts['contact_resnames_single'].isin(group)]
                #print(gene_GE_group_contacts)

                num_gene_GE_group_contacts = len(gene_GE_group_contacts)
                print(f'{group_label} num_gene_GE_group_contacts: {num_gene_GE_group_contacts}')
                outdf[group_label] += [num_gene_GE_group_contacts]


        outdf = pd.DataFrame(outdf)
        print(outdf)
        outfile = os.path.join(self.Outpath, f'LoopformingContactsClassification.csv')
        outdf.to_csv(outfile, sep='|', index=False)
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
    parser.add_argument("-c", "--contact_files", type=str, required=True, help="path to native contact files")
    parser.add_argument("-e", "--uent_files", type=str, required=True, help="path to unique entanglement files")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="path to output directory. will be made if doesnt exist")
    args = parser.parse_args()

    contact_files = args.contact_files
    uent_files = args.uent_files
    log_file = args.log_file
    outpath = args.outpath

    # Setup logging configuration
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s %(message)s') 
    logging.info(f'{"#"*100}\nNEW RUN')

    # Make outpath if it doesnt exists
    if not os.path.exists(outpath):
        os.makedirs(outpath)
        print(f'Made directory: {outpath}')

    # Initalize the FrequencyGenerator class object
    analysis = Analyzer(outpath, contact_files, uent_files)

    # Get
    All_contacts = analysis.get_contact_lib('All')
    print(f'All_contacts:\n{All_contacts}')

    # Get loop closing contacts
    All_GE_contacts = analysis.get_GE_contact_lib(All_contacts, 'All')
    print(f'All_GE_contacts:\n{All_GE_contacts}')

    # Get gene groups 
    analysis.get_protein_groups(All_GE_contacts)


if __name__ == "__main__":
    main()

print('NORMAL TERMINATION')
logging.info('NORMAL TERMINATION')
