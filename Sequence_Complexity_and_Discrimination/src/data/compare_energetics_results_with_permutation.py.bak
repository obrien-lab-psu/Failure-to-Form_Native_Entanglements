import sys, os, re, time, logging
from functools import reduce
from sklearn.utils import shuffle
import ast
from collections import defaultdict
import warnings
import multiprocessing as mp 
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
from scipy.stats import permutation_test, ttest_ind, false_discovery_control
import matplotlib.pyplot as plt
from math import comb
pd.set_option('display.max_rows', 500)
np.set_printoptions(linewidth=500)
# Suppress the specific RuntimeWarning
warnings.filterwarnings("ignore", message="divide by zero encountered in log10")

class FrequencyGenerator:
    """
    This class is meant to calculate the three frequencies in Travatos paper.
    (1) F_ab = N(a,b)/N 
        this is the frequency of a given pair of amino acids (a) and (b) being in present across the proteome
        where N(a,b) is the number of pairs of amino acids (a) and (b), and N is the total pairs of all amino acids across the proteome

    (2) Fc_ab = Nc(a,b)/Nc 
        this is the frequence of a given pair of amino acids (a) and (b) being in contact across the proteome
        where Nc(a,b) is the number of time (a) and (b) where in contact across the proteome, and Nc is the total number of contacts across the proteome

    (3) FcG_ab = NcG(a,b)/NcG
        thi is the frequence of a given pair of amino acids (a) and (b) being in contact and being an entangled native contact
        where NcG(a,b) is the number of entangled native contacts involving (a) and (b), and NcG is the total number of entangled contacts

        There are several ways to define entangled contacts and we will allow the following definitions. 
        Ec1 = loop closing native contacts (i.e. the loop closing native contacts in each cluster of unique entanglements)
        Ec2 = contacts between crossing and the rest of the protein (not inclusing bufer)
        Ec3 = all contacts within the entangled region (i.e. a,b must be in the set of entangled residues)
        Ec4 = all contacts involving the entangled region (i.e. atleast a or b must be in the set of entangled residues but doesnt require both like Ec3 does)
    """

    def __init__(self, outpath, fasta_files, contact_files, uent_files, floor):
        """
        Initializing the FrequencyGenerator object and ensure atleast the top level outpath is present and makes it if not. 
        """

        # Make outpath if it doesnt exists
        self.outpath = outpath
        self.FrequencyGeneratorOutpath = os.path.join(self.outpath, 'FrequencyGeneratorOutput/')
        if not os.path.exists(self.FrequencyGeneratorOutpath):
            os.makedirs(self.FrequencyGeneratorOutpath)
            print(f'Made directory: {self.FrequencyGeneratorOutpath}')
        self.FrequencyGeneratorContactsOutpath = os.path.join(self.outpath, 'FrequencyGeneratorOutput/Contacts/')
        if not os.path.exists(self.FrequencyGeneratorContactsOutpath):
            os.makedirs(self.FrequencyGeneratorContactsOutpath)
            print(f'Made directory: {self.FrequencyGeneratorContactsOutpath}')

        # Define the indexing for amino acids in the 2D matrix to keep things consistent
        #self.AAs2idx = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19}
        self.AAs2idx = {'W': 0, 'F': 1, 'Y': 2, 'C': 3, 'I': 4, 'V': 5, 'L': 6, 'M': 7, 'H': 8, 'T': 9, 'R': 10, 'P': 11, 'Q': 12, 'N': 13, 'A': 14, 'S': 15, 'K': 16, 'G': 17, 'E': 18, 'D': 19}

        self.three_to_one_letter = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'MSE': 'M', 'PHE': 'F', 
        'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 
        'VAL': 'V'}
        
        self.fasta_files = fasta_files
        self.contact_files = contact_files
        self.uent_files = uent_files
        self.floor = floor

    #################################################################################################################
    def F_ab(self, fasta_files, gene_mask):
        """
        Calculate 
            F_ab = N(a,b)/N 
            this is the frequency of a given pair of amino acids (a) and (b) being in present across the proteome
            where N(a,b) is the number of pairs of amino acids (a) and (b), and N is the total pairs of all amino acids across the proteome
        """ 

        # Load only those FASTA files in the gene_mask
        #self.fasta_files = [f for f in glob(fasta_files) if f.split('/')[-1].replace('.fasta','') in gene_mask]
        self.fasta_files = [f for f in glob(fasta_files) if f.split('/')[-1].split('_')[0].replace('.fasta','') in gene_mask]
        self.proteome_seq = ''
        for fi, f in enumerate(self.fasta_files):
            gene_seq = ''.join([l.strip('\n') for i,l in enumerate(open(f, 'r').readlines()) if i != 0])
            #print('\n', fi, f, gene_seq, len(gene_seq))
            self.proteome_seq += gene_seq
        print(f'proteome_seq length: {len(self.proteome_seq)}')
        logger.info(f'proteome_seq length: {len(self.proteome_seq)}')
        self.proteome_size = len(self.proteome_seq)

        # Initialize the F_ab outpath
        self.F_ab_outfile = f'{self.FrequencyGeneratorOutpath}F_ab.csv'

        # Check if the file already exists. If so tell user and skip.
        # Else make all pairs
        if not os.path.exists(self.F_ab_outfile):
            print(f'Generating F_ab... This may take some time')
            self.N_ab = self.count_pairs(self.proteome_seq, self.AAs2idx)
            #print(self.N_ab)
            N = np.triu(self.N_ab)
            N = np.sum(N)

            # Calculate N = number of total pairs of AA across the sequence
            self.N = math.comb(len(self.proteome_seq), 2)
            #self.N = N
            #print(f'N: {self.N}')

            # Calculate F_ab = N_ab/N and make it into a pandas df
            self.F_ab = self.N_ab/self.N
            AA_categories = self.AAs2idx.keys()
            self.F_ab_df = pd.DataFrame(self.F_ab, columns = AA_categories)
            self.F_ab_df['AA'] = AA_categories
            #print(f'self.F_ab_df:\n{self.F_ab_df}')

            self.F_ab_df.to_csv(self.F_ab_outfile, sep='|', index=False)
            print(f'SAVED: {self.F_ab_outfile}')
            logger.info(f'SAVED: {self.F_ab_outfile}')

        else:
            print(f'WARNING: {self.F_ab_outfile} Already exists. Delete if you wish to remake the file.')
            logger.info(f'WARNING: {self.F_ab_outfile} Already exists. Delete if you wish to remake the file.')
            print(f'Loading: {self.F_ab_outfile}')
            logger.info(f'Loading: {self.F_ab_outfile}')
            self.F_ab_df = pd.read_csv(self.F_ab_outfile, sep='|')
            #print(f'self.F_ab_df:\n{self.F_ab_df}')
    #################################################################################################################

    #################################################################################################################
    def count_pairs(self, s, char_index):
        """
        To efficiently estimate the number of pairs of amino acids across the proteome use the combenitorics shortcut
        (1) count the frequence of each character
        (2) compute thenumber nuber of pairs 
            For each pair of characters (i, j), the number of such pairs can be calculated by multiplying the frequency of i by the frequency of j. 
            If i is the same as j, you need to calculate the number of ways to pick 2 items from the frequency of i, which is given by the formula for combinations: Comb(n,2)=(n(n-1))/2.
        """
        
        # Initialize frequency dictionary
        freq = defaultdict(int)
        
        # Count frequencies of each character
        for char in s:
            freq[char] += 1
        
        # Initialize the matrix
        matrix = np.zeros((20, 20), dtype=int)
        
        # Compute number of pairs
        for char_i in char_index:
            idx_i = char_index[char_i]
            for char_j in char_index:
                idx_j = char_index[char_j]
                if idx_i == idx_j:
                    # Count combinations within the same character group
                    #matrix[idx_i][idx_j] = freq[char_i] * (freq[char_i] - 1) // 2
                    matrix[idx_i][idx_j] = freq[char_i] * (freq[char_i] - 1) 
                else:
                    # Count combinations between different character groups
                    matrix[idx_i][idx_j] = freq[char_i] * freq[char_j]
        return matrix
    #################################################################################################################

    #################################################################################################################
    def get_contact_lib(self, gene_mask, tag):
        contact_outfile = f'{self.FrequencyGeneratorContactsOutpath}contacts_{tag}.csv'
        # Load only those contact files in the gene_mask
        contact_files = [f for f in glob(self.contact_files) if f.split('/')[-1].split('_')[0] in gene_mask]
        #print(f'# contact_files: {len(contact_files)}')

        contact_df = {'gene':[], 'contact_resids':[], 'contact_resnames':[], 'contact_mapping2matrix':[]}

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
                        contact_df['contact_mapping2matrix'] += [(x, y)]
                        #print(fi, f, i, iAA, j, jAA, x, y)
            

            contact_df = pd.DataFrame(contact_df)
            contact_df.to_csv(contact_outfile, sep='|', index=False)
            print(f'SAVED: {contact_outfile}')

        print(f'LOADING: {contact_outfile}')
        logger.info(f'LOADING: {contact_outfile}')
        contact_df = pd.read_csv(contact_outfile, sep='|')

        return contact_df
    #################################################################################################################

    #################################################################################################################
    def get_GE_contact_lib(self, gene_mask, uent_files, contact_df, rep_genes, tag):
        GE_contact_outfile = f'{self.FrequencyGeneratorContactsOutpath}GE_contacts_{tag}.csv'

        if not os.path.exists(GE_contact_outfile):
            # Load unique ent files
            uent_files = [f for f in glob(uent_files) if f.split('/')[-1].split('_')[0] in gene_mask]
            print(f'Number of uent files: {len(self.uent_files)}')

            GE_contact_df = []
            # Check if the file already exists. If so tell user and skip.
            for fi, f in enumerate(uent_files):
                gene = f.split('/')[-1].split('_')[0]
                gene_rep_df = rep_genes[rep_genes['gene'] == gene]
                pdb = gene_rep_df['pdb'].values[0]
                chain = gene_rep_df['chain'].values[0]
                print(fi, f, gene, pdb, chain)
                gene_contacts = contact_df[contact_df['gene'] == gene]
                print(gene_contacts)

                #uent_df = pd.read_csv(f, sep='|')['contacts'].values
                uent_df = pd.read_csv(f, sep='|')
                uent_df = uent_df[uent_df['CCBond'] == False]
                uent_df = uent_df['contacts'].values
                uent_df = ';'.join(uent_df)
                #print(uent_df)
                
                loop_contacts = [[x.split('-')[0], x.split('-')[1]] for x in uent_df.split(';') if len(x.split('-')) == 2]
                #print(loop_contacts)
                
                loop_contact_strs = []
                for loop_contact in loop_contacts:
                    i = int(loop_contact[0])
                    j = int(loop_contact[1])
                    #print(loop_contact, i, j)

                    loop_contact_str1 = f'({loop_contact[0]}, {loop_contact[1]})'
                    loop_contact_str2 = f'({loop_contact[1]}, {loop_contact[0]})'
                    loop_contact_strs += [loop_contact_str1, loop_contact_str2]
                    #print(loop_contact, loop_contact_str1, loop_contact_str2)
                loop_contact_info = gene_contacts[gene_contacts['contact_resids'].isin(loop_contact_strs)]
                #print(loop_contact_info)
                GE_contact_df += [loop_contact_info]
            
            GE_contact_df = pd.concat(GE_contact_df)
            GE_contact_df.to_csv(GE_contact_outfile, sep='|', index=False)
            print(f'SAVED: {GE_contact_outfile}')

        print(f'LOADING: {GE_contact_outfile}')
        logger.info(f'LOADING: {GE_contact_outfile}')
        GE_contact_df = pd.read_csv(GE_contact_outfile, sep='|')

        return GE_contact_df
    #################################################################################################################

    #################################################################################################################
    def Fc_ab(self, contact_df, tag):
        """
        Fc_ab = Nc(a,b)/Nc 
        this is the frequence of a given pair of amino acids (a) and (b) being in contact across the proteome
        where Nc(a,b) is the number of time (a) and (b) where in contact across the proteome, and Nc is the total number of contacts across the proteome
        """

        # Initialize the Fc_ab outpath
        Fc_ab_outfile = f'{self.FrequencyGeneratorOutpath}Fc_ab_{tag}.csv'
        #Nc_ab = np.zeros((20,20))
        Nc_ab = np.full((20,20), self.floor)

        # Check if the file already exists. If so tell user and skip.
        # Else make all pairs
        #Ess_contacts:        contact_resids contact_resnames contact_mapping2matrix
        #print(f'Generating Fc_ab... This may take some time')
        for m in contact_df['contact_mapping2matrix'].values:
            x, y = ast.literal_eval(m)

            if x != y:
                Nc_ab[x,y] += 1 
                Nc_ab[y,x] += 1 
            elif x == y:
                Nc_ab[x,y] += 1 

        # Create a mask for the lower triangle including the diagonal
        mask = np.tril(np.ones(Nc_ab.shape, dtype=bool))
        lower_triangle = Nc_ab[mask]
        zero_indices = np.argwhere(lower_triangle == 0)
        row_col_indices = np.argwhere(mask)
        zeros = [row_col_indices[i] for i in zero_indices.flatten()]

        # Calculate Fc_ab = N_ab/N and make it into a pandas df
        Nc = len(contact_df['contact_mapping2matrix'].values) 

        Fc_ab = Nc_ab/Nc
        AA_categories = self.AAs2idx.keys()
        Fc_ab_df = pd.DataFrame(Fc_ab, columns = AA_categories)
        Fc_ab_df['AA'] = AA_categories

        if tag != 'None':
            print(f'Nc_ab:\n{Nc_ab}')
            logger.info(f'Nc_ab:\n{Nc_ab}')

            # get zeros locations
            print('zeros: {zeros}')

            print(f'{tag} Nc: {Nc}')
            logger.info(f'{tag} Nc: {Nc}')
            print(f'Fc_ab_df:\n{Fc_ab_df}')
            logger.info(f'Fc_ab_df:\n{Fc_ab_df}')

            Fc_ab_df.to_csv(Fc_ab_outfile, sep='|', index=False)
            #print(f'SAVED: {Fc_ab_outfile}')
            logger.info(f'SAVED: {Fc_ab_outfile}')

        return Fc_ab_df, zeros
    #################################################################################################################

    #################################################################################################################
    def Get_contact_energy_maps(self, F_ab, Fc_ab, FcG_ab, tag):
        #print(f'{"#"*100}\nGet_contact_energy_maps')

        self.Enorm_ab_outfile = f'{self.FrequencyGeneratorOutpath}Enorm_ab_{tag}.csv'
        self.Ege_ab_outfile = f'{self.FrequencyGeneratorOutpath}Ege_ab_{tag}.csv'
        self.deltaE_outfile = f'{self.FrequencyGeneratorOutpath}deltaE_{tag}.csv'

        F_ab_sub = F_ab.loc[:, F_ab.columns != 'AA']
        Fc_ab_sub = Fc_ab.loc[:, Fc_ab.columns != 'AA']
        FcG_ab_sub = FcG_ab.loc[:, FcG_ab.columns != 'AA']
        if tag != 'None':
            print(f'Fc_ab_sub:\n{Fc_ab_sub}')
            print(f'FcG_ab_sub:\n{FcG_ab_sub}')

        Enorm_ab = -100*np.log10(Fc_ab_sub/F_ab_sub)
        Ege_ab = -100*np.log10(FcG_ab_sub/F_ab_sub)
        #deltaE = Ege_ab - Enorm_ab
        deltaE = -100*np.log10(FcG_ab_sub/Fc_ab_sub)
        
        # add AA column to all df and print
        Enorm_ab['AA'] = F_ab['AA']
        Ege_ab['AA'] = F_ab['AA']
        deltaE['AA'] = F_ab['AA']

        if tag != 'None':
            pd.set_option('display.float_format', lambda x: '%.2f' % x)
            #print(f'Enorm_ab:\n{Enorm_ab}')
            Enorm_ab.to_csv(self.Enorm_ab_outfile, sep='|', index=False)
            print(f'SAVED: {self.Enorm_ab_outfile}')
            logger.info(f'SAVED: {self.Enorm_ab_outfile}')

            #print(f'Ege_ab:\n{Ege_ab}')
            Ege_ab.to_csv(self.Ege_ab_outfile, sep='|', index=False)
            print(f'SAVED: {self.Ege_ab_outfile}')
            logger.info(f'SAVED: {self.Ege_ab_outfile}')

            #print(f'deltaE:\n{deltaE}')
            deltaE.to_csv(self.deltaE_outfile, sep='|', index=False)
            print(f'SAVED: {self.deltaE_outfile}')
            logger.info(f'SAVED: {self.deltaE_outfile}')

        return deltaE
    #################################################################################################################

    #################################################################################################################
    def deltaDeltaE(self, Ess_GT_FcG_ab, NonEss_GT_FcG_ab, tag):
        
        deltaDeltaE_outfile = f'{self.FrequencyGeneratorOutpath}deltaDeltaE_{tag}.csv'
        #print(f'Ess_deltaE:\n{Ess_deltaE}')
        #print(f'NonEss_deltaE:\n{NonEss_deltaE}')

        Ess_sub = Ess_GT_FcG_ab.loc[:, Ess_GT_FcG_ab.columns != 'AA']
        NonEss_sub = NonEss_GT_FcG_ab.loc[:, NonEss_GT_FcG_ab.columns != 'AA']
        
        deltaDeltaE = -100*np.log10(Ess_sub/NonEss_sub)
        #print(f'deltaDeltaE:\n{deltaDeltaE}')

        deltaDeltaE['AA'] = Ess_GT_FcG_ab['AA'] 
        deltaDeltaE_full = deltaDeltaE.copy()
        #print(f'deltaDeltaE_full:\n{deltaDeltaE_full}')

        if tag != 'None':
            deltaDeltaE.to_csv(deltaDeltaE_outfile, sep='|', index=False)
            print(f'SAVED: {deltaDeltaE_outfile}')
            logger.info(f'SAVED: {deltaDeltaE_outfile}')

        return deltaDeltaE
    #################################################################################################################

    #################################################################################################################
    def permute(self, p, Ess_GE_contacts, NonEss_GE_contacts, All_Fc_ab, Ess_GT_deltaE, NonEss_GT_deltaE, deltaDeltaE_GT):
        #print(f'\nPERMUTATION {p}')
        #logger.info(f'PERMUTATION {p}')

        # Combine the DataFrames
        combined_df = pd.concat([Ess_GE_contacts, NonEss_GE_contacts], ignore_index=True)

        # Shuffle the combined DataFrame
        shuffled_df = shuffle(combined_df).reset_index(drop=True)

        # Split the shuffled DataFrame back into two DataFrames with the original sizes
        Ess_GE_contacts_permuted = shuffled_df.iloc[:len(Ess_GE_contacts)]
        NonEss_GE_contacts_permuted = shuffled_df.iloc[len(Ess_GE_contacts):]
        #print(f'Ess_GE_contacts_permuted:\n{Ess_GE_contacts_permuted}')
        #print(f'NonEss_GE_contacts_permuted:\n{NonEss_GE_contacts_permuted}')

        # save the permutated contacts
        #Ess_GE_contacts_permuted.to_csv(Ess_GE_contact_outfile, sep='|', index=False)
        #print(f'SAVED: {Ess_GE_contact_outfile}')
        #NonEss_GE_contacts_permuted.to_csv(NonEss_GE_contact_outfile, sep='|', index=False)
        #print(f'SAVED: {NonEss_GE_contact_outfile}')

        # calculate the permuted Essential gene deltaE
        Ess_FcG_ab_permute, _ = self.Fc_ab(Ess_GE_contacts_permuted, 'None')
        #Ess_FcG_ab_permute = self.Fc_ab(Ess_GE_contacts_permuted, f'Ess_GE_permute{p}')
        #print(f'Ess_FcG_ab_permute:\n{Ess_FcG_ab_permute}')
        
        Ess_deltaE_permute = self.Get_contact_energy_maps(self.F_ab_df, All_Fc_ab, Ess_FcG_ab_permute, 'None')
        #Ess_deltaE_permute = self.Get_contact_energy_maps(self.F_ab_df, All_Fc_ab, Ess_FcG_ab_permute, f'Ess_permute{p}')
        #print(f'Ess_deltaE_permute:\n{Ess_deltaE_permute}')

        # calculate the permuted NonEssential gene deltaE
        NonEss_FcG_ab_permute, _ = self.Fc_ab(NonEss_GE_contacts_permuted, 'None')
        #NonEss_FcG_ab_permute = self.Fc_ab(NonEss_GE_contacts_permuted, f'NonEss_GE_permute{p}')
        #print(f'NonEss_FcG_ab_permute:\n{NonEss_FcG_ab_permute}')
        
        NonEss_deltaE_permute = self.Get_contact_energy_maps(self.F_ab_df, All_Fc_ab, NonEss_FcG_ab_permute, 'None')
        #NonEss_deltaE_permute = self.Get_contact_energy_maps(self.F_ab_df, All_Fc_ab, NonEss_FcG_ab_permute, f'NonEss_permute{p}')
        #print(f'NonEss_deltaE_permute:\n{NonEss_deltaE_permute}')

        # Get permuted deltaDeltaE
        deltaDeltaE_permute = self.deltaDeltaE(Ess_FcG_ab_permute, NonEss_FcG_ab_permute, 'None')
        #deltaDeltaE_permute = self.deltaDeltaE(Ess_deltaE_permute, NonEss_deltaE_permute, f'permute{p}')
        #print(f'deltaDeltaE_permute {p}:\n{deltaDeltaE_permute}')

        # caclualte deltaDeltaE pvalues
        premuted_dfs = [Ess_deltaE_permute.loc[:, Ess_deltaE_permute.columns != 'AA']]
        Ess_deltaE_pvalue = self.count_conditionally(Ess_GT_deltaE.loc[:, Ess_GT_deltaE.columns != 'AA'], premuted_dfs)

        # caclualte deltaDeltaE pvalues
        premuted_dfs = [NonEss_deltaE_permute.loc[:, NonEss_deltaE_permute.columns != 'AA']]
        NonEss_deltaE_pvalue = self.count_conditionally(NonEss_GT_deltaE.loc[:, NonEss_GT_deltaE.columns != 'AA'], premuted_dfs)

        # caclualte deltaDeltaE pvalues
        premuted_dfs = [deltaDeltaE_permute.loc[:, deltaDeltaE_permute.columns != 'AA']]
        deltaDeltaE_pvalue = self.count_conditionally(deltaDeltaE_GT.loc[:, deltaDeltaE_GT.columns != 'AA'], premuted_dfs)

        return Ess_deltaE_pvalue, NonEss_deltaE_pvalue, deltaDeltaE_pvalue
    #################################################################################################################

    #################################################################################################################
    def pvalues(self, result_df, df_type, AA):
        print(f'Calculate Pvalues {df_type}')
        #print(f'GT_df:\n{GT_df}')
        pvalue_df_outfile = f'{self.FrequencyGeneratorOutpath}{df_type}_pvalues.csv'
        FDR_df_outfile = f'{self.FrequencyGeneratorOutpath}{df_type}_FDR_pvalues.csv'

        # Zero the upper triangle
        result_df = self.zero_upper_triangle(result_df)
        result_df['AA'] = AA
        result_df.to_csv(pvalue_df_outfile, sep='|', index=False)
        print(f'SAVED: {pvalue_df_outfile}')
        logger.info(f'SAVED: {pvalue_df_outfile}')

        # Apply FDR correction
        corrected_df = self.apply_fdr_correction(result_df.loc[:, result_df.columns != 'AA'], alpha=0.05)

        corrected_df['AA'] = AA
        print(f'corrected_df:\n{corrected_df.to_string()}')
        corrected_df.to_csv(FDR_df_outfile, sep='|', index=False)
        print(f'SAVED: {FDR_df_outfile}')
        logger.info(f'SAVED: {FDR_df_outfile}')
    #################################################################################################################

    #################################################################################################################
    def count_conditionally(self, ground_truth_df, resampled_dfs):
        # Initialize a DataFrame to hold the counts, same shape as ground_truth_df
        result_df = pd.DataFrame(0, index=ground_truth_df.index, columns=ground_truth_df.columns)

        for resampled_df in resampled_dfs:
            # Increment the result_df based on the conditions
            #print(ground_truth_df)
            #print(resampled_df)

            # Iterate through each element of the DataFrame
            for row in ground_truth_df.index:
                for col in ground_truth_df.columns:
                    ground_value = ground_truth_df.at[row, col]
                    resampled_value = resampled_df.at[row, col]

                    if pd.notna(ground_value) and pd.notna(resampled_value):  # Check for non-NaN values
                        if ground_value > 0 and resampled_value > ground_value:
                            result_df.at[row, col] += 1
                        elif ground_value < 0 and resampled_value < ground_value:
                            result_df.at[row, col] += 1
                        elif ground_value == 0 and resampled_value != ground_value:
                            result_df.at[row, col] += 1

        return result_df
    #################################################################################################################

    #################################################################################################################
    # Function to zero the upper triangle of the matrix
    def zero_upper_triangle(self, matrix):
        matrix = matrix.copy()
        for i in range(len(matrix)):
            for j in range(i+1, len(matrix)):
                matrix.iat[i, j] = 0
        return matrix
    #################################################################################################################

    #################################################################################################################
    # Function to apply FDR correction
    def apply_fdr_correction(self, matrix, alpha=0.05):
        # Flatten the lower triangle of the matrix and keep non-zero values

        lower_triangle = matrix.where(np.tril(np.ones(matrix.shape), 0).astype(bool))
        p_values = lower_triangle.stack().values
        print(f'p_values: {p_values}')
        p_values = np.where(p_values == 0, 0.000000000000001, p_values)
        print(f'p_values: {p_values}')

        # Apply FDR correction
        corrected_pvals = false_discovery_control(p_values, method='bh')
        #_, corrected_pvals = fdrcorrection(p_values, alpha=alpha)
        print(f'corrected_pvals: {corrected_pvals}')

        # Reshape corrected p-values back to the lower triangle format
        corrected_matrix = pd.DataFrame(np.zeros(matrix.shape), index=matrix.index, columns=matrix.columns)
        corrected_matrix = corrected_matrix.where(np.triu(np.ones(corrected_matrix.shape), 0).astype(bool))  # Upper triangle should remain zero
        idx = 0
        for i in range(0, len(matrix)):
            for j in range(i + 1):
                corrected_matrix.iat[i, j] = corrected_pvals[idx]
                idx += 1

        return corrected_matrix
    #################################################################################################################


#################################################################################################################

def main():
    """
    Main function to control workflow. 
    (1) parse user arguments 
    (2) making logger file
    (3) attempt to make F_ab, Fc_ab, FcG_ab(Ec1,...Ec4)
    """

    # Parse the user supplied arguments
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-f", "--fasta_files", type=str, required=True, help=f"path to fasta files for all genes in dataset")
    parser.add_argument("-g", "--rep_gene_list", type=str, required=True, help=f"path to representative gene list used in the calculation of F_ab")
    parser.add_argument("-Eg", "--Ess_gene_list", type=str, required=True, help=f"path to Essential gene list used for mask in Fc_ab and FcG_ab calcs")
    parser.add_argument("-NEg", "--NonEss_gene_list", type=str, required=True, help=f"path to Non-Essential gene list used for mask in Fc_ab and FcG_ab calcs")
    parser.add_argument("-l", "--log_file", type=str, required=True, help="Path to logger file")
    parser.add_argument("-c", "--contact_files", type=str, required=True, help="path to native contact files")
    parser.add_argument("-r", "--resFeat_files", type=str, required=True, help="path to residue Feature files")
    parser.add_argument("-e", "--uent_files", type=str, required=True, help="path to unique entanglement files")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="path to output directory. will be made if doesnt exist")
    parser.add_argument("-p", "--num_permute", type=int, required=True, help="Number of permutations")
    parser.add_argument("--floor", type=int, required=True, help="Floor of Fc_ab contact matrix")
    args = parser.parse_args()


    fasta_files = args.fasta_files
    contact_files = args.contact_files
    resFeat_files = args.resFeat_files
    uent_files = args.uent_files
    log_file = args.log_file
    rep_gene_list = args.rep_gene_list
    outpath = args.outpath
    Ess_gene_list = args.Ess_gene_list
    NonEss_gene_list = args.NonEss_gene_list
    num_permute = args.num_permute
    floor = args.floor

    # Setup logger configuration
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s %(message)s') 
    global logger
    logger = logging.getLogger(__name__)
    logger.info(f'{"#"*100}\nNEW RUN')

    # Make outpath if it doesnt exists
    if not os.path.exists(outpath):
        os.makedirs(outpath)
        print(f'Made directory: {outpath}')

    # Initalize the FrequencyGenerator class object
    Fgen = FrequencyGenerator(outpath, fasta_files, contact_files, uent_files, floor)
    print(f'Fgen: {Fgen}')

    # Load the representative genes across the protome
    rep_genes = pd.read_csv(rep_gene_list, sep=' ')
    print(f'rep_genes: {rep_genes} {len(rep_genes)}')

    # Get F_ab for all 1294 representative genes
    Fgen.F_ab(fasta_files, rep_genes['gene'].values)

    # Get Fc_ab for all genes
    All_contacts = Fgen.get_contact_lib(rep_genes['gene'].values, 'All')
    print(f'All_contacts:\n{All_contacts}')

    All_Fc_ab, All_Fc_ab_zeros = Fgen.Fc_ab(All_contacts, 'All')
    print(f'All_Fc_ab:\n{All_Fc_ab.to_string()}')

    # Load Essential gene info and get ground truth (GT) 
    # (2) get Ess contacts
    # (3) get Ess contact frequency matrix Fc_ab
    # (4) get Ess GE contact frequency matrix FcG_ab
    # (5) get the final Ess Ege matrix
    Ess_gene_mask = np.loadtxt(Ess_gene_list, dtype=str)
    print(f'Ess_gene_mask: {Ess_gene_mask} {len(Ess_gene_mask)}')

    Ess_contacts = Fgen.get_contact_lib(Ess_gene_mask, 'Ess')
    print(f'Ess_contacts:\n{Ess_contacts}')

    #Ess_GT_Fc_ab = Fgen.Fc_ab(Ess_contacts, 'Ess_GT')
    #print(f'Ess_GT_Fc_ab:\n{Ess_GT_Fc_ab}')

    Ess_GE_contacts = Fgen.get_GE_contact_lib(Ess_gene_mask, uent_files, Ess_contacts, rep_genes, 'Ess_GE_GT')
    print(f'Ess_GE_contacts:\n{Ess_GE_contacts}')

    Ess_GT_FcG_ab, Ess_GT_FcG_ab_zeros = Fgen.Fc_ab(Ess_GE_contacts, 'Ess_GE_GT')
    print(f'Ess_GT_FcG_ab:\n{Ess_GT_FcG_ab.to_string()}')
    print(f'Ess_GT_FcG_ab_zeros: {Ess_GT_FcG_ab_zeros}')

    # Calulate Enorm_ab, Ege_ab, and deltaE 
    #Ess_GT_deltaE = Fgen.Get_contact_energy_maps(Fgen.F_ab_df, Ess_GT_Fc_ab, Ess_GT_FcG_ab, 'Ess_GT')
    Ess_GT_deltaE = Fgen.Get_contact_energy_maps(Fgen.F_ab_df, All_Fc_ab, Ess_GT_FcG_ab, 'Ess_GT')
    print(f'Ess_GT_deltaE:\n{Ess_GT_deltaE.to_string()}')

    # Load Non-Essential gene info and get ground truth (GT) 
    # (2) get Ess contacts
    # (3) get Ess contact frequency matrix Fc_ab
    # (4) get Ess GE contact frequency matrix FcG_ab
    # (5) get the final Ess Ege matrix
    NonEss_gene_mask = np.loadtxt(NonEss_gene_list, dtype=str)
    print(f'NonEss_gene_mask: {NonEss_gene_mask} {len(NonEss_gene_mask)}')

    NonEss_contacts = Fgen.get_contact_lib(NonEss_gene_mask, 'NonEss')
    print(f'NonEss_contacts:\n{NonEss_contacts}')

    #NonEss_GT_Fc_ab = Fgen.Fc_ab(NonEss_contacts, 'NonEss_GT')
    #print(f'NonEss_GT_Fc_ab:\n{NonEss_GT_Fc_ab}')

    NonEss_GE_contacts = Fgen.get_GE_contact_lib(NonEss_gene_mask, uent_files, NonEss_contacts, rep_genes, 'NonEss_GE_GT')
    print(f'NonEss_GE_contacts:\n{NonEss_GE_contacts}')

    NonEss_GT_FcG_ab, NonEss_GT_FcG_ab_zeros = Fgen.Fc_ab(NonEss_GE_contacts, 'NonEss_GE_GT')
    print(f'NonEss_GT_FcG_ab:\n{NonEss_GT_FcG_ab.to_string()}')
    print(f'NonEss_GT_FcG_ab_zeros: {NonEss_GT_FcG_ab_zeros}')

    # Calulate Enorm_ab, Ege_ab, and deltaE 
    #NonEss_GT_deltaE = Fgen.Get_contact_energy_maps(Fgen.F_ab_df, NonEss_GT_Fc_ab, NonEss_GT_FcG_ab, 'NonEss_GT')
    NonEss_GT_deltaE = Fgen.Get_contact_energy_maps(Fgen.F_ab_df, All_Fc_ab, NonEss_GT_FcG_ab, 'NonEss_GT')
    print(f'NonEss_GT_deltaE:\n{NonEss_GT_deltaE.to_string()}')

    # Get GT deltaDeltaE
    #deltaDeltaE = Fgen.deltaDeltaE(Ess_GT_deltaE, NonEss_GT_deltaE, 'GT')
    deltaDeltaE = Fgen.deltaDeltaE(Ess_GT_FcG_ab, NonEss_GT_FcG_ab, 'GT')
    print(f'deltaDeltaE:\n{deltaDeltaE.to_string()}')

    # add all the permuted contact matrixes together for each system
    pvalues = parallel_permute(Fgen, num_permute, Ess_GE_contacts, NonEss_GE_contacts, All_Fc_ab, Ess_GT_deltaE, NonEss_GT_deltaE, deltaDeltaE)
    Ess_deltaE_pvalues = reduce(lambda df1, df2: df1.add(df2, fill_value=0), [d[0] for d in pvalues])
    #Ess_deltaE_pvalues[Ess_GT_FcG_ab_zeros] = 0

    NonEss_deltaE_pvalues = reduce(lambda df1, df2: df1.add(df2, fill_value=0), [d[1] for d in pvalues])
    #NonEss_deltaE_pvalues[NonEss_GT_FcG_ab_zeros] = 0

    deltaDeltaE_pvalues = reduce(lambda df1, df2: df1.add(df2, fill_value=0), [d[2] for d in pvalues])
    #deltaDeltaE_pvalues[Ess_GT_FcG_ab_zeros] = 0
    #deltaDeltaE_pvalues[NonEss_GT_FcG_ab_zeros] = 0

    print(f'\n{"#"*50}\nEss_deltaE_pvalues:\n{Ess_deltaE_pvalues.to_string()}')
    Ess_deltaE_pvalues /= num_permute
    print(f'Ess_deltaE_pvalues:\n{Ess_deltaE_pvalues.to_string()}')
    AA = Ess_GT_deltaE['AA'].values
    Fgen.pvalues(Ess_deltaE_pvalues, 'deltaE_Ess', AA)

    print(f'\n{"#"*50}\nNonEss_deltaE_pvalues:\n{NonEss_deltaE_pvalues.to_string()}')
    NonEss_deltaE_pvalues /= num_permute
    print(f'NonEss_deltaE_pvalues:\n{NonEss_deltaE_pvalues.to_string()}')
    Fgen.pvalues(NonEss_deltaE_pvalues, 'deltaE_NonEss', AA)

    print(f'\n{"#"*50}\ndeltaDeltaE_pvalues:\n{deltaDeltaE_pvalues.to_string()}')
    deltaDeltaE_pvalues /= num_permute
    print(f'deltaDeltaE_pvalues:\n{deltaDeltaE_pvalues.to_string()}')
    Fgen.pvalues(deltaDeltaE_pvalues, 'deltaDeltaE', AA)



def parallel_permute(Fgen, num_permute, Ess_GE_contacts, NonEss_GE_contacts, All_Fc_ab, Ess_GT_deltaE, NonEss_GT_deltaE, deltaDeltaE):
    pool = mp.Pool(processes=10)  # Create a pool with 4 processes

    # Create the arguments to pass to the function in parallel
    args = [(p, Ess_GE_contacts, NonEss_GE_contacts, All_Fc_ab, Ess_GT_deltaE, NonEss_GT_deltaE, deltaDeltaE) for p in range(num_permute)]

    # Use pool.starmap to execute the function in parallel
    pvalues = pool.starmap(Fgen.permute, args)

    pool.close()  # Close the pool to free resources
    pool.join()   # Ensure all processes have completed


    return pvalues


start_time = time.time()
if __name__ == "__main__":
    main()
print(f'NORMAL TERMINATION: {time.time() - start_time}')
