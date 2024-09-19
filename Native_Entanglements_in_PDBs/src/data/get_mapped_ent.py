import sys, os, re, time, logging
import glob
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

class Analyzer:
    """
    """

    def __init__(self, outpath, uent_files, mapping):
        """
        """

        # Make outpath if it doesnt exists
        self.outpath = outpath
        self.uent_files = uent_files
        self.mapping = mapping
        if self.mapping != 'None':
            self.mapping_files = glob.glob(os.path.join(self.mapping, '*_resid_mapping.txt'))
        else:
            self.mapping_files = 'None'

    #################################################################################################################
    #################################################################################################################
    def get_GE_contact_lib(self,):

        # Load unique ent files
        uent_files = glob.glob(os.path.join(self.uent_files, '*'))
        print(f'Number of uent files: {len(self.uent_files)}')


        # Check if the file already exists. If so tell user and skip.
        for f_i, f in enumerate(uent_files):
            #Ecoli/EXP/unmapped_GE/P00350-3FWN_A_GE.txt
            gene = f.split('/')[-1].split('-')[0]
            pdb = f.split('/')[-1].split('-')[1].split('_')[0]
            chain = f.split('/')[-1].split('-')[1].split('_')[1]
            print(f'\n{"#"*50}\n', f_i, gene, pdb, chain, f)

            # Get mapping file
            if self.mapping != 'None':
                map_file = [f for f in self.mapping_files if f'{gene}-{pdb}_{chain}' in f]
                print(f'map_file: {map_file}')
                if len(map_file) != 0:
                    mapping = np.loadtxt(map_file[0], dtype='O')
                    mapping = np.vstack([x[1:] for x in mapping if ('Mapped' in x[0] or 'Modifed_Residue' in x[0] or 'Missense' in x[0])]).astype(int)
                    mapping_pdb2uniprot = {pdb:uni for pdb, uni in mapping}

            #ent_df = pd.read_csv(f, sep='|')
            ent_df = [x.strip('\n') for x in open(f, 'r').readlines()]

            outlines = []
            for line_i, line in enumerate(ent_df):
                #print(line)
                
                _, ijr, _, _, _ = line.split(' | ')
                ijr = ijr.strip('()')
                parts = ijr.split(', ', 2)
                i = int(parts[0])
                j = int(parts[1])
                r = ast.literal_eval(parts[2].replace(" ", ", "))
                crossings = []
                for cross in r:
                    for sub in ['--', '+-', '-+']:
                        if sub in cross:
                            print(f'FOUND unique crossing: {cross}')
                    cross = cross.lstrip('-')
                    cross = cross.lstrip('+')
                    crossings += [int(cross)]
                key_res = [i, j] + crossings
                #print(gene, ijr, i, j, r, crossings, key_res)
                
                mapped = True
                for res in key_res:
                    if res not in mapping_pdb2uniprot:
                        print(f'Res: {res} not mapped! this entanglement will be discarded from {gene} {pdb} {chain}')
                        mapped = False
                                       
                if mapped:
                    outlines += [line]
                else:
                    print(f'Entanglement {line_i} {line} was not mapped for {gene} {pdb} {chain}')
           
            if len(outlines) != 0:
                outfile = os.path.join(self.outpath, f'{gene}_{pdb}_{chain}_GE.txt')
                np.savetxt(outfile, outlines, fmt='%s')
                print(f'SAVED: {outfile}')
            else:
                print(f'No entanglements left for {gene} {pdb} {chain}')
                
#################################################################################################################

def main():
    """
    Main function to control workflow. 
    check unmapped raw entanglements for those that do not map
    """

    # Parse the user supplied arguments
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-e", "--uent_files", type=str, required=True, help="path to unique entanglement files")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="path to output directory. will be made if doesnt exist")
    parser.add_argument("-m", "--mapping", type=str, required=True, help="mapping files or None")
    args = parser.parse_args()

    uent_files = args.uent_files
    outpath = args.outpath
    mapping = args.mapping

    # Make outpath if it doesnt exists
    if not os.path.exists(outpath):
        os.makedirs(outpath)
        print(f'Made directory: {outpath}')

    # Initalize the FrequencyGenerator class object
    analysis = Analyzer(outpath, uent_files, mapping)
    print(f'analysis: {analysis}')

    # check entanglement files for nonmapped entanglements
    analysis.get_GE_contact_lib()






start_time = time.time()
if __name__ == "__main__":
    main()
print(f'NORMAL TERMINATION: {time.time() - start_time}')
