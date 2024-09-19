import sys, os, re, time, logging
import glob
import shutil
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
from Bio.PDB import PDBParser
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

    def __init__(self, outpath, uent_files, AFpdbs):
        """
        """

        # Make outpath if it doesnt exists
        self.outpath = outpath
        self.uent_files = uent_files
        self.AFpdbs = AFpdbs


    #################################################################################################################
    #################################################################################################################
    def get_HQ_AF_ents(self,):

        # Load unique ent files
        uent_files = glob.glob(os.path.join(self.uent_files, '*'))
        print(f'Number of uent files: {len(self.uent_files)}')

        self.AFpdbs = glob.glob(os.path.join(self.AFpdbs, '*'))
        for pdb in self.AFpdbs:
            gene = pdb.split('/')[-1].replace('.pdb', '')
            #if gene != 'P22523':
            #    continue
            avg_pLDDT, pLDDT_df = self.average_pLDDT(pdb)
            print(f'\n{"#"*50}\n', pdb, gene, avg_pLDDT)
            #print(pLDDT_df)

            # check if the AF structure has a <pLDDT> >= 70
            if avg_pLDDT >= 70:
                print(f'HQ AF structure found {gene} with {avg_pLDDT}')
                #check for uent file
                ent_file = [f for f in uent_files if f'{gene}_GE.txt' in f]

                if len(ent_file) != 0:
                    ent_file = ent_file[0]
                    print(f'ENT file found: {ent_file}')

                    ent_data = self.load_ent_file(ent_file)
                    outlines = []
                    for ent_i, ent in ent_data.items():
                        print(ent_i, ent)
                        line = ent['line']
                        i = ent['NC'][0]
                        j = ent['NC'][1]
                        # (1) check if both i and j have pLDDt >= 70. if so continue else completely ignore the ent
                        NC_pLDDT = pLDDT_df[pLDDT_df['resid'].isin(ent['NC'])]['pLDDT'].values
                        if all(NC_pLDDT >= 70):
                            print(f'Native contact pLDDT are >= 70 {NC_pLDDT}')
                        else:
                            print(f'Native contact pLDDT are < 70 {NC_pLDDT}')
                            continue
                         

                        # (2) starting from the loop base get the set of ordered crossings that have pLDDT > 70 and discard any after the first crossings that fails this. 
                        print(f'Getting HQ N-terminal entanglements')
                        Ncrossings_resids = ent['Ncrossings_resids']
                        Ncrossings = ent['Ncrossings']
                        Ncrossings_signs = ent['Ncrossings_signs']
                        Ncrossings_resids, Ncrossings, Ncrossings_signs, Ndup_flag = self.remove_duplicates(Ncrossings_resids, Ncrossings, Ncrossings_signs)
                        HQ_Ncrossings = []
                        HQ_Ncrossings_resids = []
                        HQ_Ncrossings_signs = []
                        if len(Ncrossings_resids) != 0:
                            sorted_indices = np.argsort(Ncrossings_resids)[::-1]  # [::-1] reverses the order
                            Ncrossings_resids = [Ncrossings_resids[i] for i in sorted_indices]
                            Ncrossings = [Ncrossings[i] for i in sorted_indices]
                            Ncrossings_signs = [Ncrossings_signs[i] for i in sorted_indices]
                            Ncrossings_pLDDTs = pLDDT_df[pLDDT_df['resid'].isin(Ncrossings_resids)]['pLDDT'].values
                            print(Ncrossings_resids, Ncrossings, Ncrossings_signs, Ncrossings_pLDDTs)
                            for cross_i, cross in enumerate(Ncrossings_resids):
                                if Ncrossings_pLDDTs[cross_i] >=70:
                                    HQ_Ncrossings += [Ncrossings[cross_i]]
                                    HQ_Ncrossings_resids += [cross]
                                    HQ_Ncrossings_signs += [Ncrossings_signs[cross_i]]
                                else:
                                    break
                            
                            # check that the remaining HQ Nterminal crossigns are not a slipknot
                            if sum(HQ_Ncrossings_signs) == 0:
                                print(f'SlipKnot foundin N terminus after HQ search')
                                HQ_Ncrossings = []
                                HQ_Ncrossings_resids = []
                                HQ_Ncrossings_signs = []
                        
                        print(f'Getting HQ C-terminal entanglements')
                        Ccrossings_resids = ent['Ccrossings_resids']
                        Ccrossings = ent['Ccrossings']
                        Ccrossings_signs = ent['Ccrossings_signs']
                        Ccrossings_resids, Ccrossings, Ccrossings_signs, Cdup_flag = self.remove_duplicates(Ccrossings_resids, Ccrossings, Ccrossings_signs)
                        HQ_Ccrossings = []
                        HQ_Ccrossings_resids = []
                        HQ_Ccrossings_signs = []
                        if len(Ccrossings_resids) != 0:
                            sorted_indices = np.argsort(Ccrossings_resids) # [::-1] reverses the order
                            Ccrossings_resids = [Ccrossings_resids[i] for i in sorted_indices]
                            Ccrossings = [Ccrossings[i] for i in sorted_indices]
                            Ccrossings_signs = [Ccrossings_signs[i] for i in sorted_indices]
                            Ccrossings_pLDDTs = pLDDT_df[pLDDT_df['resid'].isin(Ccrossings_resids)]['pLDDT'].values
                            print(Ccrossings_resids, Ccrossings, Ccrossings_signs, Ccrossings_pLDDTs)
                            for cross_i, cross in enumerate(Ccrossings_resids):
                                if Ccrossings_pLDDTs[cross_i] >=70:
                                    HQ_Ccrossings += [Ccrossings[cross_i]]
                                    HQ_Ccrossings_resids += [cross]
                                    HQ_Ccrossings_signs += [Ccrossings_signs[cross_i]]
                                else:
                                    break

                            # check that the remaining HQ Cterminal crossigns are not a slipknot
                            if sum(HQ_Ccrossings_signs) == 0:
                                print(f'SlipKnot foundin N terminus after HQ search')
                                HQ_Ccrossings = []
                                HQ_Ccrossings_resids = []
                                HQ_Ccrossings_signs = []

                        HQ_crossings = HQ_Ccrossings + HQ_Ncrossings
                        HQ_crossings = np.asarray(HQ_crossings, dtype=str)
                        print(f'HQ_crossings: {HQ_crossings}')
                        if len(HQ_crossings) != 0:
                            new_line = line.split('|')
                            new_ijr = f' ({i}, {j}, {HQ_crossings}) '
                            new_line[1] = new_ijr
                            new_line = '|'.join(new_line)
                            print(new_line)
                            outlines += [new_line]
                        else:
                            print(f'No crossing remaining in either the N or C terminus after removals. ENT will be ignored')

                    # Save outlines
                    if len(outlines) != 0:
                        outfile = os.path.join(self.outpath, f'{gene}_GE.txt')
                        np.savetxt(outfile, outlines, fmt='%s')
                        print(f'SAVED: {outfile}')
                    else:
                        print(f'No entanglements left for {gene}')
                    
    def remove_duplicates(self, crossings_resids, crossings, crossings_signs):
        # Create a list to keep track of unique elements in A and their corresponding elements in B
        unique_crossings_resids = []
        unique_crossings = []
        unique_crossings_signs = []

        # Dictionary to count occurrences of elements in A
        counts_A = {item: crossings_resids.count(item) for item in crossings_resids}

        # Remove elements from both lists where duplicates are found in A
        dup_flag = False
        for i in range(len(crossings_resids)):
            if counts_A[crossings_resids[i]] == 1:  # If the element in A is unique
                unique_crossings_resids.append(crossings_resids[i])
                unique_crossings.append(crossings[i])
                unique_crossings_signs.append(crossings_signs[i])
            else:
                unique_crossings_resids = []
                unique_crossings = []
                unique_crossings_signs = []
                dup_flag = True
                break

        return unique_crossings_resids, unique_crossings, unique_crossings_signs, dup_flag


    def load_ent_file(self, ent_file):
                  
        ent_df = [x.strip('\n') for x in open(ent_file, 'r').readlines()]

        #  preprocess to ensure there are no lines that got wrapped weird
        ent_lines = []
        i = 0
        while i < len(ent_df):
            #for line_i, line in enumerate(ent_df):
            #print(line, line[-4:])
            line = ent_df[i]

            if line.endswith('True') or line.endswith('False'):
                ent_lines += [line]
                i += 1
            else:
                #print(line)
                # find how many lines ahead dont start with chain
                offset = 0
                for future_line in ent_df[i+1:]:
                    if not future_line.startswith('Chain'):
                        offset += 1
                    else:
                        break
                #print(offset)
                new_line = line + ''.join(ent_df[i+1:i+offset+1])
                #print(new_line)
                ent_lines += [new_line]
                i += offset + 1


        # Loop through the lines of the entanglement file and get the key residue information
        out_ents = {}
        for line_i, line in enumerate(ent_lines):
            print(line, line[-4:])

            _, ijr, _, _, _ = line.split(' | ')
            ijr = ijr.strip('()')
            parts = ijr.split(', ', 2)
            i = int(parts[0])
            j = int(parts[1])

            r = ast.literal_eval(parts[2].replace(" ", ", "))
            Ncrossings = []
            Ncrossings_resids = []
            Ncrossings_signs = []
            Ccrossings = []
            Ccrossings_resids = []
            Ccrossings_signs = []

            for cross in r:
                for sub in ['--', '+-', '-+']:
                    if sub in cross:
                        print(f'FOUND unique crossing: {cross}')
                cross_sign = cross[0]
                cross_int = int(cross[1:])

                # check for N terminal crossings
                if cross_int < i:
                    Ncrossings += [cross]
                    Ncrossings_resids += [cross_int]
                    if cross_sign == '+':
                        Ncrossings_signs += [1]
                    elif cross_sign == '-':
                        Ncrossings_signs += [-1]

                # check for C terminal crossings
                elif cross_int > j:
                    Ccrossings += [cross]
                    Ccrossings_resids += [cross_int]
                    if cross_sign == '+':
                        Ccrossings_signs += [1]
                    elif cross_sign == '-':
                        Ccrossings_signs += [-1]

            out_ents[line_i] = {'line':line, 
                                'NC':[i,j], 
                                'Ncrossings':Ncrossings, 
                                'Ncrossings_resids':Ncrossings_resids,
                                'Ncrossings_signs':Ncrossings_signs,
                                'Ccrossings':Ccrossings,
                                'Ccrossings_resids':Ccrossings_resids,
                                'Ccrossings_signs':Ccrossings_signs}
        return out_ents   


    def average_pLDDT(self, pdb_filename):
        
        # Create a PDB parser
        parser = PDBParser(QUIET=True)
        
        # Parse the PDB structure
        structure = parser.get_structure('PDB_structure', pdb_filename)
        
        # List to hold pLDDT
        pLDDTs = []
        pLDDT_df = {'resid':[], 'pLDDT':[]}
        # Iterate over all atoms in the structure to extract pLDDT
        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        if atom.get_name() == 'CA':
                            pLDDTs.append(atom.get_bfactor())
                            pLDDT_df['resid'] += [residue.get_id()[1]]
                            pLDDT_df['pLDDT'] += [atom.get_bfactor()]
        
        # Calculate the average pLDDT
        if len(pLDDTs) > 0:
            avg_pLDDT = sum(pLDDTs) / len(pLDDTs)
            pLDDT_df = pd.DataFrame(pLDDT_df)
            return avg_pLDDT, pLDDT_df
        else:
            return None   
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
    parser.add_argument("-p", "--AFpdbs", type=str, required=True, help="path to AF pdbs files or None")
    args = parser.parse_args()

    uent_files = args.uent_files
    outpath = args.outpath
    AFpdbs = args.AFpdbs

    # Make outpath if it doesnt exists
    if not os.path.exists(outpath):
        os.makedirs(outpath)
        print(f'Made directory: {outpath}')

    # Initalize the FrequencyGenerator class object
    analysis = Analyzer(outpath, uent_files, AFpdbs)
    print(f'analysis: {analysis}')

    # check entanglement files for nonmapped entanglements
    analysis.get_HQ_AF_ents()






start_time = time.time()
if __name__ == "__main__":
    main()
print(f'NORMAL TERMINATION: {time.time() - start_time}')
