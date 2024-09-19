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

    def __init__(self, outpath, uent_files):
        """
        """

        # Make outpath if it doesnt exists
        self.outpath = outpath
        self.uent_files = uent_files

    #################################################################################################################
    #################################################################################################################
    def remove_slipknots(self,):

        # Load unique ent files
        uent_files = glob.glob(os.path.join(self.uent_files, '*'))
        print(f'Number of uent files: {len(self.uent_files)}')


        # Check if the file already exists. If so tell user and skip.
        for f_i, f in enumerate(uent_files):
            #Ecoli/EXP/mapped_GE/P00350_3FWN_A_GE.txt
            try:
                # assuming EXP PDB at first
                gene, pdb, chain, _ = f.split('/')[-1].split('_')
            except:
                # if failed assume AF structure
                gene = f.split('/')[-1].split('_')[0]
                pdb = 'AF'
                chain = 'A'

            #if gene != 'Q8NB78':
            #    continue
            print(f'\n{"#"*50}\n', f_i, gene, pdb, chain, f)

            #ent_df = pd.read_csv(f, sep='|')
            ent_df = [x.strip('\n') for x in open(f, 'r').readlines()]

            outlines = []
            for line_i, line in enumerate(ent_df):
                print(f'\n{line}')
                
                _, ijr, _, _, _ = line.split(' | ')
                ijr = ijr.strip('()')
                parts = ijr.split(', ', 2)
                i = int(parts[0])
                j = int(parts[1])
                r = ast.literal_eval(parts[2].replace(" ", ", "))

                # check if crossings signs add to 0 if there is more than 1 crossing
                if len(r) > 1:
                    N_crossings_signs = []
                    N_crossings = []
                    N_crossings_resids = []
                    C_crossings_signs = []
                    C_crossings = []
                    C_crossings_resids = []
                    for cross in r:
                        #print(cross)
                        cross_sign = cross[0]
                        cross_int = int(cross[1:])
                        #print(i, j, cross)

                        #check if the crossing in N terminal
                        if cross_int < i:
                            N_crossings += [cross]
                            N_crossings_resids += [cross_int]
                            if cross_sign == '+':
                                N_crossings_signs += [1]
                            elif cross_sign == '-':
                                N_crossings_signs += [-1]
                            else:
                                raise ValueError(f'The crossing sign was not + or - {cross}')

                        #check if the crossing is C terminal
                        if cross_int > j:
                            C_crossings += [cross]
                            C_crossings_resids += [cross_int]
                            if cross_sign == '+':
                                C_crossings_signs += [1]
                            elif cross_sign == '-':
                                C_crossings_signs += [-1]
                            else:
                                raise ValueError(f'The crossing sign was not + or - {cross}')                           
                    
                    # check if either termini has duplicate crossings and empty out those crossings lists as we will not have any confidence
                    N_duplicate_flag = False
                    if len(N_crossings_resids) != len(set(N_crossings_resids)):
                        N_duplicate_flag = True
                    
                    C_duplicate_flag = False
                    if len(C_crossings_resids) != len(set(C_crossings_resids)):
                        C_duplicate_flag = True

                    # check if both tails have duplicates. if so discard whole ent
                    if C_duplicate_flag and N_duplicate_flag:
                        print(f'Both N and C terminus have duplicates: {N_crossings_resids} | {C_crossings_resids}')
                        continue 

                    # get the sum of the N and C terminal crossings
                    if len(N_crossings_signs) != 0:
                        sum_N_crossings_signs = np.sum(N_crossings_signs)
                    else: 
                        sum_N_crossings_signs = 9999

                    if len(C_crossings_signs) != 0:
                        sum_C_crossings_signs = np.sum(C_crossings_signs)
                    else:
                        sum_C_crossings_signs = 9999

                    print(gene, ijr, i, j, r, N_crossings, N_crossings_signs, sum_N_crossings_signs, C_crossings, C_crossings_signs, sum_C_crossings_signs)
    
                    ## No slipknots found case
                    if sum_N_crossings_signs != 0 and sum_C_crossings_signs != 0:
                        print(f'No SlipKnots found')
                        if N_duplicate_flag == False and C_duplicate_flag == False:
                            outlines += [line]
                        elif N_duplicate_flag == True and C_duplicate_flag == False:
                            if len(C_crossings) != 0:
                                print(f'Duplicates found in N terminus: {N_crossings_resids} and C terminius still has crossings')
                                new_ijr = f' ({i}, {j}, {C_crossings}) '
                                new_line = line.split('|')
                                new_line[1] = new_ijr
                                new_line = '|'.join(new_line)
                                print(f'New line: {new_line}')
                                outlines += [new_line]
                            else:
                                print(f'Duplicates found in N terminus: {N_crossings_resids} and C terminius has NO crossings')
                                continue

                        elif N_duplicate_flag == False and C_duplicate_flag == True:
                            if len(N_crossings) != 0:
                                print(f'Duplicates found in C terminus: {C_crossings_resids} and N terminius still has crossings')
                                new_ijr = f' ({i}, {j}, {N_crossings}) '
                                new_line = line.split('|')
                                new_line[1] = new_ijr
                                new_line = '|'.join(new_line)
                                print(f'New line: {new_line}')
                                outlines += [new_line]
                            else:
                                print(f'Duplicates found in C terminus: {C_crossings_resids} and N terminius has NO crossings')
                                continue

                    
                    ## N terminal slipknot found and there is still a C terminal ENT make a new ent line
                    elif sum_N_crossings_signs == 0 and sum_C_crossings_signs != 0:
                        print(f'SlipKnot FOUND in N terminus for {line} in {gene} {pdb} {chain}')
                        if len(C_crossings) == 0 or C_duplicate_flag == True:
                            print(f'No C terminal entanglement or duplicate crossings found so this whole ent will be discarded')
                            continue

                        C_crossings = np.asarray(C_crossings, dtype=str)
                        new_ijr = f' ({i}, {j}, {C_crossings}) '
                        new_line = line.split('|')
                        new_line[1] = new_ijr
                        new_line = '|'.join(new_line)
                        print(f'New line: {new_line}')
                        outlines += [new_line]
                        
                    ## C terminal slipknot found and there is still a N terminal ENT make a new ent line
                    elif sum_N_crossings_signs != 0 and sum_C_crossings_signs == 0:
                        print(f'SlipKnot FOUND in C terminus for {line} in {gene} {pdb} {chain}')
                        if len(N_crossings) == 0 or N_duplicate_flag == True:
                            print(f'No N terminal entanglement or duplicate crossings found so this whole ent will be discarded')
                            continue

                        N_crossings = np.asarray(N_crossings, dtype=str)
                        new_ijr = f' ({i}, {j}, {N_crossings}) '
                        new_line = line.split('|')
                        new_line[1] = new_ijr
                        new_line = '|'.join(new_line)
                        print(f'New line: {new_line}')
                        outlines += [new_line]
                        
                    elif sum_N_crossings_signs == 0 and sum_C_crossings_signs == 0:
                        # both are slipknots and entanglemet will be discarded
                        print(f'SlipKnot FOUND in both N and C terminus for {line} in {gene} {pdb} {chain}')
                    
                else:
                    #else if only one crossing was idetified then output the line as it cannot have a slipknot
                    outlines += [line]
           
            #if gene == 'P16615':
            #    quit()
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
    args = parser.parse_args()

    uent_files = args.uent_files
    outpath = args.outpath

    # Make outpath if it doesnt exists
    if not os.path.exists(outpath):
        os.makedirs(outpath)
        print(f'Made directory: {outpath}')

    # Initalize the FrequencyGenerator class object
    analysis = Analyzer(outpath, uent_files)
    print(f'analysis: {analysis}')

    # check entanglement files for nonmapped entanglements
    analysis.remove_slipknots()






start_time = time.time()
if __name__ == "__main__":
    main()
print(f'NORMAL TERMINATION: {time.time() - start_time}')
