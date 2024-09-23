import sys,os
import ast
import math
import logging
import argparse
import re
import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import permutation_test, bootstrap
import glob
# Set display options to prevent summarization
pd.set_option('display.max_rows', None)
#pd.set_option('display.max_columns', None)

class MotifScanner:
    def __init__(self, outpath, fasta_dir, rep_genes_file, Ess_genes_file, NonEss_genes_file, clustered_ent_dir, buff, spa):
        """
        Initializes the MotifScanner class with necessary paths and parameters.
        "--outpath", type=str, required=True, help="Path to output directory"
        "--clustered_ent_dir", type=str, required=True, help="Tag for output filenames"
        "--rep_genes_files", type=str, required=True, help="Load style (True: load by gene, False: load a single file with all genes present)"
        "--Ess_genes_file", type=str, required=True, help="Path to the Essential genes list"
        "--NonEss_genes_file", type=str, required=True, help="Path to the NonEssential genes list"
        """
        self.outpath = outpath
        self.buff = buff
        self.spa = spa

        self.fasta_dir = fasta_dir
        self.fasta_files = glob.glob(f'{self.fasta_dir}*') 
        self.clustered_ent_files = glob.glob(f'{clustered_ent_dir}*')

        self.rep_genes = np.loadtxt(rep_genes_file, dtype=str)
        self.Ess_genes = np.loadtxt(Ess_genes_file, dtype=str)
        self.NonEss_genes = np.loadtxt(NonEss_genes_file, dtype=str)
        #self.logger = self.setup_logging(log_file)

        if not os.path.exists(f'{self.outpath}'):
            os.makedirs(f'{self.outpath}')
            print(f'Made output directories {self.outpath}')

        self.outpathGeneScans = f'{self.outpath}GeneScans/'
        if not os.path.exists(self.outpathGeneScans):
            os.makedirs(self.outpathGeneScans)
            print(f'Made output directories {self.outpathGeneScans}')

        ## Define motif patterns to search 
        self.Bukau4 = {0: ['R', 'K'],
                1: ['I', 'L', 'V', 'F', 'Y'],
                2: ['I', 'L', 'V', 'F', 'Y'],
                3: ['I', 'L', 'V', 'F', 'Y'],
                4: ['I', 'L', 'V', 'F', 'Y'],
                5: ['R', 'K']}
        
        self.Bukau4RF = {0: ['I', 'L', 'V', 'F', 'Y'],
                1: ['I', 'L', 'V', 'F', 'Y'],
                2: ['I', 'L', 'V', 'F', 'Y'],
                3: ['I', 'L', 'V', 'F', 'Y'],
                4: ['R', 'K']}

        self.Bukau4LF = {0: ['R', 'K'],
                1: ['I', 'L', 'V', 'F', 'Y'],
                2: ['I', 'L', 'V', 'F', 'Y'],
                3: ['I', 'L', 'V', 'F', 'Y'],
                4: ['I', 'L', 'V', 'F', 'Y']}

        self.Bukau5 = {0: ['R', 'K'],
                1: ['I', 'L', 'V', 'F', 'Y'],
                2: ['I', 'L', 'V', 'F', 'Y'],
                3: ['I', 'L', 'V', 'F', 'Y'],
                4: ['I', 'L', 'V', 'F', 'Y'],
                5: ['I', 'L', 'V', 'F', 'Y'],
                6: ['R', 'K']}

        self.Bukau5RF = {0: ['I', 'L', 'V', 'F', 'Y'],
                1: ['I', 'L', 'V', 'F', 'Y'],
                2: ['I', 'L', 'V', 'F', 'Y'],
                3: ['I', 'L', 'V', 'F', 'Y'],
                4: ['I', 'L', 'V', 'F', 'Y'],
                5: ['R', 'K']}

        self.Bukau5LF = {0: ['R', 'K'],
                1: ['I', 'L', 'V', 'F', 'Y'],
                2: ['I', 'L', 'V', 'F', 'Y'],
                3: ['I', 'L', 'V', 'F', 'Y'],
                4: ['I', 'L', 'V', 'F', 'Y'],
                5: ['I', 'L', 'V', 'F', 'Y']}

        self.Schymkowitz = {0: ['E', 'K', 'Q', 'W', 'Y'],
                1: ['F', 'I', 'K', 'L', 'R', 'V', 'Y'],
                2: ['F', 'L', 'R', 'V', 'W', 'Y'],
                3: ['I', 'L', 'M', 'T', 'V'],
                4: ['F', 'L', 'M', 'P', 'R', 'V', 'Y'],
                5: ['F', 'I', 'L', 'W', 'Y'],
                6: ['F', 'N', 'P', 'R', 'Y']}

        self.Emperical4 = {0: ['R', 'K', 'H'],
                1: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                2: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                3: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                4: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                5: ['R', 'K', 'H']}

        self.Emperical4RF = {0: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                1: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                2: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                3: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                4: ['R', 'K', 'H']}

        self.Emperical4LF = {0: ['R', 'K', 'H'],
                1: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                2: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                3: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                4: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y']}

        self.Emperical5 = {0: ['R', 'K', 'H'],
                1: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                2: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                3: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                4: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                5: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                6: ['R', 'K', 'H']}

        self.Emperical5RF = {0: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                1: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                2: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                3: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                4: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                5: ['R', 'K', 'H']}

        self.Emperical5LF = {0: ['R', 'K', 'H'],
                1: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                2: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                3: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                4: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'],
                5: ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y']}

    def setup_logging(self, log_file):
        """
        Sets up the logging configuration.

        Returns:
        - logger (logging.Logger): Configured logger.
        """
        # Setup logging configuration
        logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s %(message)s') 
        logger = logging.getLogger(__name__)
        logger.info(f'{"#"*15}NewRun{"#"*15}')
        return logger


    def GetResiduePools(self, ):
        """
        For each gene in the pool of Essential and NonEssential genes get the pool of residues to search for the motifs in.
        'All' -- all residues in the protein
        'threads' -- threads of unique ents +/- 10 residues 
        'loops' -- minimal loops of the unique entanglements  

        residue feature files being loaded have the following columns
        gene
        pdb
        chain
        essential
        pdb_resid
        resname
        AA
        nearest_neighbors
        num_nearest_neighbors
        region
        ent_idx
        res_sasa
        median_sasa
        NC
        crossing
        mapped_resid
        secondary_struct
        SCOP_class
        IDR
        cut_str
        cut_C_Rall
        cut_CD_Rall
        cut_CG_Rall
        buried

        """

        pool = {}
        for gene, pdb, chain in self.rep_genes:

            print(f'\n{"#"*30}\ngene: {gene} | pdb: {pdb} | chain: {chain}')
            
            #if gene != 'P0A7G6':
            #    continue

            # get FASTA seq
            gene_fasta_file = [f for f in self.fasta_files if gene in f]
            #print(f'gene: {gene} | gene_fasta_file: {gene_fasta_file}')
            Fasta_seq = [line.strip('\n') for line in open(gene_fasta_file[0]).readlines() if '>' not in line]
            Fasta_seq = ''.join(Fasta_seq)
            print(f'Fasta_seq: {Fasta_seq}')
            prot_length = len(Fasta_seq)
            print(f'prot_length: {prot_length}')
            #total_seq = ''.join(filled_df['AA'].values)
            #print(f'Total_seq: {total_seq}')
  
            # get the mapping file for mapping PDB resid to uniprot resid
            # User can change this dictionary if they require custome mapping but in this case 
            mapping_pdb2uniprot = {i+1:i+1 for i,x in enumerate(Fasta_seq)}
            #print(f'mapping_pdb2uniprot: {mapping_pdb2uniprot}')

            # get the clustered ent file if it exists
            gene_cent_file = [f for f in self.clustered_ent_files if gene in f]
            print(f'gene_cent_file: {gene_cent_file}')

            # loop through entangled residues and get residues within +/- 10 of each crossing and the loop residues
            if gene_cent_file:
                ent_present = True

                # map clustered entangled residues to uniprot mapping so i can slice fasta file
                cent_ijr = [ast.literal_eval(ijr) for ijr in pd.read_csv(gene_cent_file[0], sep='|')['ijr'].values]
                print(f'cent_ijr: {cent_ijr}')
                
                threads_seq = []
                threads_seq_start = []
                loops_seq = []
                loops_seq_start = []

                for idx, ijr in enumerate(cent_ijr):
                    print(f'\n{"#"*30}\nidx: {idx} | ijr: {ijr}')

                    # get mapping for start of loop
                    i = ijr[0]
                    if i in mapping_pdb2uniprot:
                        i = mapping_pdb2uniprot[i]
                    else:
                        print(f'PDB residue i={i} is not mappable')
                        continue

                    # get mapping for end of loop
                    j = ijr[1]
                    if j in mapping_pdb2uniprot:
                        j = mapping_pdb2uniprot[j]
                    else:
                        print(f'PDB residue j={j} is not mappable')
                        continue

                    # get mapping for crossings
                    r = [abs(ast.literal_eval(c)) for c in ijr[2:]]
                    crossings = []
                    unmappable_crossing = False
                    for cross in r:
                        if cross in mapping_pdb2uniprot:
                            crossings += [mapping_pdb2uniprot[cross]]
                        else:
                            print(f'PDB residue r={cross} is not mappable')
                            unmappable_crossing = True
                            break
                    if unmappable_crossing == True:
                        continue

                    mapped_ijr = (i, j, *crossings)
                    print(f'mapped_ijr: {mapped_ijr}')
                    for cross in mapped_ijr[2:]:
                        expanded_cross = np.arange(cross - 10, cross + 11)
                        expanded_cross = expanded_cross[expanded_cross >= 1] 
                        expanded_cross = expanded_cross[expanded_cross <= prot_length] 
                        print(f'cross: {cross} | expanded_cross: {expanded_cross}')
                        thread_seq = Fasta_seq[min(expanded_cross) - 1 : max(expanded_cross)]
                        print(f'thread_seq: {thread_seq}')
                        threads_seq += [thread_seq]
                        threads_seq_start += [min(expanded_cross)]

                    #loop = np.arange(i, j + 1)
                    loop_seq = Fasta_seq[i - 1: j]
                    loops_seq_start += [i]
                    loops_seq += [loop_seq]
                    print(f'loop_seq: {loop_seq}')

            else:
                threads_seq = ['']
                threads_seq_start = ['']
                loops_seq  = ['']
                loops_seq_start = ['']
                ent_present = False

            # add information to output dictionary    
            pool[gene] = {'All':Fasta_seq, 'Threads':threads_seq, 'Threads_start':threads_seq_start, 'Loops':loops_seq, 'Loops_start':loops_seq_start, 'ent_present':ent_present}

        self.pool = pool
        

    def CalcD(self, motif):
        """"
        for each gene calculate D(X) = Number of patterns found / protein_size
        where X = {'All':total_seq, 'Threads':threads_seq, 'Loops':loops_seq}
        ensure there is no double counting of motifs
        """

        if motif == 'Bukau4':
            pattern_dict = self.Bukau4
        elif motif == 'Bukau4RF':
            pattern_dict = self.Bukau4RF
        elif motif == 'Bukau4LF':
            pattern_dict = self.Bukau4LF

        elif motif == 'Bukau5':
            pattern_dict = self.Bukau5
        elif motif == 'Bukau5RF':
            pattern_dict = self.Bukau5RF
        elif motif == 'Bukau5LF':
            pattern_dict = self.Bukau5LF

        elif motif == 'Schymkowitz':
            pattern_dict = self.Schymkowitz

        elif motif == 'Emperical4':
            pattern_dict = self.Emperical4
        elif motif == 'Emperical4RF':
            pattern_dict = self.Emperical4RF
        elif motif == 'Emperical4LF':
            pattern_dict = self.Emperical4LF

        elif motif == 'Emperical5':
            pattern_dict = self.Emperical5
        elif motif == 'Emperical5RF':
            pattern_dict = self.Emperical5RF
        elif motif == 'Emperical5LF':
            pattern_dict = self.Emperical5LF
        print(f'CalcD with {motif} motif')
        print(f'pattern_dict: {pattern_dict}')

        ### scane sequence for pattern
        gene_scan = {'gene': [], 'motif':[], 'num_matches_All':[], 'D_All':[], 'num_matches_Thread':[], 'D_Thread':[], 'num_matches_Loop':[], 'D_Loop':[], 'ent_present':[]}

        for gene, seq_data in self.pool.items():
            
            self.gene = gene

            print(f'{"#"*30}\nGene: {gene} {motif} search')

            ###############################
            # calculate D_All
            L = len(seq_data['All'])
            seq = seq_data['All']
            #print(f'seq: {seq} {L}')
            ent_present = seq_data['ent_present']
                
            ### Scan sequence
            All_match_found = self.scan_seq(seq, pattern_dict, 'All', 0)
            num_All_match_found = len(All_match_found)
            print(f'All_match_found: {All_match_found} {num_All_match_found}\n')
            D_All = num_All_match_found/L
            
            #if gene == 'P0AFI7':
            #    quit()

            if D_All == 0:
                gene_scan['gene'] += [gene]
                gene_scan['motif'] += [motif]
                gene_scan['num_matches_All'] += [0]
                gene_scan['D_All'] += [0.0]
                gene_scan['num_matches_Thread'] += [0]
                gene_scan['D_Thread'] += [0.0]
                gene_scan['num_matches_Loop'] += [0]
                gene_scan['D_Loop'] += [0.0]
                gene_scan['ent_present'] += [ent_present]
                continue

            ###############################
            # calculate D_Thread
            Thread_match_found = set()
            for seq_i, seq in enumerate(seq_data['Threads']):
                start = seq_data['Threads_start'][seq_i]
                Thread_match_found = Thread_match_found.union(self.scan_seq(seq, pattern_dict, 'Thread', start))
                #print(f'Thread_match_found: {Thread_match_found} in seq: {seq} {seq_i}')
            num_Thread_match_found = len(Thread_match_found)
            print(f'ALL Thread_match_found: {Thread_match_found} {num_Thread_match_found}\n')
            D_Thread = num_Thread_match_found/L
            

            ###############################
            # calculate D_Loop
            Loop_match_found = set()
            for seq_i, seq in enumerate(seq_data['Loops']):
                start = seq_data['Loops_start'][seq_i]
                Loop_match_found = Loop_match_found.union(self.scan_seq(seq, pattern_dict, 'Loop', start))
                #print(f'Loop_match_found: {Loop_match_found} in seq: {seq} {seq_i}')
            num_Loop_match_found = len(Loop_match_found)
            print(f'ALL Loop_match_found: {Loop_match_found} {num_Loop_match_found}\n')
            D_Loop = num_Loop_match_found/L

            ##############################
            # update output dataframe
            gene_scan['gene'] += [gene]
            gene_scan['motif'] += [motif]
            gene_scan['num_matches_All'] += [num_All_match_found]
            gene_scan['D_All'] += [D_All]
            gene_scan['num_matches_Thread'] += [num_Thread_match_found]
            gene_scan['D_Thread'] += [D_Thread]
            gene_scan['num_matches_Loop'] += [num_Loop_match_found]
            gene_scan['D_Loop'] += [D_Loop]
            gene_scan['ent_present'] += [ent_present]

        gene_scan = pd.DataFrame(gene_scan)
        print(f'{motif} gene_scan:\n{gene_scan}')

        gene_scan_outfile = f'{self.outpathGeneScans}{motif}_scan_{self.buff}_{self.spa}.csv'
        gene_scan.to_csv(gene_scan_outfile, sep='|', index=False)
        print(f'SAVED: {gene_scan_outfile}')
        return gene_scan
        
    def scan_seq(self, seq, pattern_dict, tag, offset):
        print(f'{"#"*30}\nGene: {self.gene}\nseq: {seq}')
        pattern_length = len(pattern_dict)
        seq_length = len(seq)
        match_found = 0
        unique_matches = []
        for start in range(seq_length - pattern_length):
            end = start + pattern_length
            sub_seq = seq[start:end]
 
            ## determine if match to pattern
            matchs = []
            for i,p in enumerate(sub_seq):
                if p in pattern_dict[i]:
                    matchs += [True]
                else:
                    matchs += [False]
                
            #match = sum(matchs)
            #match = match == pattern_length
            match = all(matchs)
            if match == True:
                match_found += 1
                unique_matches += [(start + offset, end + offset, sub_seq)]
            print(tag, self.gene, start + offset, end + offset, sub_seq, matchs, match, match_found)
        unique_matches = set(unique_matches)
        print(f'Matches found in {tag} seq: {seq} {unique_matches}')
        return unique_matches

    def Dstats(self, gene_scan, motif='Emperical', OnlyEnt=False):
        """
        Using the gene scane of all representative PBDs now estiamte the ratio of <D_X(Ess)>/<D_X(NonEss)>
        and get a pvalue for permutation differences
        """

        Ess_gene_scan = gene_scan[gene_scan['gene'].isin(self.Ess_genes)]
        NonEss_gene_scan = gene_scan[gene_scan['gene'].isin(self.NonEss_genes)]
        if OnlyEnt == True:
            Ess_gene_scan = Ess_gene_scan[Ess_gene_scan['ent_present'] == True]
            NonEss_gene_scan = NonEss_gene_scan[NonEss_gene_scan['ent_present'] == True]
        #print(f'Ess_gene_scan:\n{Ess_gene_scan}')
        #print(f'NonEss_gene_scan:\n{NonEss_gene_scan}')

        Dstats_results = {'motif':[], 
                            'D_type':[], 
                            'Ess_n':[], 
                            'NonEss_n':[], 
                            'Ess_hits':[], 
                            'NonEss_hits':[], 
                            'Ess_mean':[], 
                            'Ess_lower_ci':[], 
                            'Ess_upper_ci':[], 
                            'NonEss_mean':[], 
                            'NonEss_lower_ci':[], 
                            'NonEss_upper_ci':[], 
                            'ratio':[], 
                            'proportion':[], 
                            'abs_diff':[], 
                            'p_value':[], 
                            'OnlyEnt':[]}

        columns_to_test = ['D_All', 'D_Thread', 'D_Loop']
        for col in columns_to_test:
            ess_col = Ess_gene_scan[col].values
            num_ess_hits = len([x for x in ess_col if x != 0])
            mean_ess = np.mean(ess_col)
            noness_col = NonEss_gene_scan[col].values
            num_noness_hits = len([x for x in noness_col if x != 0])
            mean_noness = np.mean(noness_col)

            ## bootstrap confidence intervals
            ess_res = bootstrap((ess_col,), np.mean, vectorized=True)
            ess_ci = (ess_res.confidence_interval.low, ess_res.confidence_interval.high)
            noness_res = bootstrap((noness_col,), np.mean, vectorized=True)
            noness_ci = (noness_res.confidence_interval.low, noness_res.confidence_interval.high)
            
            stat = abs(mean_ess - mean_noness)
            proportion = mean_ess / (mean_ess + mean_noness)
            ratio = mean_ess / mean_noness
            Dstats_results['motif'] += [motif]
            Dstats_results['D_type'] += [col]
            Dstats_results['Ess_n'] += [len(Ess_gene_scan)]
            Dstats_results['NonEss_n'] += [len(NonEss_gene_scan)]
            Dstats_results['Ess_hits'] += [num_ess_hits]
            Dstats_results['NonEss_hits'] += [num_noness_hits]
            Dstats_results['Ess_mean'] += [mean_ess]
            Dstats_results['Ess_lower_ci'] += [ess_res.confidence_interval.low]
            Dstats_results['Ess_upper_ci'] += [ess_res.confidence_interval.high]
            Dstats_results['NonEss_mean'] += [mean_noness]
            Dstats_results['NonEss_lower_ci'] += [noness_res.confidence_interval.low]
            Dstats_results['NonEss_upper_ci'] += [noness_res.confidence_interval.high]
            Dstats_results['ratio'] += [ratio]
            Dstats_results['proportion'] += [proportion]
            Dstats_results['abs_diff'] += [stat]

            # Permutation test
            combined = np.hstack([ess_col, noness_col])
            #print(f'len ess_col: {len(ess_col)} | len noness_col: {len(noness_col)} | combined: {len(combined)}')
            #observed_diff = mean_ess - mean_noness
            p_value = 0
            for _ in range(10000):
                np.random.shuffle(combined)
                perm_mean_ess = np.mean(combined[:len(Ess_gene_scan)])
                perm_mean_noness = np.mean(combined[len(Ess_gene_scan):])
                perm_proportion = perm_mean_ess / (perm_mean_ess + perm_mean_noness)
                perm_stat = abs(perm_mean_ess - perm_mean_noness)
                
                if perm_stat >= stat:
                    p_value += 1

            p_value /= 10000
            #print(f'p_value: {p_value}')
            Dstats_results['p_value'] += [p_value]
            Dstats_results['OnlyEnt'] += [OnlyEnt]

        Dstats_results = pd.DataFrame(Dstats_results)
        return Dstats_results

###############################################################
### MAIN ######################################################
###############################################################
pattern_usages = """
The pattern can be any length and the code must fit the following parameters
B = basic amino acid: R K 
H = hydrophobic amino acids: F L I Y V   

Other amino acids not recognized in motif as they are severly depleted in DnaK binders: D E W P N Q S T H C M G A 
"""

def main():
    """
    Main function to parse arguments and run the MotifScanner class.
    """
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("-f", "--fasta_dir", type=str, required=True, help="Path to gene lists to use")
    parser.add_argument("-c", "--clustered_ent_dir", type=str, required=True, help="Tag for output filenames")
    parser.add_argument("-r", "--rep_genes_file", type=str, required=True, help="Load style (True: load by gene, False: load a single file with all genes present)")
    parser.add_argument("-Eg", "--Ess_genes_file", type=str, required=True, help="Path to the Essential genes list")
    parser.add_argument("-NEg", "--NonEss_genes_file", type=str, required=True, help="Path to the NonEssential genes list")
    parser.add_argument("-b", "--buff", type=str, required=True, help="buffer: C CD CG")
    parser.add_argument("-s", "--spa", type=str, required=True, help="spa; 0 10 20 30 ... 90")
    args = parser.parse_args()

    motifscanner = MotifScanner(
        outpath=args.outpath,
        fasta_dir=args.fasta_dir,
        rep_genes_file=args.rep_genes_file,
        Ess_genes_file=args.Ess_genes_file,
        NonEss_genes_file=args.NonEss_genes_file, 
        clustered_ent_dir=args.clustered_ent_dir,
        buff=args.buff,
        spa=args.spa)

    # for each gene get the sets of residues for potential investigation
    # (1) all residues in the protein
    # (2) all residues in the threads of the protein 
    # (3) all residues in the loops of the protein 
    motifscanner.GetResiduePools()

    # for the Essential genes calculate D = #hits/L protein
    #['Bukau4', 'Schymkowitz', 'Emperical4', 'Emperical5']:
    Bukau4_gene_scan = motifscanner.CalcD('Bukau4')
    Bukau4RF_gene_scan = motifscanner.CalcD('Bukau4RF')
    Bukau4LF_gene_scan = motifscanner.CalcD('Bukau4LF')

    Bukau5_gene_scan = motifscanner.CalcD('Bukau5')
    Bukau5RF_gene_scan = motifscanner.CalcD('Bukau5RF')
    Bukau5LF_gene_scan = motifscanner.CalcD('Bukau5LF')

    Schymkowitz_gene_scan = motifscanner.CalcD('Schymkowitz')

    Emperical4_gene_scan = motifscanner.CalcD('Emperical4')
    Emperical4RF_gene_scan = motifscanner.CalcD('Emperical4RF')
    Emperical4LF_gene_scan = motifscanner.CalcD('Emperical4LF')

    Emperical5_gene_scan = motifscanner.CalcD('Emperical5')
    Emperical5RF_gene_scan = motifscanner.CalcD('Emperical5RF')
    Emperical5LF_gene_scan = motifscanner.CalcD('Emperical5LF')

    data = [motifscanner.Dstats(Bukau4_gene_scan, motif='Bukau4', OnlyEnt=False),
             motifscanner.Dstats(Bukau4RF_gene_scan, motif='Bukau4RF', OnlyEnt=False),
             motifscanner.Dstats(Bukau4LF_gene_scan, motif='Bukau4LF', OnlyEnt=False),
             motifscanner.Dstats(Bukau4_gene_scan, motif='Bukau4', OnlyEnt=True),
             motifscanner.Dstats(Bukau4RF_gene_scan, motif='Bukau4RF', OnlyEnt=True),
             motifscanner.Dstats(Bukau4LF_gene_scan, motif='Bukau4LF', OnlyEnt=True),
             motifscanner.Dstats(Bukau5_gene_scan, motif='Bukau5', OnlyEnt=False),
             motifscanner.Dstats(Bukau5RF_gene_scan, motif='Bukau5RF', OnlyEnt=False),
             motifscanner.Dstats(Bukau5LF_gene_scan, motif='Bukau5LF', OnlyEnt=False),
             motifscanner.Dstats(Bukau5_gene_scan, motif='Bukau5', OnlyEnt=True),
             motifscanner.Dstats(Bukau5RF_gene_scan, motif='Bukau5RF', OnlyEnt=True),
             motifscanner.Dstats(Bukau5LF_gene_scan, motif='Bukau5LF', OnlyEnt=True),
             motifscanner.Dstats(Schymkowitz_gene_scan, motif='Schymkowitz', OnlyEnt=False),
             motifscanner.Dstats(Schymkowitz_gene_scan, motif='Schymkowitz', OnlyEnt=True),
             motifscanner.Dstats(Emperical4_gene_scan, motif='Emperical4', OnlyEnt=False),
             motifscanner.Dstats(Emperical4RF_gene_scan, motif='Emperical4RF', OnlyEnt=False),
             motifscanner.Dstats(Emperical4LF_gene_scan, motif='Emperical4LF', OnlyEnt=False),
             motifscanner.Dstats(Emperical4_gene_scan, motif='Emperical4', OnlyEnt=True),
             motifscanner.Dstats(Emperical4RF_gene_scan, motif='Emperical4RF', OnlyEnt=True),
             motifscanner.Dstats(Emperical4LF_gene_scan, motif='Emperical4LF', OnlyEnt=True),
             motifscanner.Dstats(Emperical5_gene_scan, motif='Emperical5', OnlyEnt=False),
             motifscanner.Dstats(Emperical5RF_gene_scan, motif='Emperical5RF', OnlyEnt=False),
             motifscanner.Dstats(Emperical5LF_gene_scan, motif='Emperical5LF', OnlyEnt=False),
             motifscanner.Dstats(Emperical5_gene_scan, motif='Emperical5', OnlyEnt=True),
             motifscanner.Dstats(Emperical5RF_gene_scan, motif='Emperical5RF', OnlyEnt=True),
             motifscanner.Dstats(Emperical5LF_gene_scan, motif='Emperical5LF', OnlyEnt=True)]

    total_df = pd.concat(data)
    total_df['buff'] = args.buff
    total_df['spa'] = args.spa
    total_df['q_value'] = stats.false_discovery_control(total_df['p_value'].values)
    print(f'total_df:\n{total_df}')
    
    total_df_outfile = f'{motifscanner.outpath}{args.buff}_{args.spa}_DnaK_binding_motif_stats.csv'
    total_df.to_csv(total_df_outfile, sep='|', index=False)
    print(f'SAVED: {total_df_outfile}')

if __name__ == "__main__":
    main()

print('NORMAL TERMINATION')
