import time, sys
import multiprocessing as mp
from scipy.stats import bootstrap
import logging
import argparse
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import poisson, binom, fisher_exact, chi2, norm
import scipy.stats as st
from matplotlib.ticker import MultipleLocator

#pd.set_option('display.max_rows', 4000)

class DataAnalysis:
    """
    A class to handle the data analysis process including encoding, regression, and statistical tests.
    """

    def __init__(self, resFeat_files, outpath, tag, buffer):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
        - resFeat_files (str): Path to residue feature files.
        - outpath (str): Path to the output directory.
        - gene_lists (str): Path to gene lists to use.
        - tag (str): Tag for output filenames.
        - load_style (str): Load style (True: load by gene, False: load a single file with all genes present).
        - buffer (str): Buffer system to use.
        - spa (str): SPA threshold.
        - cov (str): LiPMS cov threshold.
        - reg_formula (str): Regression formula.
        """
        self.resFeat_files = resFeat_files
        self.outpath = outpath
        self.tag = tag
        self.data = {}
        self.buffer = buffer
        #self.logger = self.setup_logging()
        #self.gene_list_files = glob.glob(self.gene_list)

        if not os.path.exists(f'{self.outpath}'):
            os.makedirs(f'{self.outpath}')
            print(f'Made output directories {self.outpath}')



    ###################################################################
    def load_data(self, ):
        """
        Loads the residue feature files and filters the data for analysis.
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/EXP/EXP_0.6g_C_Rall_spa50_LiPMScov50_essential_ent_genes.txt
        """

        ent_genes = np.loadtxt(f'../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/EXP/EXP_0.6g_{self.buffer}_Rall_spa50_LiPMScov50_ent_genes.txt', dtype=str)

        res_files = glob.glob(self.resFeat_files)
        print(f'Number of res_files: {len(res_files)}')
        for i, f in enumerate(res_files):
            gene = f.split('/')[-1].split('_')[0]
            #print(f, gene)
            if gene not in ent_genes:
                continue
            if len(self.data) == 0:
                self.data = pd.read_csv(f, sep='|')
            else:
                self.data = pd.concat((self.data, pd.read_csv(f, sep='|')))

        self.data = self.data[self.data['AA'] != 'NC']
        self.data = self.data[self.data['mapped_resid'].notna()]
        self.data = self.data[self.data['AA'].notna()]
        keys = ['gene', 'uniprot_length', 'region', f'cut_{self.buffer}_Rall']
        self.data = self.data[keys]
        self.data = self.data.reset_index()
        print(f"Data loaded and filtered. Number of unique genes: {len(self.data['gene'].unique())}")
        print(self.data)
        print(self.data.keys())


        ## condense data into one row per unique gene that has the following columns
        ## 1. the gene name
        ## 2. the length of the protein
        ## 3. the number of cut sites in the cut_C_Rall that are also in the region = 1 case
        ## 4. the number of cut sites in the cut_C_Rall that are also in the region = 0 case
        ## 9. the total number of cut sites in the cut_C_Rall
        ## 12. the total number of region = 1 sites
        ## 13. the total number of region = 0 sites

        condensed_data = pd.DataFrame(columns=['gene', 'uniprot_length', f'cut_{self.buffer}_Rall_1', f'cut_{self.buffer}_Rall_0', 
                                               f'total_cut_{self.buffer}_Rall', 'total_region_1', 'total_region_0', 'total_res'])
        for gene in self.data['gene'].unique():
            gene_data = self.data[self.data['gene'] == gene]
            uniprot_length = gene_data['uniprot_length'].unique()[0]
            cut_Rall_1 = len(gene_data[(gene_data[f'cut_{self.buffer}_Rall'] == 1) & (gene_data['region'] == 1)])
            cut_Rall_0 = len(gene_data[(gene_data[f'cut_{self.buffer}_Rall'] == 1) & (gene_data['region'] == 0)])
            total_cut_Rall = len(gene_data[gene_data[f'cut_{self.buffer}_Rall'] == 1])
            total_region_1 = len(gene_data[gene_data['region'] == 1])
            total_region_0 = len(gene_data[gene_data['region'] == 0])
            total_res = len(gene_data)
            condensed_data = condensed_data._append({'gene': gene, 'uniprot_length': uniprot_length, f'cut_{self.buffer}_Rall_1': cut_Rall_1, f'cut_{self.buffer}_Rall_0': cut_Rall_0,
                                                    f'total_cut_{self.buffer}_Rall': total_cut_Rall, 'total_region_1': total_region_1, 'total_region_0': total_region_0, 'total_res': total_res}, ignore_index=True)    
        print(f'condensed_data:\n{condensed_data}')
        
        ## save the dataframe to the output directory
        condensed_data_outfile = os.path.join(self.outpath, f'{self.buffer}_condensed_CutDist_data.csv')
        condensed_data.to_csv(condensed_data_outfile, index=False)
        print(f"SAVED: {condensed_data_outfile}")
        return condensed_data
    ###################################################################

    ###################################################################
    def calc_avg_cuts(self, condensed_df):
        """
        Calculate the average number of cuts in the total_C_Rall, total_CD_Rall, and total_CG_Rall as a function of the total_res
        but bin the total_res into n bins with width 10. 
        Also calculate the 95% confidence intervals for the average number of cuts in the total_C_Rall, total_CD_Rall, and total_CG_Rall as a function of the total_res
        and if there is no data in a bin then the average number of cuts is 0 and the 95% confidence intervals are (0,0)
        """
        ## bin the total_res into n bins with width 10
        n = 100
        condensed_df = condensed_df[['gene', 'total_res', f'total_cut_{self.buffer}_Rall']]
        condensed_df['total_res_bin'] = pd.cut(condensed_df['total_res'], bins=np.arange(0, condensed_df['total_res'].max()+n, n), right=False)
        print(f'condensed_df:\n{condensed_df}')
        print(f'condensed_df keys:\n{condensed_df.keys()}')
        avg_cuts = pd.DataFrame(columns=['total_res_bin', f'avg_cut_{self.buffer}_Rall', f'cut_{self.buffer}_Rall_CI', 'n'])
        for total_res_bin in condensed_df['total_res_bin'].unique():
            print(f'total_res_bin: {total_res_bin}')
            total_res_bin_data = condensed_df[condensed_df['total_res_bin'] == total_res_bin]
            print(f'total_res_bin_data:\n{total_res_bin_data}')
            if len(total_res_bin_data) == 0:
                continue
            avg_cut_Rall = total_res_bin_data[f'total_cut_{self.buffer}_Rall'].mean()
            std_cut_Rall = total_res_bin_data[f'total_cut_{self.buffer}_Rall'].std()
            print(f'avg_cut_Rall: {avg_cut_Rall}')
            print(f'std_cut_Rall: {std_cut_Rall}')
            if avg_cut_Rall == 0 or std_cut_Rall == 0 or len(total_res_bin_data[f'total_cut_{self.buffer}_Rall'].values) == 1:
                cut_Rall_CI = (0, 0)
            else:
                cut_Rall_CI = bootstrap(data=(total_res_bin_data[f'total_cut_{self.buffer}_Rall'].values,), statistic=np.mean, n_resamples=10000).confidence_interval
            avg_cuts = avg_cuts._append({'total_res_bin': total_res_bin, f'avg_cut_{self.buffer}_Rall': avg_cut_Rall, f'cut_{self.buffer}_Rall_CI': cut_Rall_CI, 'n':len(total_res_bin_data)}, ignore_index=True)
        print(f'avg_cuts:\n{avg_cuts}')

        ## order the dataframe by the start of each bin
        avg_cuts['total_res_bin'] = avg_cuts['total_res_bin'].apply(lambda x: x.left)
        avg_cuts = avg_cuts.sort_values(by='total_res_bin')
        print(f'avg_cuts:\n{avg_cuts}')

        ## save the dataframe to the output directory
        avg_cuts_outfile = os.path.join(self.outpath, f'{self.buffer}_avg_cuts.csv')
        avg_cuts.to_csv(avg_cuts_outfile, index=False)
        print(f"SAVED: {avg_cuts_outfile}")

        return avg_cuts



    #################################################################
    def run(self):
        """
        Orchestrates the workflow by loading data, performing regression, and saving results.
        """
        start_time = time.time()


        # Load data
        condensed_df = self.load_data()

        # call function to calculate the average per protein number of cuts in the total_C_Rall, total_CD_Rall, annd total_CG_Rall as a function of the total_res
        calc_avg_cuts = self.calc_avg_cuts(condensed_df)

        print('NORMAL TERMINATION')

def main():
    """
    Main function to parse arguments and run the DataAnalysis class.
    """
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-f", "--resFeat_files", type=str, required=True, help="Path to residue feature files")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("-t", "--tag", type=str, required=True, help="Tag for output filenames")
    parser.add_argument("-b", "--buffer", type=str, required=True, help="Buffer system to use")

    args = parser.parse_args()

    analysis = DataAnalysis(
        resFeat_files=args.resFeat_files,
        outpath=args.outpath,
        tag=args.tag, 
        buffer=args.buffer)
    analysis.run()

if __name__ == "__main__":
    main()

