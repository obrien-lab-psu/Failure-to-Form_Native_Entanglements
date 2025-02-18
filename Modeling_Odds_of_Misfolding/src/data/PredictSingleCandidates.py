import time, sys
import multiprocessing as mp
from scipy.stats import bootstrap
import logging
import argparse
import glob
import numpy as np
import pandas as pd
from sklearn.preprocessing import OneHotEncoder, LabelEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.model_selection import train_test_split, StratifiedKFold, KFold, cross_validate, GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, average_precision_score, f1_score, recall_score, precision_score, roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from scipy.spatial.distance import euclidean
import matplotlib.pyplot as plt
import os
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import poisson, binom, fisher_exact, chi2, norm, permutation_test, false_discovery_control
import scipy.stats as st
from matplotlib.ticker import MultipleLocator

#pd.set_option('display.max_rows', 4000)

class DataAnalysis:
    """
    A class to handle the data analysis process including encoding, regression, and statistical tests.
    """

    def __init__(self, resFeat_files, outpath, gene_list, tag, load_style, buffer, spa, cov):
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
        """
        self.resFeat_files = resFeat_files
        self.outpath = outpath
        self.gene_list = gene_list
        self.tag = tag
        self.load_style = load_style
        self.buffer = buffer
        self.spa = spa
        self.cov = cov
        self.data = {}
        #self.logger = self.setup_logging()
        #self.gene_list_files = glob.glob(self.gene_list)

        if not os.path.exists(f'{self.outpath}'):
            os.makedirs(f'{self.outpath}')
            print(f'Made output directories {self.outpath}')

    #########################################################################################################
    def setup_logging(self):
        """
        Sets up the logging configuration.

        Returns:
        - logger (logging.Logger): Configured logger.
        """
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)
        return logger
    #########################################################################################################

    #########################################################################################################
    def encode_boolean_columns(self, df: pd.DataFrame, boolean_columns: list) -> pd.DataFrame:
        """
        Encodes boolean-like columns in a DataFrame to binary 0 and 1.

        Parameters:
        - df (pd.DataFrame): The input DataFrame.
        - boolean_columns (list): A list of column names to be encoded.

        Returns:
        - pd.DataFrame: The DataFrame with encoded columns.
        """
        label_encoder = LabelEncoder()
        
        for column in boolean_columns:
            if column in df.columns:
                df[column] = label_encoder.fit_transform(df[column])
            else:
                print(f"Column '{column}' does not exist in the DataFrame.")
        
        return df
    #########################################################################################################

    #########################################################################################################
    def regression(self, df, formula):
        """
        Performs quasi-binomial regression analysis on the provided DataFrame.

        Parameters:
        - df (pd.DataFrame): DataFrame containing the data for regression.
        - formula (str): The formula specifying the regression model.

        Returns:
        - table_1_df (pd.DataFrame): DataFrame containing the regression results with p-values.
        """
        model = sm.GLM.from_formula(formula, family=sm.families.Binomial(), data=df)
        #model = smf.logit(formula=formula, data=df)
        result = model.fit()


        ## recalculate the pvalue to add more digits as statsmodels truncates it to 0 if it is below 0.0001 for some reason. 
        print(result.summary())
        table = result.summary().tables[1]
        table_df = pd.DataFrame(table.data[1:], columns=table.data[0])
        pvalues = []
        for z in table_df['z']:
            z = float(z)
            if z < 0:
                p = st.norm.cdf(z)
            else:
                p = 1 - st.norm.cdf(z)
            pvalues += [p*2]

        table_df['P>|z|'] = pvalues
        table_df = table_df.rename(columns={"": "var"})
        return table_df, odds, cov_matrix
    #########################################################################################################

    #########################################################################################################
    def load_data(self):
        """
        Loads the residue feature files and filters the data for analysis.
        """
        if self.load_style == 'True':
            res_files = glob.glob(self.resFeat_files)
            print(f'Number of res_files: {len(res_files)}')
            for i, gene in enumerate(self.reg_genes):
                gene_resFeat = [f for f in res_files if gene in f]
                if len(gene_resFeat) == 0:
                    print(f"No residue feature file found for gene {gene}")
                    continue
                elif len(gene_resFeat) > 1:
                    print(f"More than 1 residue feature file found for gene {gene}")
                    continue
                gene_resFeat_file = gene_resFeat[0]
                #print(f'gene_resFeat_file: {gene_resFeat_file} {i}')
                if len(self.data) == 0:
                    self.data = pd.read_csv(gene_resFeat_file, sep='|')
                else:
                    self.data = pd.concat((self.data, pd.read_csv(gene_resFeat_file, sep='|')))
        else:
            self.data = pd.read_csv(self.resFeat_files, sep='|')

        self.data = self.data[self.data['gene'].isin(self.reg_genes)]
        self.data = self.data[self.data['AA'] != 'NC']
        self.data = self.data[self.data['mapped_resid'].notna()]
        self.data = self.data[self.data['AA'].notna()]
        self.data = self.data.reset_index()
        self.data = self.data[self.data['gene'].isin(self.reg_genes)]
        keys = ['index', 'gene', 'pdb', 'chain', 'uniprot_length', 'essential', 'ent_present', 'pdb_resid', 'resname', 'AA', 'region', 'mapped_resid', 'cut_str', 'cut_C_Rall', 'cut_CD_Rall', 'cut_CG_Rall']
        self.data = self.data[keys]
        print(f"Data loaded and filtered. Number of unique genes: {len(self.data['gene'].unique())}")
        print(f'self.data\n{self.data}')

        ## for each protein determine proportion of entangled region cut and same for non-entangled region
        protein_level_df = {'gene':[], 'uniprot_length':[], 'buffer':[], 
                            'EntR_N':[], 'EntR_cut_N':[], 'FracEntR_cut':[], 
                            'NonEntR_N':[], 'NonEntR_cut_N':[], 'FracNonEntR_cut':[], 'statistic':[], 'pvalue':[],
                            'NumTimepointsEntRCut':[], 'NumTimepointsNonEntRCut':[], 
                            'Both5min2hrEntRCut':[], 'Both5min2hrNonEntRCut':[]}

        protein_level_outfile = os.path.join(self.outpath, 'protein_level_summary.csv')        
        if not os.path.exists(protein_level_outfile):
            for buff in ['C', 'CD', 'CG']:
                for gene, gene_df in self.data.groupby('gene'):
                    #print(gene_df)
                    uniprot_length = gene_df['uniprot_length'].values[0]

                    EntR_N = len(gene_df[gene_df['region'] == 1])
                    NonEntR_N = len(gene_df[gene_df['region'] == 0])


                    ## get the number of timepoitns the entangled region was cut
                    EntR_cut = gene_df[(gene_df[f'cut_{buff}_Rall'] == True) & (gene_df['region'] == 1)]
                    EntR_cut_N = len(EntR_cut)
                    #print(EntR_cut, EntR_cut_N)

                    EntR_cut_strs = EntR_cut['cut_str'].values
                    R1min, R5min, R2hr = False, False, False
                    for cut_str in EntR_cut_strs:
                        if f'{buff},R1min' in cut_str:
                            R1min = True
                        if f'{buff},R5min' in cut_str:
                            R5min = True
                        if f'{buff},R2hr' in cut_str:
                            R2hr = True
                    NumTimepointsEntRCut = np.sum([R1min, R5min, R2hr])

                    if R5min and R2hr:
                        Both5min2hrEntRCut = True
                    else:
                        Both5min2hrEntRCut = False

                    ## get the number of timepoitns the non-entangled region was cut
                    NonEntR_cut = gene_df[(gene_df[f'cut_{buff}_Rall'] == True) & (gene_df['region'] == 0)]
                    NonEntR_cut_N = len(NonEntR_cut)
                    #print(NonEntR_cut, NonEntR_cut_N)

                    NonEntR_cut_strs = NonEntR_cut['cut_str'].values
                    R1min, R5min, R2hr = False, False, False
                    for cut_str in NonEntR_cut_strs:
                        if f'{buff},R1min' in cut_str:
                            R1min = True
                        if f'{buff},R5min' in cut_str:
                            R5min = True
                        if f'{buff},R2hr' in cut_str:
                            R2hr = True
                    NumTimepointsNonEntRCut = np.sum([R1min, R5min, R2hr])

                    if R5min and R2hr:
                        Both5min2hrNonEntRCut = True
                    else:
                        Both5min2hrNonEntRCut = False

                    ## calculate the permutation value for this gene for the probability of observing ratio=(EntR_cut_N/EntR_N)/(NonEntR_cut_N/NonEntR_N) greater than 1 by random chance
                    EntR_arr = [1]*EntR_cut_N + [0]*EntR_N
                    NonEntR_arr = [1]*NonEntR_cut_N + [0]*NonEntR_N
                    res = permutation_test((EntR_arr, NonEntR_arr), statistic, vectorized=True, n_resamples=10000, alternative='greater')

                    protein_level_df['gene'] += [gene]
                    protein_level_df['uniprot_length'] += [uniprot_length]
                    protein_level_df['buffer'] += [buff]
                    protein_level_df['EntR_N'] += [EntR_N]
                    protein_level_df['EntR_cut_N'] += [EntR_cut_N]                                
                    protein_level_df['FracEntR_cut'] += [EntR_cut_N/EntR_N]
                    protein_level_df['NonEntR_N'] += [NonEntR_N]
                    protein_level_df['NonEntR_cut_N'] += [NonEntR_cut_N]
                    protein_level_df['FracNonEntR_cut'] += [NonEntR_cut_N/NonEntR_N]                
                    protein_level_df['statistic'] += [res.statistic] 
                    protein_level_df['pvalue'] += [res.pvalue] 
                    protein_level_df['NumTimepointsEntRCut'] += [NumTimepointsEntRCut]
                    protein_level_df['NumTimepointsNonEntRCut'] += [NumTimepointsNonEntRCut]
                    protein_level_df['Both5min2hrEntRCut'] += [Both5min2hrEntRCut]
                    protein_level_df['Both5min2hrNonEntRCut'] += [Both5min2hrNonEntRCut]                   

            protein_level_df = pd.DataFrame(protein_level_df)
            all_qvalues = []
            for buff, buff_df in protein_level_df.groupby('buffer'):
                pvalues = buff_df['pvalue'].values
                qvalues = false_discovery_control(pvalues)
                all_qvalues += [qvalues]
            protein_level_df['qvalues'] = np.hstack(all_qvalues)

            protein_level_df.to_csv(protein_level_outfile, index=False)
            print(f'SAVED: {protein_level_outfile}')
        else:
            protein_level_df = pd.read_csv(protein_level_outfile)
            print(f'LOADED: {protein_level_outfile}')
        print(f'protein_level_df:\n{protein_level_df}')
        print(protein_level_df[protein_level_df['gene'] == 'P0AD61'])
        self.protein_level_df = protein_level_df
    #########################################################################################################

    #########################################################################################################
    def filter(self,):
        """
        Filter proteins and rank order them based on 
        1. avg Fraction of entangled region cut across all three buffer conditions (max = 1)
        2. number of buffer conditions with entangled region cut (max = 3)
        3. number of timepoints cut (max = 9)
        """

        filtered_df = {'gene':[], 'avgFractEntR_cut':[], 'avgFractNonEntR_cut':[], 'avgDiff':[], 'NumBuffEntR_cut':[]}

        ## inital filter for those that have cuts at 5min and 2hr in entR atleast
        df = self.protein_level_df[self.protein_level_df['Both5min2hrEntRCut'] == True]
        print(df)

        for gene, gene_df in df.groupby('gene'):
            #print(gene_df)

            avgFractEntR_cut = np.mean(gene_df['FracEntR_cut'].values)
            avgFractNonEntR_cut = np.mean(gene_df['FracNonEntR_cut'].values)
            avgDiff = np.mean(gene_df['statistic'].values)
            NumBuffEntR_cut = len(gene_df[gene_df['EntR_cut_N'] != 0])
            #Score = avgFractEntR_cut + (NumBuffEntR_cut/3) + (NumTimepointsEntR_cut/9)

            filtered_df['gene'] += [gene]
            filtered_df['avgFractEntR_cut'] += [avgFractEntR_cut]
            filtered_df['avgFractNonEntR_cut'] += [avgFractNonEntR_cut]
            filtered_df['avgDiff'] += [avgDiff]
            filtered_df['NumBuffEntR_cut'] += [NumBuffEntR_cut]
            #filtered_df['Score'] = Score

        filtered_df = pd.DataFrame(filtered_df)

        # Sort by 'Score' in descending order
        sorted_df = filtered_df.sort_values(by=['NumBuffEntR_cut', 'avgFractEntR_cut'], ascending=False, ignore_index=True)
        print(f'sorted_df:\n{sorted_df.to_string()}')
        sorted_df_outfile = os.path.join(self.outpath, 'filtered_protein_level_summary.csv')
        sorted_df.to_csv(sorted_df_outfile, index=False)
        print(f'SAVED: {sorted_df_outfile}')

    #########################################################################################################

    #########################################################################################################
    def score(self,):
        """
        Filter proteins and rank order them based on 
        1. avg Fraction of entangled region cut across all three buffer conditions (max = 1)
        2. number of buffer conditions with entangled region cut (max = 3)
        3. number of timepoints cut (max = 9)
        """

        #filtered_df = {'gene':[], 'avgFractEntR_cut':[], 'FractBuffEntR_cut':[], 'FractTimepointsEntR_cut':[]}
        filtered_df = {'gene':[], 'avgDiffRegionFractCut':[], 'FractBuffEntR_cut':[], 'FractTimepointsEntR_cut':[]}

        for gene, gene_df in self.protein_level_df.groupby('gene'):
            #print(gene_df)

            #avgFractEntR_cut = np.mean(gene_df['FracEntR_cut'].values)
            avgDiffRegionFractCut = np.mean(gene_df['statistic'].values)
            NumBuffEntR_cut = len(gene_df[gene_df['EntR_cut_N'] != 0])
            NumTimepointsEntR_cut = np.sum(gene_df['NumTimepointsEntRCut'].values)
            #Score = avgFractEntR_cut + (NumBuffEntR_cut/3) + (NumTimepointsEntR_cut/9)

            filtered_df['gene'] += [gene]
            #filtered_df['avgFractEntR_cut'] += [avgFractEntR_cut]
            filtered_df['avgDiffRegionFractCut'] += [avgDiffRegionFractCut]
            filtered_df['FractBuffEntR_cut'] += [NumBuffEntR_cut/3]
            filtered_df['FractTimepointsEntR_cut'] += [NumTimepointsEntR_cut/9]

        filtered_df = pd.DataFrame(filtered_df)
        print(f'filtered_df:\n{filtered_df}')

        ## standardize the dataframe
        scaler = StandardScaler()
        # Selecting the columns to be standardized
        #cols_to_standardize = ['avgFractEntR_cut', 'FractBuffEntR_cut', 'FractTimepointsEntR_cut']
        cols_to_standardize = ['avgDiffRegionFractCut', 'FractBuffEntR_cut', 'FractTimepointsEntR_cut']
        
        # Standardizing the selected columns
        filtered_df[cols_to_standardize] = scaler.fit_transform(filtered_df[cols_to_standardize])   
        print(f'filtered_df:\n{filtered_df}')  

        #filtered_df['Score'] = filtered_df['avgFractEntR_cut'] + filtered_df['FractBuffEntR_cut'] + filtered_df['FractTimepointsEntR_cut']
        filtered_df['Score'] = filtered_df['avgDiffRegionFractCut'] + filtered_df['FractBuffEntR_cut'] + filtered_df['FractTimepointsEntR_cut']
        print(f'filtered_df:\n{filtered_df}')  

        sorted_df = filtered_df.sort_values(by=['Score'], ascending=False, ignore_index=True)
        print(sorted_df.to_string())
        quit()
        scores = sorted_df['Score'].values
        #print(f'scores: {scores} {scores.shape}')

        rng = np.random.default_rng()
        pvalues = np.zeros(len(scores))
        for p in range(10000):
            arr = sorted_df[cols_to_standardize].values
            #print(arr)
            p_arr = rng.permutation(arr)
            #print(p_arr)
            p_scores = np.sum(p_arr, axis=1)
            #print(f'p_scores: {p_scores} {p_scores.shape}')

            pvalues += np.where(p_scores >= scores, 1, 0)
        pvalues /= 10000
        qvalues = false_discovery_control(pvalues)
        sorted_df['pvalues'] = pvalues
        sorted_df['qvalues'] = qvalues
        print(sorted_df.to_string())
        quit()

        # Sort by 'Score' in descending order
        sorted_df = filtered_df.sort_values(by=['NumBuffEntR_cut', 'NumTimepointsEntR_cut', 'avgFractEntR_cut'], ascending=False, ignore_index=True)
        sorted_df_outfile = os.path.join(self.outpath, 'filtered_protein_level_summary.csv')
        sorted_df.to_csv(sorted_df_outfile, index=False)
        print(f'SAVED: {sorted_df_outfile}')
    #########################################################################################################

    #########################################################################################################
    def run(self):
        """
        Orchestrates the workflow by loading data, performing regression, and saving results.
        """
        start_time = time.time()

        # Get list of genes to select for in regression
        print(f'reg_gene_list: {self.gene_list}')
        reg_genes = np.loadtxt(self.gene_list, dtype=str)
        self.reg_genes = reg_genes
        try:
            num_genes = len(reg_genes)
        except:
            print(f'Failed to get number of genes in {self.gene_list}')
            quit()
        print(f'Number of genes in list: {len(reg_genes)}')

        if len(reg_genes) == 0:
            print(f"Gene list is empty: {self.gene_list} {len(reg_genes)}")
            quit()

        # Load data
        self.load_data()

        self.filter()

        #self.score()
        print('NORMAL TERMINATION')

def statistic(x, y, axis):
    return np.mean(x, axis=axis) - np.mean(y, axis=axis)

def main():
    """
    Main function to parse arguments and run the DataAnalysis class.
    """
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-f", "--resFeat_files", type=str, required=True, help="Path to residue feature files")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("-g", "--gene_list", type=str, required=True, help="Path to gene list to use")
    parser.add_argument("-t", "--tag", type=str, required=True, help="Tag for output filenames")
    parser.add_argument("-l", "--load_style", type=str, required=True, help="Load style (True: load by gene, False: load a single file with all genes present)")
    parser.add_argument("-b", "--buffer", type=str, required=True, help="Buffer system to use")
    parser.add_argument("-s", "--spa", type=str, required=True, help="SPA threshold")
    parser.add_argument("-c", "--cov", type=str, required=True, help="LiPMS cov threshold")

    args = parser.parse_args()

    analysis = DataAnalysis(
        resFeat_files=args.resFeat_files,
        outpath=args.outpath,
        gene_list=args.gene_list,
        tag=args.tag,
        load_style=args.load_style,
        buffer=args.buffer,
        spa=args.spa,
        cov=args.cov)
    analysis.run()

if __name__ == "__main__":
    main()

