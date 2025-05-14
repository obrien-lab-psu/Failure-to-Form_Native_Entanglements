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
from scipy.stats import poisson, binom, fisher_exact, chi2, norm
import scipy.stats as st
from matplotlib.ticker import MultipleLocator
from itertools import permutations
#pd.set_option('display.max_rows', 4000)

class DataAnalysis:
    """
    A class to handle the data analysis process including encoding, regression, and statistical tests.
    """

    def __init__(self, resFeat_files, outpath, gene_list_path, tag):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
        - resFeat_files (str): Path to residue feature files.
        - outpath (str): Path to the output directory.
        - gene_list_path (str): Path to gene lists to use.
        - tag (str): Tag for output filenames.
        """
        self.resFeat_files = resFeat_files
        self.outpath = outpath
        self.gene_list_path = gene_list_path
        self.tag = tag

        if not os.path.exists(f'{self.outpath}'):
            os.makedirs(f'{self.outpath}')
            print(f'Made output directories {self.outpath}')

    def setup_logging(self):
        """
        Sets up the logging configuration.

        Returns:
        - logger (logging.Logger): Configured logger.
        """
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)
        return logger

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
                df.loc[:,column] = label_encoder.fit_transform(df[column])
            else:
                print(f"Column '{column}' does not exist in the DataFrame.")

        return df[boolean_columns + ['AA', 'gene']]


    def regression(self, data):
        """
        Performs quasi-binomial regression analysis on the provided DataFrame.

        Parameters:
        - df (pd.DataFrame): DataFrame containing the data for regression.
        - formula (str): The formula specifying the regression model.

        Returns:
        - table_1_df (pd.DataFrame): DataFrame containing the regression results with p-values.
        """
        outdf = {'buffer': [], 'essential': [], 'OR': [], 'OR_lb':[], 'OR_ub':[], 'pvalue': []}
        for buffer in data:
            for essential in data[buffer]:
                loc_df = data[buffer][essential]
                #print(f'buffer: {buffer} essential: {essential}')
                formula = f'cut_{buffer}_Rall ~ AA + region'
                model = sm.GLM.from_formula(formula, family=sm.families.Binomial(), data=loc_df)
                result = model.fit()

                ## recalculate the pvalue to add more digits as statsmodels truncates it to 0 if it is below 0.0001 for some reason. 
                table = result.summary().tables[1]
                table_df = pd.DataFrame(table.data[1:], columns=table.data[0])

                OR = np.exp(float(table_df[table_df[''] == 'region']['coef'].values[0]))
                OR_lb = np.exp(float(table_df[table_df[''] == 'region']['[0.025'].values[0]))
                OR_ub = np.exp(float(table_df[table_df[''] == 'region']['0.975]'].values[0]))
                pvalue = float(table_df[table_df[''] == 'region']['P>|z|'].values[0])
                z = float(table_df[table_df[''] == 'region']['z'].values[0])
                if z < 0:
                    pvalue = st.norm.cdf(z)*2
                else:
                    pvalue = (1 - st.norm.cdf(z))*2
                #print(f'OR: {OR} OR_lb: {OR_lb} OR_ub: {OR_ub} pvalue: {pvalue}')

                ## populate the outdf
                outdf['buffer'] += [buffer]
                outdf['essential'] += [essential]
                outdf['OR'] += [OR]
                outdf['OR_lb'] += [OR_lb]
                outdf['OR_ub'] += [OR_ub]
                outdf['pvalue'] += [pvalue]

        outdf = pd.DataFrame(outdf)
        return outdf

    def load_data(self):
        """
        Loads the residue feature files and filters the data for analysis.
        """
        gene_lists = glob.glob(self.gene_list_path + '/*.txt')
        print(f'Number of gene_lists: {len(gene_lists)}')

        ## Make a dictionary that contains the residue feature dfs for each gene list in self.gene_list_path
        res_files = glob.glob(self.resFeat_files + '/*.csv')
        print(f'Number of residue feature files: {len(res_files)}')
        if len(res_files) == 0:
            print(f"No residue feature files found in {self.resFeat_files}")
        
        self.data = {}
        for gene_list in gene_lists:

            genes = np.loadtxt(gene_list, dtype=str)
            print(f'\nNumber of genes in gene_list: {gene_list} {len(genes)}')

            ## determine the buffer and essentiality of the list
            buffer = gene_list.split('/')[-1].split('_')[2]
            essential = gene_list.split('/')[-1].split('_')[6]
            print(f'buffer: {buffer} essential: {essential}')
            if buffer not in self.data:
                self.data[buffer] = {}
            if essential not in self.data[buffer]:
                self.data[buffer][essential] = {}
            
            ## get the residue feature files for the genes in the gene list
            self.gene_list_resFeat = pd.DataFrame()
            for i, gene in enumerate(genes):
                gene_resFeat = [f for f in res_files if gene in f]
                if len(gene_resFeat) == 0:
                    print(f"No residue feature file found for gene {gene}")
                    quit()
                elif len(gene_resFeat) > 1:
                    print(f"More than 1 residue feature file found for gene {gene}")
                    quit()
                gene_resFeat_file = gene_resFeat[0]
                #print(f'gene_resFeat_file: {gene_resFeat_file} {i}')
                if len(self.gene_list_resFeat) == 0:
                    self.gene_list_resFeat = pd.read_csv(gene_resFeat_file, sep='|')
                else:
                    self.gene_list_resFeat = pd.concat((self.gene_list_resFeat, pd.read_csv(gene_resFeat_file, sep='|')))
            
            ## Remove non-canonical amino acids and unmapped residues
            self.gene_list_resFeat = self.gene_list_resFeat[self.gene_list_resFeat['AA'] != 'NC']
            self.gene_list_resFeat = self.gene_list_resFeat[self.gene_list_resFeat['mapped_resid'].notna()]
            self.gene_list_resFeat = self.gene_list_resFeat[self.gene_list_resFeat['AA'].notna()]
            self.gene_list_resFeat = self.gene_list_resFeat.reset_index()

            ## encode the columns for cut_{buff}_Rall, region, AA
            self.gene_list_resFeat = self.encode_boolean_columns(self.gene_list_resFeat, [f'cut_{buffer}_Rall', 'region'])
            self.data[buffer][essential] = self.gene_list_resFeat
            print(self.gene_list_resFeat)
            print(f'Number of genes in self.gene_list_resFeat: {len(self.gene_list_resFeat["gene"].unique())}')

    def SignifantRescue(self, data, GT_regression, essential='essential', num_permutes=10000):
        """
        Determine within the set of essential genes whether there was a significant decrease in the OR beyond random chance in the presece of chaperones
        """
        print(f'\nDetermining if there is a significant rescue for the set of {essential} genes')
        C_df = data['C'][essential].copy()
        C_df.rename(columns={'cut_C_Rall': 'cuts'}, inplace=True)
        C_df['label'] = C_df['gene']+'_C'
        num_C_genes = len(C_df['gene'].unique())

        CD_df = data['CD'][essential].copy()
        CD_df.rename(columns={'cut_CD_Rall': 'cuts'}, inplace=True)
        CD_df['label'] = CD_df['gene']+'_CD'
        num_CD_genes = len(CD_df['gene'].unique())

        CG_df = data['CG'][essential].copy()
        CG_df.rename(columns={'cut_CG_Rall': 'cuts'}, inplace=True)
        CG_df['label'] = CG_df['gene']+'_CG'
        num_CG_genes = len(CG_df['gene'].unique())

        print(f'C_df: {len(C_df)} {num_C_genes} | CD_df: {len(CD_df)} {num_CD_genes} | CG_df: {len(CG_df)} {num_CG_genes}')

        ## get the Ground truth differences between the C-CD buffers for the set of essenital results and the C-CG
        GT_regression = GT_regression[GT_regression['essential'] == essential]
        print(f'GT_regression:\n{GT_regression}')
        
        GT_C2CD = GT_regression.loc[GT_regression['buffer'] == 'C', 'OR'].values[0] - GT_regression.loc[GT_regression['buffer'] == 'CD', 'OR'].values[0]
        GT_C2CG = GT_regression.loc[GT_regression['buffer'] == 'C', 'OR'].values[0] - GT_regression.loc[GT_regression['buffer'] == 'CG', 'OR'].values[0]
        print(f'GT diffs: delta(C - CD): {GT_C2CD} | delta(C - CG): {GT_C2CG}')

        ## do random permutations and calcualte the probability that the difference is greater
        combined_CnCD = pd.concat((C_df, CD_df))
        #print(combined_CnCD)
        combined_CnCD_labels = combined_CnCD['label'].unique()
        #print(f'combined_CnCD_labels: {combined_CnCD_labels}')

        combined_CnCG = pd.concat((C_df, CG_df))
        #print(combined_CnCG)
        combined_CnCG_labels = combined_CnCG['label'].unique()
        #print(f'combined_CnCG_labels: {combined_CnCG_labels}')

        formula = f'cuts ~ AA + region'
        CnCD_permute_diffs = []
        CnCG_permute_diffs = []
        for p in range(num_permutes):

            ## permutation for the C - CD case
            #p_CnCD_labels = np.random.permutation(combined_CnCD_labels)
            #p_C = combined_CnCD[combined_CnCD['label'].isin(p_CnCD_labels[:num_C_genes])]
            #p_CD = combined_CnCD[combined_CnCD['label'].isin(p_CnCD_labels[num_C_genes:])]

            p_CnCD_df = combined_CnCD.sample(frac=1).reset_index(drop=True)
            p_C = p_CnCD_df.iloc[:len(C_df)]
            p_CD = p_CnCD_df.iloc[len(C_df):]
            
            model = sm.GLM.from_formula(formula, family=sm.families.Binomial(), data=p_C)
            result = model.fit()
            table = result.summary().tables[1]
            table_df = pd.DataFrame(table.data[1:], columns=table.data[0])
            C_OR = np.exp(float(table_df[table_df[''] == 'region']['coef'].values[0]))

            model = sm.GLM.from_formula(formula, family=sm.families.Binomial(), data=p_CD)
            result = model.fit()
            table = result.summary().tables[1]
            table_df = pd.DataFrame(table.data[1:], columns=table.data[0])
            CD_OR = np.exp(float(table_df[table_df[''] == 'region']['coef'].values[0]))         

            p_C2CD = C_OR - CD_OR
            #print(f'p: {p} | C_OR: {C_OR} | CD_OR: {CD_OR} | delta(C - CD): {p_C2CD}')
            CnCD_permute_diffs += [p_C2CD]


            ## permutation for the C - CG case
            #p_CnCG_labels = np.random.permutation(combined_CnCG_labels)
            #p_C = combined_CnCG[combined_CnCG['label'].isin(p_CnCG_labels[:num_C_genes])]
            #p_CG = combined_CnCG[combined_CnCG['label'].isin(p_CnCG_labels[num_C_genes:])]

            p_CnCG_df = combined_CnCG.sample(frac=1).reset_index(drop=True)
            p_C = p_CnCG_df.iloc[:len(C_df)]
            p_CG = p_CnCG_df.iloc[len(C_df):]
            
            model = sm.GLM.from_formula(formula, family=sm.families.Binomial(), data=p_C)
            result = model.fit()
            table = result.summary().tables[1]
            table_df = pd.DataFrame(table.data[1:], columns=table.data[0])
            C_OR = np.exp(float(table_df[table_df[''] == 'region']['coef'].values[0]))

            model = sm.GLM.from_formula(formula, family=sm.families.Binomial(), data=p_CG)
            result = model.fit()
            table = result.summary().tables[1]
            table_df = pd.DataFrame(table.data[1:], columns=table.data[0])
            CG_OR = np.exp(float(table_df[table_df[''] == 'region']['coef'].values[0]))

            p_C2CG = C_OR - CG_OR
            #print(f'p: {p} | C_OR: {C_OR} | CG_OR: {CG_OR} | delta(C - CG): {p_C2CG}')
            CnCG_permute_diffs += [p_C2CG]

        CnCD_permute_diffs = np.array(CnCD_permute_diffs)
        CnCG_permute_diffs = np.array(CnCG_permute_diffs)

        ## calculate the p-values
        p_CnCD = len(CnCD_permute_diffs[CnCD_permute_diffs > GT_C2CD])/num_permutes
        p_CnCG = len(CnCG_permute_diffs[CnCG_permute_diffs > GT_C2CG])/num_permutes
        #print(f'p_CnCD: {p_CnCD} | p_CnCG: {p_CnCG}')
        return p_CnCD, p_CnCG

    def run(self):
        """
        Orchestrates the workflow by loading data, performing regression, and saving results.
        """
        start_time = time.time()

        # Load data
        self.load_data()

        # Perform ground truth regressions 
        GT_regression = self.regression(self.data)
        print(f'GT_regression:\n{GT_regression}')
        GT_regression.to_csv(f'{self.outpath}/{self.tag}_GT_regression.csv', index=False)
        print(f'SAVED: {self.outpath}/{self.tag}_GT_regression.csv')

        # determine if there is a significant drop in the OR in the presence of chaperones for essential and nonessential proteins
        Ess_C2CD_pvalue, Ess_C2CG_pvalue = self.SignifantRescue(self.data, GT_regression, essential='essential')
        print(f'{self.tag}| Essential: C2CD_pvalue: {Ess_C2CD_pvalue} | C2CG_pvalue: {Ess_C2CG_pvalue}')

        NonEss_C2CD_pvalue, NonEss_C2CG_pvalue = self.SignifantRescue(self.data, GT_regression, essential='nonessential')
        print(f'{self.tag}| NonEssential: C2CD_pvalue: {NonEss_C2CD_pvalue} | C2CG_pvalue: {NonEss_C2CG_pvalue}')
        print(f'NORMAL TERMINATION {time.time() - start_time} seconds')

def main():
    """
    Main function to parse arguments and run the DataAnalysis class.
    """
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-f", "--resFeat_files", type=str, required=True, help="Path to residue feature files")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("-g", "--gene_list_path", type=str, required=True, help="Path to gene lists to use")
    parser.add_argument("-t", "--tag", type=str, required=True, help="Tag for output filenames")

    args = parser.parse_args()

    analysis = DataAnalysis(
        resFeat_files=args.resFeat_files,
        outpath=args.outpath,
        gene_list_path=args.gene_list_path,
        tag=args.tag)
    analysis.run()

if __name__ == "__main__":
    main()

