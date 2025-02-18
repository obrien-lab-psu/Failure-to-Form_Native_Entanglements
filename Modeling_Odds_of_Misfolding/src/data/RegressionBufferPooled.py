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

#pd.set_option('display.max_rows', 4000)

class DataAnalysis:
    """
    A class to handle the data analysis process including encoding, regression, and statistical tests.
    """

    def __init__(self, resFeat_files, outpath, tag):
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
        #self.logger = self.setup_logging()
        #self.gene_list_files = glob.glob(self.gene_list)

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


    ###################################################################
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
    ###################################################################

    ###################################################################
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
        return table_df
    ###################################################################


    ###################################################################
    def load_data(self, ):
        """
        Loads the residue feature files and filters the data for analysis.
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/EXP/EXP_0.6g_C_Rall_spa50_LiPMScov50_essential_ent_genes.txt
        """

        res_files = glob.glob(self.resFeat_files)
        print(f'Number of res_files: {len(res_files)}')
        for i, f in enumerate(res_files):
            gene = f.split('/')[-1].split('_')[0]
            #print(f, gene)

            if len(self.data) == 0:
                self.data = pd.read_csv(f, sep='|')
            else:
                self.data = pd.concat((self.data, pd.read_csv(f, sep='|')))

        self.data = self.data[self.data['AA'] != 'NC']
        self.data = self.data[self.data['mapped_resid'].notna()]
        self.data = self.data[self.data['AA'].notna()]
        keys = ['gene', 'region', 'AA', 'cut_C_Rall', 'cut_CD_Rall', 'cut_CG_Rall']
        self.data = self.data[keys]
        self.data = self.data.reset_index()
        print(f"Data loaded and filtered. Number of unique genes: {len(self.data['gene'].unique())}")
        print(self.data)
    ###################################################################

    ###################################################################
    def pool_data(self, essential):
        """
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/EXP/EXP_0.6g_C_Rall_spa50_LiPMScov50_essential_ent_genes.txt
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/EXP/EXP_0.6g_CD_Rall_spa50_LiPMScov50_essential_ent_genes.txt
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/EXP/EXP_0.6g_CG_Rall_spa50_LiPMScov50_essential_ent_genes.txt
        """
        all_genes = []
        #all_genes = np.loadtxt(f'../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/EXP/EXP_0.6g_C_Rall_spa50_LiPMScov50_{essential}_ent_genes.txt', dtype=str)
        dfs = []
        for buff in ['C', 'CD', 'CG']:
        #for buff in ['C']:
            genes = np.loadtxt(f'../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/EXP/EXP_0.6g_{buff}_Rall_spa50_LiPMScov50_{essential}_ent_genes.txt', dtype=str)
            print('Num genes:', essential, buff, genes, len(genes))
            all_genes += [genes]

            buff_df = self.data[self.data['gene'].isin(genes)]
            buff_df['cut'] = buff_df[f'cut_{buff}_Rall'].astype(int)

            if buff == 'C':
                buff_df['CD'] = 0
                buff_df['CG'] = 0
            elif buff == 'CD':
                buff_df['CD'] = 1
                buff_df['CG'] = 0
            elif buff == 'CG':
                buff_df['CD'] = 0
                buff_df['CG'] = 1  
            buff_df = buff_df[['gene','cut', 'region', 'AA', 'CD', 'CG']]              
            print(buff_df)
            dfs += [buff_df]

        pool_df = pd.concat(dfs)
        print(pool_df)

        common_genes = set(all_genes[0])
        for genes in all_genes[1:]:
            common_genes = common_genes.intersection(genes)
        #print(common_genes, len(common_genes))

        #pool_df = pool_df[pool_df['gene'].isin(common_genes)]
        #print(pool_df)
        
        #print(pool_df[['cut','region', 'CD', 'CG']].sum())
        return pool_df

    ###################################################################


    def run(self):
        """
        Orchestrates the workflow by loading data, performing regression, and saving results.
        """
        start_time = time.time()


        # Load data
        self.load_data()

        # Pool data
        essential_df = self.pool_data('essential')
        nonessential_df = self.pool_data('nonessential')
        
        # Get distribution of cuts per protein
        #self.get_cut_dist()

        # Perform regression
        #reg_formula = 'cut ~ region*(cut_CD_Rall + cut_CG_Rall)'
        reg_formula = 'cut ~ AA + region*(CD + CG)'
        #reg_formula = 'cut ~ AA + region*(CG)'
        #reg_formula = 'cut ~ region*(CD + CG)'
        #reg_formula = 'cut ~ region + AA'
        reg = self.regression(essential_df, reg_formula)
        reg['essential'] = True
        print(reg)
        dfs = [reg]
        reg = self.regression(nonessential_df, reg_formula)
        print(reg)
        reg['essential'] = False
        reg = pd.concat(dfs + [reg])

        reg['coef'] = reg['coef'].astype(float)
        reg['OR'] = np.exp(reg['coef'].astype(float))
        reg['std err'] = reg['std err'].astype(float)
        reg['z'] = reg['z'].astype(float)
        reg['P>|z|'] = reg['P>|z|'].astype(float)
        reg['[0.025'] = np.exp(reg['[0.025'].astype(float))
        reg['0.975]'] = np.exp(reg['0.975]'].astype(float))
        print(f'Regression Results:\n{reg.to_string()}')

        reg_outfile = os.path.join(self.outpath, 'RegressionPooled_v2.1.csv')
        reg.to_csv(reg_outfile, index=False)
        print(f"SAVED: {reg_outfile}")

        print('NORMAL TERMINATION')

def main():
    """
    Main function to parse arguments and run the DataAnalysis class.
    """
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-f", "--resFeat_files", type=str, required=True, help="Path to residue feature files")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("-t", "--tag", type=str, required=True, help="Tag for output filenames")

    args = parser.parse_args()

    analysis = DataAnalysis(
        resFeat_files=args.resFeat_files,
        outpath=args.outpath,
        tag=args.tag)
    analysis.run()

if __name__ == "__main__":
    main()

