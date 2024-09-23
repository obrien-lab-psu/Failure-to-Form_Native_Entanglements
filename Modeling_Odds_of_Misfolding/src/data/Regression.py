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

    def __init__(self, resFeat_files, outpath, gene_list, tag, load_style, buffer, spa, cov, reg_formula, hold_var):
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
        self.gene_list = gene_list
        self.tag = tag
        self.load_style = load_style
        self.buffer = buffer
        self.spa = spa
        self.cov = cov
        self.reg_formula = reg_formula
        self.data = {}
        #self.logger = self.setup_logging()
        #self.gene_list_files = glob.glob(self.gene_list)
        self.hold_var = hold_var

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
                df[column] = label_encoder.fit_transform(df[column])
            else:
                print(f"Column '{column}' does not exist in the DataFrame.")
        
        return df

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

        # Get the cov_params
        cov_matrix = result.cov_params()
        print(f'cov_matrix:\n{cov_matrix}')
    
        # Get the coefficients
        print("Coefficients:")
        coefficients = {'A': 0}
        for k,v in result.params.items():
            if 'AA' in k:
                k = k.replace('AA[T.', '').replace(']', '')
            coefficients[k] = v
        for k,v in coefficients.items():
            print(k,v)
        self.coefficients = coefficients

        ## calculate the odds from the fit holding either region or ent_present constant
        if 'AA' in formula: 
            odds = []
            for value in df[self.hold_var].unique():
                odds += [self.calculate_overall_odds_with_confidence_intervals(df, coefficients, cov_matrix, self.hold_var, value=value, confidence_level=0.95)]
            odds = pd.DataFrame(odds)
        else:
            odds = None
        print(f'odds:\n{odds}')

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
        print(f"Data loaded and filtered. Number of unique genes: {len(self.data['gene'].unique())}")


    def calculate_overall_odds_with_confidence_intervals(self, dataframe, coefficients, cov_matrix, hold, value=1, confidence_level=0.95):
        # Calculate the frequency of each AA in the sample data
        freq_AA = dataframe['AA'].value_counts(normalize=True)

        # Initialize variables to store weighted log-odds and variance
        weighted_log_odds = 0
        weighted_var_log_odds = 0

        # Calculate z-score for the desired confidence level
        z = norm.ppf(1 - (1 - confidence_level) / 2)
        print(z, 1 - (1 - confidence_level) / 2)

        #var_intercept = cov_matrix.loc['Intercept', 'Intercept']
        for aa, freq in freq_AA.items():
            # Intercept + coefficient for AA + coefficient for hold
            coef_intercept = coefficients['Intercept']
            coef_aa = coefficients.get(aa, 0)
            coef_hold = coefficients[hold]

            log_odds = coef_intercept + coef_aa + coef_hold * value

            # Calculate standard error for log-odds
            var_intercept = cov_matrix.loc['Intercept', 'Intercept']
            var_aa = cov_matrix.loc[f'AA[T.{aa}]', f'AA[T.{aa}]'] if f'AA[T.{aa}]' in cov_matrix.index else 0
            var_hold = cov_matrix.loc[hold, hold]
            cov_ia = cov_matrix.loc['Intercept', f'AA[T.{aa}]'] if f'AA[T.{aa}]' in cov_matrix.index else 0
            cov_ir = cov_matrix.loc['Intercept', hold]
            cov_ar = cov_matrix.loc[f'AA[T.{aa}]', hold] if f'AA[T.{aa}]' in cov_matrix.index else 0

            # Total variance for the log-odds
            var_log_odds = var_intercept + var_aa + var_hold + 2 * (cov_ia + cov_ir + cov_ar)
            if abs(var_log_odds) > 10:
                var_log_odds = 10
            #print(f'{aa} {var_intercept} {var_aa} {var_hold} {cov_ia} {cov_ir} {cov_ar}')
            print(f'{aa} var_log_odds = {var_log_odds}')

            # Weighted log-odds and variance
            weighted_log_odds += freq * log_odds
            weighted_var_log_odds += freq ** 2 * var_log_odds

        # Standard error
        se_weighted_log_odds = np.sqrt(weighted_var_log_odds)

        # Confidence intervals for weighted log-odds
        weighted_log_odds_lower = weighted_log_odds - z * se_weighted_log_odds
        weighted_log_odds_upper = weighted_log_odds + z * se_weighted_log_odds

        # Convert to odds
        overall_odds = np.exp(weighted_log_odds)
        overall_odds_lower = np.exp(weighted_log_odds_lower)
        overall_odds_upper = np.exp(weighted_log_odds_upper)

        return {
            'var': hold,
            'hold_value': value,
            'overall_odds': overall_odds,
            'CI_lower': overall_odds_lower,
            'CI_upper': overall_odds_upper}

    def get_cut_dist(self, ):
        print(f'Getting distribution of cut sites for this set of proteins')
        #outfile = os.path.join(self.outpath, f"CutDist_{self.hold_var}_{self.tag}_{self.buffer}_{self.timepoint}_spa{self.spa}_LiPMScov{self.cov}.png")
        outfile_csv = os.path.join(self.outpath, f"CutDist_{self.hold_var}_{self.tag}_{self.buffer}_{self.timepoint}_spa{self.spa}_LiPMScov{self.cov}.csv")
        cut_dist_df = {'gene':[], 'cuts':[]}
        for gene, gene_df in self.data.groupby('gene'):
            #print(gene_df)
            cuts = np.sum(gene_df[self.cut_key])
            cut_dist_df['gene'] += [gene]
            cut_dist_df['cuts'] += [cuts]

        cut_dist_df = pd.DataFrame(cut_dist_df)
        cut_dist_df.to_csv(outfile_csv, sep='|', index=False)
        print(f'SAVED: {outfile_csv}')
        print(cut_dist_df)


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
        
        # Get data for gene list requested by user
        self.data = self.data[self.data['gene'].isin(self.reg_genes)]

        # Quality check to ensure all columns are present in the df
        reg_vars = [v for v in self.reg_formula.split(' ') if v not in ['~', '+']]
        print(reg_vars)
        keys = []
        for reg_var in reg_vars:
            if '*' in reg_var:
                reg_var = reg_var.split('*')
            else:
                reg_var = [reg_var]
            for v in reg_var:
                keys += [v]

                if 'cut_' in v:
                    cut_key = v

        keys = list(set(keys))
        bin_keys = keys.copy()
        self.cut_key = cut_key
        if 'AA' in keys:
            bin_keys.remove('AA')
        print(f'bin_keys: {bin_keys}')

        ## Find cut key
        print(f'cut_key: {cut_key}')
        if 'Rall' not in cut_key:
            timepoint = cut_key.split('_')[-1]
            self.timepoint = timepoint
            print(f'timepoint: {timepoint}')

            cuts = []
            for cut_str in self.data['cut_str'].values:
                #C,R1min,S2,-2.51/C,R2hr,S2,-1.79/C,R5min,S2,-2.03
                if isinstance(cut_str, str):
                    #print(cut_str)
                    if f'{self.buffer},{timepoint}' in cut_str:
                        cuts += [1]
                    else:
                        cuts += [0]
                else:
                    cuts += [0]
            self.data[cut_key] = cuts
            num_cuts = sum(cuts)
        else:
            num_cuts = sum(self.data[cut_key].values)
            self.timepoint = 'Rall'
        print(f'num_cuts: {num_cuts}')

        # Define output files and get gene list
        reg_outfile = os.path.join(self.outpath, f"regression_results_{self.hold_var}_{self.tag}_{self.buffer}_{self.timepoint}_spa{self.spa}_LiPMScov{self.cov}.csv")
        cov_matrix_outfile = os.path.join(self.outpath, f"regression_cov_matrix_{self.hold_var}_{self.tag}_{self.buffer}_{self.timepoint}_spa{self.spa}_LiPMScov{self.cov}.csv")
        odds_outfile = os.path.join(self.outpath, f"regression_odds_holding-{self.hold_var}_{self.tag}_{self.buffer}_{self.timepoint}_spa{self.spa}_LiPMScov{self.cov}.csv")

        # Encode boolean columns
        encoded_df = self.encode_boolean_columns(self.data, bin_keys)
        encoded_df = encoded_df[reg_vars]
        print(f'encoded_df:\n{encoded_df}')
        self.data = self.data[bin_keys + ['gene']]
      
        # Get distribution of cuts per protein
        self.get_cut_dist()

        # Perform Fisher Exact test
        table = pd.crosstab(encoded_df[cut_key], encoded_df['region'])
        print(f'Table:\n{table}')
        if table.values.shape == (2,2):
            OR, pvalue = fisher_exact(table) 
            print(f'Fisher exact: OR {OR} with pvalue {pvalue}')

        # Perform regression
        reg, odds, cov_matrix = self.regression(encoded_df, self.reg_formula)
        reg['coef'] = reg['coef'].astype(float)
        reg['OR'] = np.exp(reg['coef'].astype(float))
        reg['std err'] = reg['std err'].astype(float)
        reg['z'] = reg['z'].astype(float)
        reg['P>|z|'] = reg['P>|z|'].astype(float)
        reg['[0.025'] = np.exp(reg['[0.025'].astype(float))
        reg['0.975]'] = np.exp(reg['0.975]'].astype(float))
        reg['tag'] = self.tag
        reg['buff'] = self.buffer
        reg['timepoint'] = self.timepoint
        reg['spa'] = self.spa
        reg['cov'] = self.cov
        reg['n'] = len(self.data['gene'].unique())
        reg['num_cuts'] = num_cuts
        reg['num_res'] = len(encoded_df)
        reg['num_ent-res'] = len(encoded_df[encoded_df['region'] == 1])
        reg['num_nonent-res'] = len(encoded_df[encoded_df['region'] == 0])
        print(f'Regression Results:\n{reg.to_string()}')

        reg.to_csv(reg_outfile, index=False, sep='|')
        print(f"SAVED: {reg_outfile}")

        cov_matrix.to_csv(cov_matrix_outfile, index=False, sep='|')
        print(f'SAVED: {cov_matrix_outfile}')

        if 'AA' in self.reg_formula:
            odds['tag'] = self.tag
            odds['buff'] = self.buffer
            odds['timepoint'] = self.timepoint
            odds['spa'] = self.spa
            odds['cov'] = self.cov
            odds['n'] = len(self.data['gene'].unique())
            pval = reg[reg['var'] == self.hold_var]['P>|z|'].values
            odds['P>|z|'] = pval[0]
            print(f'odds:\n{odds}')
            odds.to_csv(odds_outfile, index=False, sep='|')
            print(f'SAVED: {odds_outfile}')

        print('NORMAL TERMINATION')

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
    parser.add_argument("-r", "--reg_formula", type=str, required=True, help="Regression formula")
    parser.add_argument("-v", "--hold_var", type=str, required=True, help="Variable to hold constant while calculating odds")

    args = parser.parse_args()

    analysis = DataAnalysis(
        resFeat_files=args.resFeat_files,
        outpath=args.outpath,
        gene_list=args.gene_list,
        tag=args.tag,
        load_style=args.load_style,
        buffer=args.buffer,
        spa=args.spa,
        cov=args.cov,
        reg_formula=args.reg_formula,
        hold_var=args.hold_var)
    analysis.run()

if __name__ == "__main__":
    main()

