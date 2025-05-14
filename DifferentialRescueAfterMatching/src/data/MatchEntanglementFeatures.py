import time, sys
import multiprocessing as mp
from scipy.stats import bootstrap
import logging
import argparse
import glob
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import NearestNeighbors
import os
#pd.set_option('display.max_rows', 4000)

class DataAnalysis:
    """
    A class to handle the data analysis process including encoding, regression, and statistical tests.
    """

    def __init__(self, uentFeat_files, outpath, gene_list_path):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
        - uentFeat_files (str): Path to residue feature files.
        - outpath (str): Path to the output directory.
        - gene_list_path (str): Path to gene lists to use.
        """
        self.uentFeat_files = uentFeat_files
        self.outpath = outpath
        self.gene_list_path = gene_list_path

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



    def load_data(self):
        """
        Loads the residue feature files and filters the data for analysis.
        """
        gene_lists = glob.glob(self.gene_list_path + '/*.txt')
        print(f'Number of gene_lists: {len(gene_lists)}')

        ## Make a dictionary that contains the residue feature dfs for each gene list in self.gene_list_path
        uentFeat_files = glob.glob(self.uentFeat_files + '/*.csv')
        print(f'Number of uentFeat feature files: {len(uentFeat_files)}')
        if len(uentFeat_files) == 0:
            print(f"No uentFeat feature files found in {self.uentFeat_files}")
        
        self.data = {}
        keys = ['Gc', 'ent_coverage', 'min_N_prot_depth_left', 'min_C_prot_depth_right', 'min_C_thread_slippage_right', 'RCO']
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
            self.gene_list_uentFeat = pd.DataFrame()
            for i, gene in enumerate(genes):
                gene_uentFeat = [f for f in uentFeat_files if gene in f]
                if len(gene_uentFeat) == 0:
                    print(f"No residue feature file found for gene {gene}")
                    quit()
                elif len(gene_uentFeat) > 1:
                    print(f"More than 1 residue feature file found for gene {gene}")
                    quit()
                gene_uentFeat_file = gene_uentFeat[0]
                #print(f'gene_uentFeat_file: {gene_uentFeat_file} {i}')
                if len(self.gene_list_uentFeat) == 0:
                    self.gene_list_uentFeat = pd.read_csv(gene_uentFeat_file, sep='|')[keys]
                    self.gene_list_uentFeat.fillna(0, inplace=True)
                    self.gene_list_uentFeat = self.gene_list_uentFeat.max().to_frame().T
                    self.gene_list_uentFeat['gene'] = gene
                    self.gene_list_uentFeat['essential'] = essential

                else:
                    df = pd.read_csv(gene_uentFeat_file, sep='|')[keys]
                    df.fillna(0, inplace=True)
                    df = df.max().to_frame().T
                    df['gene'] = gene
                    df['essential'] = essential
                    self.gene_list_uentFeat = pd.concat((self.gene_list_uentFeat, df))
            
            ## encode the columns for cut_{buff}_Rall, region, AA
            self.data[buffer][essential] = self.gene_list_uentFeat
            print(f'Number of genes in self.gene_list_uentFeat: {len(self.gene_list_uentFeat["gene"].unique())}')

    def DataMatch(self, data, treated_str='nonessential', control_str='essential', buffer='C', k=1):
        """
        Match the data from the essential to nonessential data

        Parameters:
        - data (dict): A dictionary containing the data for essential and nonessential proteins.
        - data2matchtoo (str): The key for the data to match to.
        - data2match (str): The key for the data to match.
        - buffer (str): The buffer to match.
        """
        print(f'\nMatching {control_str} to {treated_str} for buffer {buffer}')
 
        treated_df = data[buffer][treated_str]
        treated_df['category'] = 1
        control_df = data[buffer][control_str]
        control_df['category'] = 0
        #print(f'treated_df:\n{treated_df}')
        #print(f'control_df:\n{control_df}')

        df = pd.concat((treated_df, control_df))
        treatment_col = 'category'
        feature_cols = ['Gc', 'ent_coverage', 'min_N_prot_depth_left', 'min_C_prot_depth_right', 'min_C_thread_slippage_right', 'RCO']

        # Prematching feature check
        for key in feature_cols:
            treated_mean = treated_df[key].mean()
            control_mean = control_df[key].mean()
            print(f"Quality check for pre-matching: {key} {treated_mean} {control_mean}")


        # Fit logistic regression model to estimate propensity scores
        model = LogisticRegression()
        model.fit(df[feature_cols], df[treatment_col])
        
        # Compute propensity scores
        df['propensity_score'] = model.predict_proba(df[feature_cols])[:, 1]

        # Separate treated and control groups
        treated = df[df[treatment_col] == 1].copy()
        control = df[df[treatment_col] == 0].copy()

        # Use Nearest Neighbors to match based on propensity scores
        nn = NearestNeighbors(n_neighbors=k, metric='euclidean')
        nn.fit(control[['propensity_score']])  # Fit on control group propensity scores
        distances, indices = nn.kneighbors(treated[['propensity_score']])  # Find nearest matches

        # Get matched control group
        matched_control = control.iloc[indices.flatten()].reset_index(drop=True)

        # Combine matched treated and control groups
        matched_df = pd.concat([treated.reset_index(drop=True), matched_control], axis=0)
        print(f'matched_df:\n{matched_df}')

        ## Quality check for matching
        print(f"Quality check for matching: {matched_df['category'].value_counts()}")
        for key in feature_cols:
            treated_mean = treated[key].mean()
            matched_control_mean = matched_control[key].mean()
            print(f"Quality check for matching: {key} {treated_mean} {matched_control_mean}")
        

        ## save the list of genes from the matched dataset
        matched_outpath = os.path.join(self.outpath, f'{control_str}_matched2_{treated_str}')
        if not os.path.exists(f'{matched_outpath}'):
            os.makedirs(f'{matched_outpath}')
            print(f'Made output directories {matched_outpath}')

        essential_filename = f'{matched_outpath}/EXP_0.6g_{buffer}_Rall_spa50_LiPMScov50_essential_ent_genes.txt'
        nonessential_filename = f'{matched_outpath}/EXP_0.6g_{buffer}_Rall_spa50_LiPMScov50_nonessential_ent_genes.txt'

        essential_genes = matched_df[matched_df['essential'] == 'essential']['gene'].values
        nonessential_genes = matched_df[matched_df['essential'] == 'nonessential']['gene'].values

        np.savetxt(essential_filename, essential_genes, fmt='%s')
        np.savetxt(nonessential_filename, nonessential_genes, fmt='%s')
        print(f'SAVED: {essential_filename}')
        print(f'SAVED: {nonessential_filename}')


    def run(self):
        """
        Orchestrates the workflow by loading data, performing regression, and saving results.
        """
        start_time = time.time()

        # Load data
        self.load_data()

        # Match the data Essential(control) to Nonessential(treated)
        self.DataMatch(self.data, treated_str='nonessential', control_str='essential', buffer='C', k=2)
        self.DataMatch(self.data, treated_str='nonessential', control_str='essential', buffer='CD', k=2)
        self.DataMatch(self.data, treated_str='nonessential', control_str='essential', buffer='CG', k=2)

        self.DataMatch(self.data, treated_str='essential', control_str='nonessential', buffer='C', k=2)
        self.DataMatch(self.data, treated_str='essential', control_str='nonessential', buffer='CD', k=2)
        self.DataMatch(self.data, treated_str='essential', control_str='nonessential', buffer='CG', k=2)

        print(f'NORMAL TERMINATION {time.time() - start_time} seconds')

def main():
    """
    Main function to parse arguments and run the DataAnalysis class.
    """
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-f", "--uentFeat_files", type=str, required=True, help="Path to unique entanglent files")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("-g", "--gene_list_path", type=str, required=True, help="Path to gene lists to use")

    args = parser.parse_args()

    analysis = DataAnalysis(
        uentFeat_files=args.uentFeat_files,
        outpath=args.outpath,
        gene_list_path=args.gene_list_path)
    analysis.run()

if __name__ == "__main__":
    main()

