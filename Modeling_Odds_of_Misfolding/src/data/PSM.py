import logging
import glob
import statsmodels.api as sm
import argparse
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from scipy.stats import ks_2samp
from scipy.spatial.distance import euclidean
import matplotlib.pyplot as plt
import os
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
import multiprocessing
import time

class PropensityScoreMatching:
    """
    A class to perform Propensity Score Matching (PSM) and related statistical analyses.
    This class handles loading data, performing regression analysis, PSM, and generating gene-level statistics.
    """

    def __init__(self, data_path, log_path, outpath, match_var, nmatch):
        """
        Initializes the PropensityScoreMatching class with necessary paths and parameters.

        Parameters:
        - data_path (str): Path to the residue feature files.
        - log_path (str): Path to the log file.
        - outpath (str): Path to the output directory.
        - match_var (str): Variable to match.
        - nmatch (int): Number of matches to attempt.
        """
        self.data_path = data_path
        self.log_path = log_path
        self.outpath = outpath
        self.match_var = match_var
        self.nmatch = nmatch
        self.data = None
        self.slim_df_for_matching = None
        self.logger = self.setup_logging()
        
        self.ent_df = None
        self.nonent_df = None

        if not os.path.exists(self.outpath):
            os.makedirs(self.outpath)

    def setup_logging(self):
        """
        Sets up the logging configuration.

        Returns:
        - logger (logging.Logger): Configured logger.
        """
        logging.basicConfig(filename=self.log_path, level=logging.INFO, format='%(asctime)s %(message)s') 
        logger = logging.getLogger(__name__)
        return logger

    def load_data(self):
        """
        Loads the residue feature files and filters the data for analysis.
        """
        for fi, f in enumerate(glob.glob(self.data_path)):
            if fi == 0:
                df = pd.read_csv(f, delimiter='|')
            else: 
                df = pd.concat((df, pd.read_csv(f, delimiter='|')))
        

        df = df[df['mapped_resid'].notna()].reset_index(drop=True)
        self.data = df
        print(f'self.data\n{self.data}')
        self.logger.info(f"Data loaded and filtered. Number of unique genes: {len(df['gene'].unique())}")

    def R_pmatch(self, data, caliper, treatment_variable='region', covariates=['num_nearest_neighbors', 'res_sasa']):
        """
        Performs propensity score matching using R's MatchIt package.

        Parameters:
        - data (pd.DataFrame): DataFrame containing the data for matching.
        - caliper (float): Caliper value for PSM.
        - treatment_variable (str): The treatment variable in the data.
        - covariates (list): List of covariate column names to be used in the matching.

        Returns:
        - matched_data_pd (pd.DataFrame): DataFrame containing the matched data.
        - nones (list): Indices of non-matched entries.
        """
        pandas2ri.activate()
        utils = importr('utils')

        if not robjects.packages.isinstalled("MatchIt"):
            utils.install_packages("MatchIt")
        
        matchit = importr('MatchIt')

        # Convert pandas DataFrame to R data.frame
        with localconverter(robjects.default_converter + pandas2ri.converter):
            r_data = robjects.conversion.py2rpy(data)

        # Specify the formula for propensity score matching
        formula_str = treatment_variable + ' ~ ' + ' + '.join(covariates)
        formula = robjects.Formula(formula_str)
        #print(f'formula: {formula}')

        # Perform propensity score matching
        #matched_data = matchit.matchit(formula=formula, data=r_data, method='nearest', distance = 'glm', ratio=nmatch, replace=False, caliper=caliper)
        match_result = matchit.matchit(formula=formula, data=r_data, method='nearest', distance='glm', ratio=self.nmatch, replace=False, caliper=caliper)

        # Convert matchit str matrix to numpy array
        with localconverter(robjects.default_converter + pandas2ri.converter):
            matched_data_pd = robjects.conversion.rpy2py(match_result.rx2('match.matrix'))
        
        Tidx = np.asarray(matched_data_pd.rownames).astype(int)
        match_matrix = np.asarray(matched_data_pd)
        nones = np.where(match_matrix == None)[0]
        match_matrix = match_matrix[match_matrix != None].astype(int)

        T = data.loc[Tidx]
        matchedC = data.loc[match_matrix.flatten()]
        matched_data_pd = pd.concat((T, matchedC))
        return matched_data_pd, nones

    def match_and_collect(self, caliper):
        """
        Matches data based on the specified caliper and collects the results.

        Parameters:
        - caliper (float): Caliper value for PSM.

        Returns:
        - tuple: Caliper value, matched data, and non-matched entries.
        """
        matched_data, nones = self.R_pmatch(self.slim_df_for_matching, caliper, treatment_variable='region', covariates=[self.match_var])
        return caliper, matched_data, nones

    def run(self):
        """
        Orchestrates the workflow by loading data, performing PSM, and generating gene-level statistics.
        """
        start_time = time.time()

        # Load data
        self.load_data()
        
        # Filter data for matching
        self.ent_df = self.data[self.data['region'] == 1]
        self.nonent_df = self.data[self.data['region'] == 0]
        self.slim_df_for_matching = pd.concat((self.ent_df, self.nonent_df))[['region', 'num_nearest_neighbors', 'res_sasa', 'median_sasa', 'mapped_resid', 'NC', 'crossing', 'nearest_neighbors', 'pdb_resid']]
        
        calipers = [0.01] + np.arange(0.05, 1.05, 0.05).tolist()
        with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
            results = pool.map(self.match_and_collect, calipers)
        
        matched_data_dict = {res[0]: {'df': res[1], 'nones': res[2]} for res in results}
        for caliper, data in matched_data_dict.items():
            if len(data['nones']) == 0:
                matched_data = data['df']
                opt_caliper = caliper
                break

        print(f'opt_caliper: {opt_caliper}')
        self.logger.info(f'opt_caliper: {opt_caliper}')
        df = self.data.loc[matched_data.index]
        print(f'df:\n{df}')

        outfile = os.path.join(self.outpath, f"matched_residue_features_{self.match_var}_n{self.nmatch}.csv")
        df.to_csv(outfile, sep='|', index=False)
        print(f'SAVED: {outfile}')

        print(f'NORMAL TERMINATION {time.time() - start_time}')

def main():
    """
    Main function to parse arguments and run the PropensityScoreMatching class.
    """
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-f", "--inpfiles", type=str, required=True, help="Path to residue feature files")
    parser.add_argument("-l", "--logfile", type=str, required=True, help="Path to log file")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("-m", "--match_var", type=str, required=True, help="Variable to match")
    parser.add_argument("-n", "--nmatch", type=int, required=True, help="Number of matches to attempt")
    args = parser.parse_args()

    psm = PropensityScoreMatching(data_path=args.inpfiles, log_path=args.logfile, outpath=args.outpath, match_var=args.match_var, nmatch=args.nmatch)
    psm.run()

if __name__ == "__main__":
    main()

