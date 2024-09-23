import statsmodels.api as sm
import glob
import argparse
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from scipy.spatial.distance import euclidean
import matplotlib.pyplot as plt
import os
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import scipy.stats as st
import logging

class SASADistributionPlotter:
    """
    A class to handle the process of plotting SASA distributions from pre-matched and post-matched residue feature files.
    This class handles loading data, performing any necessary calculations, and plotting the distributions.
    """

    def __init__(self, pre_inpfiles, post_inpfiles, outpath, tag):
        """
        Initializes the SASADistributionPlotter class with paths for pre-matched and post-matched residue feature files,
        output directory, and a tag for figure titles.

        Parameters:
        - pre_inpfiles (str): Path to pre-matched residue feature files.
        - post_inpfiles (str): Path to propensity score matched residue feature files.
        - outpath (str): Path to the output directory.
        - tag (str): Tag for figure titles.
        """
        self.pre_inpfiles = pre_inpfiles
        self.post_inpfiles = post_inpfiles
        self.outpath = outpath
        self.tag = tag
        self.pre_data = None
        self.post_data = None
        self.logger = self.setup_logging()

        if not os.path.exists(self.outpath):
            os.makedirs(self.outpath)
            print(f'MADE: {self.outpath}')

    def setup_logging(self):
        """
        Sets up the logging configuration.

        Returns:
        - logger (logging.Logger): Configured logger.
        """
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)
        return logger

    def load_data(self):
        """
        Loads the pre-matched and post-matched residue feature files and logs the number of files loaded.
        """
        self.pre_data = glob.glob(os.path.join(self.pre_inpfiles, '*_resfeatures.csv'))
        self.post_data = glob.glob(self.post_inpfiles)

        ##############################################################################################################################################################
        ### Load pre-matched residue features
        pre_df = {}
        for fi, f in enumerate(self.pre_data):
            temp = pd.read_csv(f, delimiter='|')

            if len(pre_df) == 0:
                pre_df = temp
            else:
                pre_df = pd.concat((pre_df, temp))

        # only get the data that was mapped onto the genome
        pre_df = pre_df[pre_df['mapped_resid'].notna()]
        pre_df = pre_df.reset_index()
        self.pre_df = pre_df
        print(f'pre_df:\n{pre_df}\n')

        ##############################################################################################################################################################
        ### Load pscore matched residue features
        post_df = {}
        for fi, f in enumerate(self.post_data):
            gene = f.split('/')[-1].split('_')[0]
            if len(post_df) == 0:
                post_df = pd.read_csv(f, delimiter='|')
            else:
                post_df = pd.concat((post_df, pd.read_csv(f, delimiter='|')))

        # only get the data that was mapped onto the genome
        post_df = post_df[post_df['mapped_resid'].notna()]
        post_df = post_df.reset_index()
        self.post_df = post_df
        print(f'post_df:\n{post_df}\n')
        ##############################################################################################################################################################

        self.logger.info(f"Loaded {len(self.pre_data)} pre-matched files from {self.pre_inpfiles}")
        self.logger.info(f"Loaded {len(self.post_data)} post-matched files from {self.post_inpfiles}")

    def plot_distributions(self):
        """
        Plots the SASA distributions for the pre-matched and post-matched residue feature files.
        """
        self.logger.info("Plotting SASA distributions")

        ### Get data to plot

        #(1) prematched entR and nonentR
        pre_entR = self.pre_df[self.pre_df['region'] == 1]
        pre_nonentR = self.pre_df[self.pre_df['region'] == 0]

        #(2) prematched entR and nonentR res_sasa distribution
        pre_entR_res_sasa = pre_entR['res_sasa'].values
        pre_nonentR_res_sasa = pre_nonentR['res_sasa'].values

        #(3) prematched entR and nonentR median_sasa 
        pre_entR_median_sasa = pre_entR['median_sasa'].values
        pre_nonentR_median_sasa = pre_nonentR['median_sasa'].values


        #(7) postmatched entR and nonentR
        post_entR = self.post_df[self.post_df['region'] == 1]
        post_nonentR = self.post_df[self.post_df['region'] == 0]

        #(8) postmatched entR and nonentR res_sasa distribution
        post_entR_res_sasa = post_entR['res_sasa'].values
        post_nonentR_res_sasa = post_nonentR['res_sasa'].values

        #(9) postmatched entR and nonentR median_sasa 
        post_entR_median_sasa = post_entR['median_sasa'].values
        post_nonentR_median_sasa = post_nonentR['median_sasa'].values

        # Example plot (actual plotting logic needs to be implemented based on the data format)
        bins = 20
        x_range = (-0.1,2)
        print(f'\n{"#"*100}\nPlotting')

        fig, axs = plt.subplots(1, 2, figsize=(6, 3))

        ### second row of figure
        axs[0].hist(pre_entR_res_sasa, bins=bins, range=x_range, alpha=0.5, label='pre-entR', color='blue', density=True, edgecolor='blue')
        axs[0].hist(pre_nonentR_res_sasa, bins=bins, range=x_range, alpha=0.5, label='pre-nonentR', color='red', density=True, edgecolor='red')

        ## (0,1) the postmatched res_sasa for entR and non-entR regions
        axs[1].hist(post_entR_res_sasa, bins=bins, range=x_range, alpha=0.5, label='post-entR', color='blue', density=True, edgecolor='blue')
        axs[1].hist(post_nonentR_res_sasa, bins=bins, range=x_range, alpha=0.5, label='post-nonentR', color='red', density=True, edgecolor='red')

        axs[0].legend()
        axs[1].legend()

        axs[0].set_xlabel(f'res_sasa, nm^2')
        axs[0].set_ylabel(f'PDF')
        axs[1].set_xlabel(f'res_sasa, nm^2')
        axs[1].set_ylabel(f'PDF')

        # Adjust layout to prevent overlap
        plt.suptitle(f'{self.tag} pscore matching')
        plt.tight_layout()
        # reg_analysis_1_total_quasi_binomial_regression_results_
        outfile = f'{self.outpath}{self.tag}_matched_residue_features.png'
        #plt.show()
        plt.savefig(outfile)
        print(f'SAVED: {outfile}')
        plt.close()

        self.logger.info(f"Plot saved to {self.outpath}")

    def run(self):
        """
        Orchestrates the workflow by loading data and plotting the distributions.
        """
        self.load_data()
        self.plot_distributions()

def main():
    """
    Main function to parse arguments and run the SASADistributionPlotter class.
    """
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-p", "--pre_inpfiles", type=str, required=True, help="Path to pre-matched residue feature files")
    parser.add_argument("-f", "--post_inpfiles", type=str, required=True, help="Path to propensity score matched residue feature files")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("-t", "--tag", type=str, required=True, help="Tag for figure title")
    args = parser.parse_args()

    plotter = SASADistributionPlotter(pre_inpfiles=args.pre_inpfiles, post_inpfiles=args.post_inpfiles, outpath=args.outpath, tag=args.tag)
    plotter.run()

if __name__ == "__main__":
    main()
print(f'NORMAL TERMINATION')
