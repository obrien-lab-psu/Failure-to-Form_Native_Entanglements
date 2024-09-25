import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, combine_pvalues
import pickle
import numpy as np
import sys,os
import glob
import re
import argparse
import warnings


class DataAnalysis:
    """
    A class to handle the data analysis process including encoding, regression, and statistical tests.
    """

    #######################################################################################
    def __init__(self, inpfiles, outpath, tag):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
        - resFeat_files (str): Path to residue feature files.
        - outpath (str): Path to the output directory.
        """
        files = {}
        for f in glob.glob(inpfiles):
            buff = f.split('/')[-1].split('_')[0]
            timepoint = f.split('/')[-1].split('_')[1]
            print(f, buff, timepoint)
            files[f] = {'buff':buff, 'timepoint':timepoint}
    
        self.inpfiles = files
        self.outpath = outpath
        self.tag = tag

        if not os.path.exists(f'{self.outpath}'):
            os.makedirs(f'{self.outpath}')
            print(f'Made output directories {self.outpath}')

        self.data_path = os.path.join(f'{self.outpath}', 'Data/')
        if not os.path.exists(f'{self.data_path}'):
            os.makedirs(f'{self.data_path}')
            print(f'Made output directories {self.data_path}')

        self.plots_path = os.path.join(f'{self.outpath}', 'Plots/')
        if not os.path.exists(f'{self.plots_path}'):
            os.makedirs(f'{self.plots_path}')
            print(f'Made output directories {self.plots_path}')

    #######################################################################################
    def setup_logging(self, logfile):
        """
        Sets up the logging configuration.

        Returns:
        - logger (logging.Logger): Configured logger.
        """
        logging.basicConfig(filename=logfile, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        logger = logging.getLogger(__name__)
        return logger


    #######################################################################################
    def threshold_data(self,):

        if self.tag == 'spa':
            col = 'N_spa'
        elif self.tag == 'cov':
            col = 'coverage'

        master_dict = {}
        ### Get the thresholding for the sum of the average peptide abundance 
        ### per gene in each condition separatly
        for f, f_params in self.inpfiles.items():
            print(f, f_params)

            data = pd.read_csv(f)
            results = self.analyze_cdf(data, f_params)
            master_dict[(f_params['buff'], f_params['timepoint'])] = results

        
        ### Get the thresholding for the sum of the avverage peptide abundance
        ### per gene across the three time points (R1min, R5min, R2hr) for a given buffer C CD CG
        for buff in ['C', 'CD', 'CG']:

            buff_files = [f for f in self.inpfiles if f'{buff}_' in f]
            buff_data = {}

            # for the timepoint files for this buffer make a dictionary where the keys are the genes and the values are the list of either SPA or LiPMS cov values
            for bf in buff_files:
                print(bf)
                data = pd.read_csv(bf)
                print(data)
                for rowi, row in data.iterrows():
                    gene = row['Accession']
                    if self.tag == 'spa':
                        value = row['N_spa']
                    elif self.tag == 'cov':
                        value = row['coverage']
                    else:
                        raise ValueError(f"tag value is invalud {self.tag}. Must be either spa or cov")

                    if gene in buff_data:
                        buff_data[gene] += [value]
                    else:
                        buff_data[gene] = [value]
            
            # Make a dataframe that is all the genes observed across the three time points where the col is the mean of their SPA or LiPMScov value
            buff_df = {'Accession':[], col:[]}
            for gene, vals in buff_data.items():
                buff_df['Accession'] += [gene]
                buff_df[col] += [np.mean(vals)]
            buff_df = pd.DataFrame(buff_df)
            print(buff_df)

            results = self.analyze_cdf(buff_df, {'buff':buff, 'timepoint':'Rall'})
            master_dict[(buff, 'Rall')] = results


        ### SAVE the master dictionary to a pickle file that  contains all the info
        master_outfile = f'{self.data_path}{self.tag}_threshold_master.pkl'
        with open(master_outfile, 'wb') as fh:
            pickle.dump(master_dict, fh)
        print(f'SAVED: {master_outfile}')


    #######################################################################################
    def analyze_cdf(self, dataframe, params):
        buff = params['buff']
        timepoint = params['timepoint']
        outfile = f'{self.plots_path}{buff}_{timepoint}_CDF_vs_{self.tag}.png'
        print(buff, timepoint, outfile)

        if self.tag == 'spa':
            col = 'N_spa'
        elif self.tag == 'cov':
            col = 'coverage'
        

        # Calculate the CDF
        sorted_values = np.sort(dataframe[col])

        cdf = np.arange(1, len(sorted_values)+1) / len(sorted_values)

        # Plot the CDF
        plt.figure(figsize=(6, 6))
        plt.plot(sorted_values, cdf, marker='.', linestyle='none')
        plt.xlabel(f'<{self.tag}>')
        plt.ylabel('CDF')
        plt.title(f'CDF of {self.tag} | {buff} {timepoint}')
        plt.grid(True)
        if self.tag == 'spa':
            plt.xscale('log')
            plt.yticks(np.arange(0, 1.1, 0.1))
        #plt.show()
        plt.savefig(outfile)
        print(f'SAVED: {outfile}')

        # Calculate the specified percentiles
        percentiles = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]
        percentile_values = np.percentile(dataframe[col], percentiles)

        # Create the result dictionary
        result = {}
        for i, percentile in enumerate(percentiles):
            threshold_value = percentile_values[i]
            threshold_df = dataframe.loc[dataframe[col] >= threshold_value]
            result[(percentile, threshold_value)] =  threshold_df
        return result
        

########################################################################################################
def main():
    """
    Main function to parse arguments and run the DataAnalysis class.
    """

    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-i", "--inpfiles", type=str, required=True, help="Path to raw LiPMS file to process")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("-t", "--tag", type=str, required=True, help="Thresholding tag: spa or cov")
    args = parser.parse_args()

    inpfiles = args.inpfiles
    outpath = args.outpath
    tag = args.tag

    Analyzer = DataAnalysis(inpfiles, outpath, tag)
    print(Analyzer)

    # threshold SPA data
    Analyzer.threshold_data()

    print('NORMAL TERMINATION')


if __name__ == "__main__":
    main()

