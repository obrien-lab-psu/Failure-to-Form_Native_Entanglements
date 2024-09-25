import pandas as pd
import numpy as np
import glob
import os
import argparse
from scipy import stats
import statsmodels.stats.multitest as smm

class LiPMSProcessor:
    """
    Class to process and correct LiPMS files for peptide abundance changes.
    """

    def __init__(self, lipms_files, threshold, outpath, alpha):
        """
        Initialize the LiPMSProcessor with file paths, threshold, output path, and alpha for FWER.

        :param lipms_files: Path pattern to LiPMS files to be corrected.
        :param threshold: Integer threshold for abundance change.
        :param outpath: Path where corrected files will be saved.
        :param alpha: Family wise error rate (alpha value).
        """
        self.lipms_files = lipms_files
        self.threshold = threshold
        self.outpath = outpath
        self.alpha = alpha
        self.create_outpath()

    def create_outpath(self):
        """
        Create the output directory if it does not exist.
        """
        full_outpath = os.path.join(self.outpath, f'logRN_{self.threshold}', f'FWER_{self.alpha}')
        self.full_outpath = full_outpath
        if not os.path.exists(self.full_outpath):
            os.makedirs(self.full_outpath)
            print(f'MADE: {self.full_outpath}')

    def process_files(self):
        """
        Process each file found by the glob pattern.
        """
        lip_files = glob.glob(self.lipms_files)
        for i, file_path in enumerate(lip_files):
            self.process_file(i, file_path)

    def process_file(self, index, file_path):
        """
        Process a single file to correct peptide data.

        :param index: Index of the file in the list.
        :param file_path: Path to the file to be processed.
        """

        ## path to outfile containing all peptids
        all_outf = os.path.basename(file_path).replace('.csv', '_corrected.csv')
        all_outf = os.path.join(self.full_outpath, all_outf)

        ## path to outfile containing only PK peptides
        sig_outf = os.path.basename(file_path).replace('.csv', '_sig_corrected.csv')
        sig_outf = os.path.join(self.full_outpath, sig_outf)
        print('\n', index, file_path, all_outf, sig_outf)

        df = pd.read_csv(file_path)

        df = df[['Accession', 'proteinaseKsite', 'Peptide Sequence', 'PeptideRatio1', 'PeptidePValue1']]
        df['pvalue_raw'] = 10 ** (-1 * df['PeptidePValue1'])
        df = df[~df["proteinaseKsite"].str.contains('-')]

        df = self.apply_corrections(df)

        # save all corrected PK peptides
        df.to_csv(all_outf, index=False)
        print(f'SAVED: {all_outf} {len(df)}')

        # save only sig corrected peptides that meet abundance threshold
        # and meet FWER threshold
        df = df[abs(df['PeptideRatio1']) >= self.threshold]
        df = df[df['FDR_pvalue'] <= self.alpha]
        df.to_csv(sig_outf, index=False)
        print(f'SAVED: {sig_outf} {len(df)}')

        # get the number of significant peptides
        num_sig = np.sum(df['FDR_pvalue'] <= self.alpha)
        print(f'num_sig: {num_sig} {len(df)}')


    def apply_corrections(self, df):
        """
        Apply corrections to DataFrame containing peptide data.

        :param df: DataFrame to apply corrections to.
        :return: Corrected DataFrame.
        """

        # Apply the Benjamini-Hochberg correction
        pvals = df['pvalue_raw'].values
        pvals_corrected = stats.false_discovery_control(pvals)
        df['FDR_pvalue'] = pvals_corrected

        return df.reset_index(drop=True)

def main():
    """
    Main function to handle command line arguments and run the processing.
    """
    parser = argparse.ArgumentParser(description="Process user specified arguments for LiPMS files correction.")
    parser.add_argument("-f", "--lipms_files", type=str, required=True, help="path to LiPMS files to correct")
    parser.add_argument("-t", "--threshold", type=int, required=True, help='abundance change threshold')
    parser.add_argument("-o", "--outpath", type=str, required=True, help="outpath")
    parser.add_argument("-a", "--alpha", type=float, required=True, help="family wise error rate: alpha")
    args = parser.parse_args()

    processor = LiPMSProcessor(args.lipms_files, args.threshold, args.outpath, args.alpha)
    processor.process_files()

    print('NORMAL TERMINATION')

if __name__ == "__main__":
    main()

