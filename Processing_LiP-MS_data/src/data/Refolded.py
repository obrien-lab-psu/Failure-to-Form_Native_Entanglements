import pandas as pd
import numpy as np
import glob
import os
import argparse
from scipy import stats
import statsmodels.stats.multitest as smm
pd.set_option('display.max_rows', 5000)
class LiPMSProcessor:
    """
    Class to process and correct LiPMS files for peptide abundance changes.
    """

    def __init__(self, lipms_files, threshold, outpath, alpha, gene_lists):
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
        self.ent_files = glob.glob(os.path.join(gene_lists, '*_Rall_spa50_LiPMScov50_ent_genes.txt'))
        print(f'Number of ent gene files loaded: {len(self.ent_files)}')
        self.all_files = glob.glob(os.path.join(gene_lists, '*_Rall_spa50_LiPMScov50_all_genes.txt'))
        print(f'Number of all gene files loaded: {len(self.ent_files)}')

    def create_outpath(self):
        """
        Create the output directory if it does not exist.
        """
        if not os.path.exists(self.outpath):
            os.makedirs(self.outpath)
            print(f'MADE: {self.outpath}')

    def process_files(self):
        """
        Process each file found by the glob pattern.
        """

        lip_files = glob.glob(self.lipms_files)
        all_refolded = []
        ent_refolded = []
        for i, file_path in enumerate(lip_files):
            allR, entR = self.process_file(i, file_path)
            all_refolded += [allR]
            ent_refolded += [entR]
        all_refolded = [l for l in all_refolded if l is not None]
        ent_refolded = [l for l in ent_refolded if l is not None]
        
        # for all proteins get the set that are refolded
        Refolded_ALL = all_refolded[0]
        for i in range(1, len(all_refolded)):
            Refolded_ALL = set(Refolded_ALL).intersection(all_refolded[i])
                
        print(Refolded_ALL)
        print(f'Number refolded across all conditions: {len(Refolded_ALL)}')
        refolded_outf = os.path.join(self.outpath, 'ALL_Refolded.csv')
        df = pd.DataFrame({'gene':list(Refolded_ALL)})
        df.to_csv(refolded_outf, index=False)
        print(f'SAVED: {refolded_outf}')


        # for ENT proteins get the set that are refolded
        Refolded_ENT = ent_refolded[0]
        for i in range(1, len(ent_refolded)):
            Refolded_ENT = set(Refolded_ENT).intersection(ent_refolded[i])
                
        print(Refolded_ENT)
        print(f'Number ENT refolded across all conditions: {len(Refolded_ENT)}')
        refolded_outf = os.path.join(self.outpath, 'ENT_Refolded.csv')
        df = pd.DataFrame({'gene':list(Refolded_ENT)})
        df.to_csv(refolded_outf, index=False)
        print(f'SAVED: {refolded_outf}')

        ## quality check. all those in Refolded_ENT should be in Refolded_ALL
        for g in Refolded_ENT:
            if g not in Refolded_ALL:
                raise ValueError(f'{g} was in Refolded_ENT but not Refolded_ALL')
       

    def process_file(self, index, file_path):
        """
        Process a single file to FDR correct peptide data and select proteins that are refoldable. 
        Returns both list of proteins that are refoldable and in the all_genes mask at the 50th SPA and COV thresholds and also entangled sets at same tresholds

        :param index: Index of the file in the list.
        :param file_path: Path to the file to be processed.
        """
        ## only examine the cyto-serum only condition as we want refoldable proteins in the absence of any external help such as chaperones
        buff = file_path.split('/')[-1].split('_')[0]
        timepoint = file_path.split('/')[-1].split('_')[1]
        print(f'Buffer: {buff} & timepoint: {timepoint}')
        if buff != 'C':
            return None, None
        if timepoint not in ['R5min', 'R2hr']:
            return None, None


        ## path to outfile containing all peptids
        outf = os.path.basename(file_path).replace('.csv', '_corrected.csv')
        outf = os.path.join(self.outpath, outf)

        ## Read the file in and do the FDR correction
        df = pd.read_csv(file_path)
        df = df[['Accession', 'proteinaseKsite', 'Peptide Sequence', 'PeptideRatio1', 'PeptidePValue1']]
        df['pvalue_raw'] = 10 ** (-1 * df['PeptidePValue1'])
            
        all_df = df.copy()
        all_df = self.apply_corrections(all_df)

        #df = df[~df["proteinaseKsite"].str.contains('-')]
        df = self.apply_corrections(df)
      
        # save all corrected PK peptides
        df.to_csv(outf, index=False)
        #print(f'SAVED: {outf} {len(df)}')

        # get genes where all their peptides are insignificant
        refolded_dfs = []
        for gene, gene_df in df.groupby('Accession'):
            # as long as the pvalue is insignificant. Refolded ratio does not matter if we cannot statistically say the ratio is not 1.
            refolded_gene_df = gene_df[gene_df['FDR_pvalue'] > self.alpha]

            # if both the refolded and original df have the same length then all peptides were refolded
            if len(refolded_gene_df) == len(gene_df):
                #print(f'\n{gene_df.to_string()}')
                print(f'{gene} is REFOLDED at {buff} and {timepoint}')
                print(gene_df.to_string())
                print(refolded_gene_df.to_string())
                print(all_df[all_df['Accession'] == gene].to_string())
                refolded_dfs += [gene_df]


        ## concactenate all refolded genes and get the list regardless of thresholds on spa or lip cov
        refolded_dfs = pd.concat(refolded_dfs)    
        refolded_genes = np.unique(refolded_dfs['Accession'].values)

        ## get the proper ent_genes file for this buffer and load it
        ent_f = [f for f in self.ent_files if f'_{buff}_' in f]
        print(f'ent_f: {ent_f}')
        ent_genes = np.loadtxt(ent_f[0], dtype=str)
        print(f'Number of ent genes loaded: {ent_genes.shape}')

        ## get the proper all_genes file for this buffer and load it
        all_f = [f for f in self.all_files if f'_{buff}_' in f]
        print(f'all_f: {all_f}')
        all_genes = np.loadtxt(all_f[0], dtype=str)
        print(f'Number of all genes loaded: {all_genes.shape}')
        
        all_refolded_genes = np.asarray([g for g in all_genes if g in refolded_genes])
        print(f'Refolded: {all_refolded_genes} {all_refolded_genes.shape}')
        ent_refolded_genes = np.asarray([g for g in ent_genes if g in refolded_genes])
        print(f'Ent + Refolded: {ent_refolded_genes} {ent_refolded_genes.shape}')

        return all_refolded_genes, ent_refolded_genes


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
    parser.add_argument("-e", "--ent_genes", type=str, required=True, help="path to list of ent genes to mask over the set of refolded to get the final set of genes")
    args = parser.parse_args()

    processor = LiPMSProcessor(args.lipms_files, args.threshold, args.outpath, args.alpha, args.ent_genes)
    processor.process_files()

    print('NORMAL TERMINATION')

if __name__ == "__main__":
    main()

