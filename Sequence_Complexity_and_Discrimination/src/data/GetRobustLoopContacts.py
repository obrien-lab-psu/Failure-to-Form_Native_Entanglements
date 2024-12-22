import sys, os, re, time, logging
import argparse
import numpy as np
import pandas as pd
from glob import glob
from scipy.stats import permutation_test, ttest_ind, false_discovery_control
pd.set_option('display.max_rows', 500)
np.set_printoptions(linewidth=500)
# Suppress the specific RuntimeWarning

class DataScout:
    """
    """

    def __init__(self, args):
        """
        Initializing the data scout to find your data
        """

        self.outpath = args.outpath
        self.slug = args.slug
        print(f'outpath: {self.outpath}')
        print(f'slug: {self.slug}')

        ## Find pvalue and OR files from Compare_OR_with_permutation
        #../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Sequence_Complexity_and_Discrimination/Compare_OR_with_permutation/EXP/C_50_p100000/FrequencyGeneratorOutput/OR_GT.csv
        self.OR_files = glob(os.path.join(self.slug, f'Sequence_Complexity_and_Discrimination/Compare_OR_with_permutation/*/*_50_p100000/FrequencyGeneratorOutput/OR_GT.csv'))
        print(f'# OR_files: {len(self.OR_files)}')

        self.pvalue_files = glob(os.path.join(self.slug, f'Sequence_Complexity_and_Discrimination/Compare_OR_with_permutation/*/*_50_p100000/FrequencyGeneratorOutput/OR_pvalues.csv'))
        print(f'# pvalue_files: {len(self.pvalue_files)}')

    #################################################################################################################
    # Function to apply FDR correction
    def apply_fdr_correction(self, alpha=0.05):

        SigPairs = {'buff':[], 'dataset':[], 'sig_contacts':[], 'sig_contacts_str':[]}
        for file in self.pvalue_files:
            print(file)
            outfile = file.replace('OR_pvalues', 'OR_FDR_pvalues')
            buff = file.split('/')[-3].split('_')[0]
            dataset = file.split('/')[-4]
            print(buff, dataset)
            SigPairs['buff'] += [buff]
            SigPairs['dataset'] += [dataset]

            df = pd.read_csv(file)
            matrix = df.iloc[:,:-1]
            print(matrix, matrix.shape)

            # Flatten the lower triangle of the matrix and keep non-zero values
            lower_triangle = matrix.where(np.tril(np.ones(matrix.shape), 0).astype(bool))
            p_values = lower_triangle.stack().values
            #print(f'p_values: {p_values}')
            p_values = np.where(p_values == 0, 0.000000000000001, p_values)
            #print(f'p_values: {p_values}')

            # Apply FDR correction
            corrected_pvals = false_discovery_control(p_values, method='bh')
            #print(f'corrected_pvals: {corrected_pvals}')

            # Reshape corrected p-values back to the lower triangle format
            corrected_matrix = pd.DataFrame(np.full(matrix.shape, np.nan), index=matrix.index, columns=matrix.columns)
            #corrected_matrix = pd.DataFrame(np.zeros(matrix.shape), index=matrix.index, columns=matrix.columns)
            corrected_matrix = corrected_matrix.where(np.triu(np.ones(corrected_matrix.shape), 0).astype(bool))  # Upper triangle should remain zero
            idx = 0
            for i in range(0, len(matrix)):
                for j in range(i + 1):
                    corrected_matrix.iat[i, j] = corrected_pvals[idx]
                    idx += 1
            
            corrected_matrix.set_index(df['AA'], inplace=True)
            print(f'corrected_matrix:\n{corrected_matrix}')
            corrected_matrix.to_csv(outfile)
            print(f'SAVED: {outfile}')

            ## get tje significant pairs
            sig_contacts = corrected_matrix.where(corrected_matrix<=0.05)

            # Find the row and column pairs where the elements are not NaN
            not_nan_pairs = sig_contacts.notna().stack()  # Stacked boolean DataFrame (True for not NaN)
            not_nan_pairs = not_nan_pairs[not_nan_pairs]  # Filter to only True values
            # Get the index as row and column label pairs
            sig_contacts = list(not_nan_pairs.index)
            sig_contacts_str = ', '.join(['-'.join(x) for x in sig_contacts])               
            print(f'sig_contacts:\n{sig_contacts} {sig_contacts_str}')
            SigPairs['sig_contacts'] += [sig_contacts]
            SigPairs['sig_contacts_str'] += [sig_contacts_str]
            
        return pd.DataFrame(SigPairs)
    #################################################################################################################

    #################################################################################################################
    #################################################################################################################
    def robust(self,df):

        ## Get the robust loop forming contacts for each datatype
        robust_results = {}
        for dataset, dataset_df in df.groupby('dataset'):
            print(dataset_df)
            sig_contacts = dataset_df['sig_contacts'].values
            for i, sig_list in enumerate(sig_contacts):
                #print(i, sig_list)

                if i == 0:
                    robust_contacts = set(sig_list)
                else:
                    robust_contacts = robust_contacts.intersection(sig_list)
            print(robust_contacts)
            robust_results[dataset] = robust_contacts
        
        ## Get the final dataframe with OR and pvalues for each of these datasets
        outdf = {'dataset':[], 'buff':[], 'contact':[], 'OR':[], 'FDR_pvalue':[]}
        for dataset, robust_contacts in robust_results.items():
            dataset_OR_files = [f for f in self.OR_files if dataset in f]
            #print(f'dataset {dataset} has {len(dataset_OR_files)} dataset_OR_files')
            dataset_pvalue_files = [f for f in self.pvalue_files if dataset in f]
            #print(f'dataset {dataset} has {len(dataset_pvalue_files)} dataset_pvalue_files')         

            for buff in ['C', 'CD', 'CG']:
                OR_file = [f for f in dataset_OR_files if f'/{buff}_' in f][0]
                pvalue_file = [f for f in dataset_pvalue_files if f'/{buff}_' in f][0]
                #print(buff, OR_file, pvalue_file)
                ORs = pd.read_csv(OR_file)
                ORs.set_index(ORs['AA'], inplace=True)
                pvals = pd.read_csv(pvalue_file)
                pvals.set_index(pvals['AA'], inplace=True)
                for contact in robust_contacts:
                    #print(contact)

                    if contact == ('C', 'C'): ## ignore covalent bonded loop closures
                        continue

                    OR = ORs.loc[contact[0], contact[1]]
                    pval = pvals.loc[contact[0], contact[1]]
                    #print(OR, pval)
                    outdf['dataset'] += [dataset]
                    outdf['buff'] += [buff]
                    outdf['contact'] += ['-'.join(contact)]
                    outdf['OR'] += [OR]
                    outdf['FDR_pvalue'] += [pval]
        outdf = pd.DataFrame(outdf)
        print(outdf)
        outfile = os.path.join(self.outpath, 'RobustLoopClosingContactEnrichment.csv')
        outdf.to_csv(outfile, index=False)
        print(f'SAVED: {outfile}')
    #################################################################################################################
#################################################################################################################

def main():
    """
    This script is meant to find the loop closing contacts that are robustly enriched across all three buffer systems in either essential or non-essenentail proteins.
    """

    # Parse the user supplied arguments
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-f", "--slug", type=str, required=True, help=f"path to slug")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="path to output directory. will be made if doesnt exist")
    args = parser.parse_args()


    # Initalize the FrequencyGenerator class object
    scout = DataScout(args)
    print(f'scout: {scout}')


    # FDR correct the Pvalue files
    SigPairs = scout.apply_fdr_correction()
    print(f'SigPairs:\n{SigPairs}')

    # Get robust pairs
    scout.robust(SigPairs)

start_time = time.time()
if __name__ == "__main__":
    main()
print(f'NORMAL TERMINATION: {time.time() - start_time}')
