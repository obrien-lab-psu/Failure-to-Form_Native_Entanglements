import pandas as pd
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
    def __init__(self, inpfile, outpath):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
        - resFeat_files (str): Path to residue feature files.
        - outpath (str): Path to the output directory.
        """
        self.inpfile = inpfile
        self.outpath = outpath

        self.buff = self.inpfile.split('/')[-1].split('_')[1]
        self.timepoint = self.inpfile.split('/')[-1].split('_')[2]

        if not os.path.exists(f'{self.outpath}'):
            os.makedirs(f'{self.outpath}')
            print(f'Made output directories {self.outpath}')

        self.data_path = os.path.join(f'{self.outpath}', 'Data/')
        if not os.path.exists(f'{self.data_path}'):
            os.makedirs(f'{self.data_path}')
            print(f'Made output directories {self.data_path}')


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
    ## parse the raw LiPMS file 
    def parse_protein_file(self, file_path):
        print(f'{"#"*50}\nParse LiPMS File')
        results = {}
        outfile = f'{self.data_path}{self.buff}_{self.timepoint}_raw_data.pkl'
        print(file_path, self.buff, self.timepoint, outfile)
        cov_data = {'Accession':[], 'coverage':[]}
        if not os.path.exists(outfile):
        
            with open(file_path, 'r') as file:
                lines = file.readlines()

            blue_headers = lines[0].split('\t')

            for line_i, line in enumerate(lines[1:]):
                line = line.split('\t')

                # Determine which region the line is in 
                blue = False
                orange = False
                gray = False

                if line[0] == '' and line[1] == '':
                    gray = True
                    #print(f'\n GRAY {line}')
                    results[uniprot][(peptide, PTM, pos)] += [line[2:-14]]
                    if line[2] == 'Checked':
                        orange_headers = line

                elif line[0] == '':
                    orange = True
                    #print(f'\n ORANGE {line}')
                    if line[1] != 'Checked':
                        peptide = line[3]
                        PTM = line[4]
                        pos = line[9]
                        results[uniprot][(peptide, PTM, pos)] = []
                        continue

                elif line[0] != '':
                    blue = True
                    #print(f'\n BLUE {line}')
                    uniprot = line[3]
                    cov = line[6]
                    cov_data['Accession'] += [uniprot]
                    cov_data['coverage'] += [cov]
                    if uniprot in results:
                        print(f'ERROR: {uniprot} already in results meaning there are dupicate lines')
                        quit()
                    else:
                        results[uniprot] = {}
                        continue 

            
            ### SAVE coverage data for all proteins in the dataset
            coverage_outfile = f'{self.data_path}{self.buff}_{self.timepoint}_coverage.csv'
            cov_data = pd.DataFrame(cov_data)
            print(cov_data)
            cov_data.to_csv(coverage_outfile, index=False)
            print(f'SAVED: {coverage_outfile}')


            ## Process raw peptide information 
            keys = ['Avg. m/z [Da]', 'Charge']
            abund_keys = [k for k in orange_headers if 'Abundances (Normalized)' in k]
            keys += abund_keys
        
            new_cols = {'Avg. m/z [Da]':'Avg. m/z [Da]', 'Charge':'Charge'} 
            mapping = {0:'N1', 1:'N2', 2:'N3', 3:'R1', 4:'R2', 5:'R3'}
            for i, k in enumerate(abund_keys):
                new_cols[k] = mapping[i]


            processed_results = {}
            for uniprot, peptides in results.items():
                processed_results[uniprot] = {}
                for pep, pep_data in peptides.items():

                    ## Get PK cut site info
                    pep = self.check_rk_conditions(pep)
                    #print(f'\n{uniprot} {pep}')

                    if len(pep_data) > 0:
                        pep_data = pd.DataFrame(pep_data[1:], columns=pep_data[0])
                        pep_data = pep_data[keys]
                        pep_data.rename(columns=new_cols, inplace=True)

                        # Convert columns to numeric, replacing blanks with NaN
                        columns_to_convert = ['N1', 'N2', 'N3', 'R1', 'R2', 'R3']
                        pep_data[columns_to_convert] = pep_data[columns_to_convert].replace('', np.nan).apply(pd.to_numeric)
                        pep_data['case'] = pep_data.apply(self.determine_case, axis=1)
                        #print(f'\n{pep_data.to_string()}')

                        processed_results[uniprot][pep] = pep_data

            # view results of processings (comment out if you dont want to)
            #for uniprot, peptides in processed_results.items():
            #    for pep, pep_data in peptides.items():
            #        print(f'\n{uniprot} {pep}')
            #        print(pep_data)

            ## Save files
            with open(outfile, 'wb') as fh:
                pickle.dump(processed_results, fh)
            print(f'SAVED: {outfile}')
        
        else:
            with open(outfile, 'rb') as fh:
                processed_results = pickle.load(fh)
            print(f'LOADED: {outfile}') 

        self.raw = processed_results


    #######################################################################################
    # Function to determine the case value for each row
    def determine_case(self, row):
        n_values = [row['N1'], row['N2'], row['N3']]
        r_values = [row['R1'], row['R2'], row['R3']]
        
        n_missing = sum(pd.isnull(n_values))
        r_missing = sum(pd.isnull(r_values))
        
        if n_missing == 0 and r_missing == 0:
            return 1
        elif (n_missing == 1 and r_missing == 0) or (n_missing == 0 and r_missing == 1):
            return 2
        elif n_missing == 3 and r_missing == 0:
            return 3
        elif r_missing == 3 and n_missing == 0:
            return 3
        else:
            return 0

    #######################################################################################
    ## check for Trypsin cut sites R and K and find any half-tryptic PK cutsites
    def check_rk_conditions(self, pep_info):
        input_string = pep_info[0]
        PTM = pep_info[1]
        pos = pep_info[2].split(' ')[1]
        PK = None

        # Use regular expression to extract the parts of the string
        match = re.match(r'\[(.)\]\.(.*)\.\[(.)\]', input_string)
        
        if not match:
            warnings.warn(f'Invalid input string: {input_string}')
            return (input_string, PTM, pos, (None, None), PK)

        A = match.group(1)
        B_to_D = match.group(2)
        D = B_to_D[-1]
        E = match.group(3)
        #print(A, B_to_D, D, E)
        
        # Check the conditions
        is_A_rk = A in ['R', 'K']
        is_D_rk = D in ['R', 'K']

        # Get PK cut site
        # Use regular expression to find all numbers in the string
        numbers = re.findall(r'\d+', pos)
        seq_pos = list(map(int, numbers))

        # Check if it is an N terminal peptide
        if A == '-':
            # Check if it is terminal tryptic or terminal half-tryptic
            if is_D_rk:
                return (input_string, PTM, pos, (is_A_rk, is_D_rk), PK)
            else:
                PK = f'{D}{seq_pos[-1]}'
                return (input_string, PTM, pos, (is_A_rk, is_D_rk), PK)

        # Check if it is an C terminal peptide
        if E == '-':
            # Check if it is terminal tryptic or terminal half-tryptic
            if is_A_rk:
                return (input_string, PTM, pos, (is_A_rk, is_D_rk), PK)
            else:
                PK = f'{B_to_D[0]}{seq_pos[0]}'
                return (input_string, PTM, pos, (is_A_rk, is_D_rk), PK)

        
        # Check if it is a full tryptic peptide
        if is_A_rk == True and is_D_rk == True:
            return (input_string, PTM, pos, (is_A_rk, is_D_rk), PK)

        
        # Check if it is a full PK peptide
        if is_A_rk == False and is_D_rk == False:
            return (input_string, PTM, pos, (is_A_rk, is_D_rk), PK)


        # Check if it is half tryptic
        if np.sum([is_A_rk, is_D_rk]) != 2:
            if is_A_rk == False:
                PK = B_to_D[0]
                PK = f'{PK}{seq_pos[0]}'

            if is_D_rk == False:
                PK = D
                PK = f'{PK}{seq_pos[-1]}'
            return (input_string, PTM, pos, (is_A_rk, is_D_rk), PK)

        
    ##############################################
    def calc_peptide_features(self,):
        """
        For each peptide fill missing values if the case is 1,2,3 (skip 0 as those are invalid peptides)
        Then calculate the abundance ratios for each peptide and the pvalue
        """

        print(f'{"#"*50}\nCalculating Peptide Features')
        outfile = f'{self.data_path}{self.buff}_{self.timepoint}_processed_peptide_features_data.pkl'
        print(self.buff, self.timepoint, outfile)

        if not os.path.exists(outfile):
            for uniprot, pep_dict in self.raw.items():
                for pep, pep_df in pep_dict.items():
                    #print('\n' ,uniprot, pep)
                    #print(pep_df.to_string())

                    # Step 0: Remove any rows where the column 'case' == 0
                    # if the df is empty move onto the next peptide
                    df = pep_df[pep_df['case'] != 0].copy()
                    if len(df) == 0:
                        self.raw[uniprot][pep] = df
                        continue

                    
                    # Step 1: For rows where 'case' == 3, fill NaN values in N1, N2, N3, R1, R2, R3
                    mask_case_3 = df['case'] == 3
                    columns_to_fill = ['N1', 'N2', 'N3', 'R1', 'R2', 'R3']
                    
                    for col in columns_to_fill:
                        df.loc[mask_case_3, col] = df.loc[mask_case_3, col].apply(
                            lambda x: np.random.normal(10000, 1000) if pd.isnull(x) else x)
                    
                    # Step 2: Create 'log2(R/N)' column
                    def calculate_log2_ratio(row):
                        R_values = row[['R1', 'R2', 'R3']].dropna()
                        N_values = row[['N1', 'N2', 'N3']].dropna()
                        if len(R_values) == 0 or len(N_values) == 0:
                            return np.nan
                        R_mean = R_values.mean()
                        N_mean = N_values.mean()
                        return np.log2(R_mean / N_mean)
                    
                    df['log2(R/N)'] = df.apply(calculate_log2_ratio, axis=1)
                    
                    # Step 3: Create 'pvalue' column using Welch's t-test
                    def calculate_ttest(row):
                        R_values = row[['R1', 'R2', 'R3']].dropna()
                        R_values = np.asarray(R_values.values, dtype=float)
                        N_values = row[['N1', 'N2', 'N3']].dropna()
                        N_values = np.asarray(N_values.values, dtype=float)
                        if len(R_values) < 2 or len(N_values) < 2:
                            return pd.Series([np.nan, np.nan], index=['t_stat', 'pvalue'])
                        t_stat, p_value = ttest_ind(R_values, N_values, equal_var=False)
                        return pd.Series([t_stat, p_value], index=['t_stat', 'pvalue'])

                    df[['t_stat', 'pvalue']] = df.apply(calculate_ttest, axis=1)
                    self.raw[uniprot][pep] = df


            ## Save files
            with open(outfile, 'wb') as fh:
                pickle.dump(self.raw, fh)
            print(f'SAVED: {outfile}')
            self.pep_features = self.raw.copy()

        else:
            with open(outfile, 'rb') as fh:
                self.pep_features = pickle.load(fh)
            print(f'LOADED: {outfile}') 


    ##############################################
    def merge_peptides(self, ):
        """
        For all the ions that can inform on the susceptibility at a given cut-site, the (ratio, P-value) pairs for those ions are collectively considered.

        If they all agree in direction (i.e., the signs of the t test statistics are all the same), then the oerall ratio for the cut-site is calculated by taking the median of the ratios of all ions that map to it, and the P-values are combined with Fishers method to provide an updated (ratio, P-value) for the cut-site.

        If there are two independent ions and they disagree (e.g., the ion is more abundant in the test condition in the 2+ charge state but more abundant in the control condition in the 3+ charge state), then a median is still calculated, but the P-value is set to 1, implying there is no confidence as to whether this cut-site was more susceptible in the test or control condition.

        These cut-sites are discounted from the tally of the total valid cut-sites.

        If there are three ions, then a majority rules heuristic is applied: the disagreeing ion is disregarded, and the (ratio, P-value)s are only combined for the agreeing ions. In practice, it is relatively rare for more than three ions to be mapped to the same cut-site, but where this occurs, if a majority (or all) of the ions agree in direction, they are combined. If there is a tie, the P-value is set to 1.

        When there are redundancies in peptides mapping back to the same cut site (Figure 2), the old data analysis discarded the data point if the log2(refolded/native) values had different signs. In FLiPPR, if the signs are different, we discard the minority sign; if there is an equal amount of disagreement, we discard that cut site
        """

        print(f'{"#"*50}\nMerging peptides')
        PK_binary_outfile = f'{self.data_path}{self.buff}_{self.timepoint}_merged_PK_peptides.pkl'
        nonPK_binary_outfile = f'{self.data_path}{self.buff}_{self.timepoint}_merged_nonPK_peptides.pkl'
        outfile = f'{self.data_path}{self.buff}_{self.timepoint}_FLiPPR_processed.csv'
        SPA_outfile = f'{self.data_path}{self.buff}_{self.timepoint}_SPA.csv'
        print(self.buff, self.timepoint)

        if not os.path.exists(outfile):
            spa_df = {'Accession':[], 'N_spa':[]}
            merged_data = {'Accession':[], 'Peptide Sequence':[], 'Peptide Residue Range':[], 'proteinaseKsite':[], 'PeptideRatio1':[], 'PeptidePValue1':[]}
            for uniprot, pep_dict in self.pep_features.items():
                PK_peptides = {}
                nonPK_peptides = {}
                spa = 0
                print(f'{uniprot} has {len(pep_dict)} unique peptides')
                for pep, pep_df in pep_dict.items():
                    print('\n' ,uniprot, pep)
                    print(pep_df.to_string())

                    if len(pep_df) != 0:
                        spa += np.nanmean(pep_df[['N1', 'N2', 'N3']].values, )
                    PK = pep[-1]
                    seq = pep[0]
                    res_range = pep[2]

                    if PK != None:
                        
                        ## Check that the PK site is not already in the PK_peptides dict and or a PK site that is close by
                        PK_resid = int(''.join(re.findall(r'\d+', PK)))
                        if len(PK_peptides) == 0:
                            PK_peptides[PK] = {'pep_df': [pep_df], 'seq':[seq], 'res_range':[res_range]}
                        else:
                            current_PK_sites = {int(''.join(re.findall(r'\d+', k))):k for k in PK_peptides.keys()}
                            print(f'current_PK_sites: {current_PK_sites}')
                            merged = False
                            for resid, PK_site in current_PK_sites.items():
                                diff = abs(PK_resid - resid)
                                print(PK_resid, resid, PK_site, diff)

                                # if the difference between the current PK cutsite (PK) and the one already in the dict (PK_site) is less than or equal to 1 merge them
                                if diff <= 1:
                                    merged = True
                                    print(f'{PK} will be merged with {PK_site}')
                                    PK_peptides[PK_site]['pep_df'] += [pep_df]
                                    PK_peptides[PK_site]['seq'] += [seq]
                                    PK_peptides[PK_site]['res_range'] += [res_range]
                            
                            # check if it wasnt merged add it to its own entry
                            if not merged:
                                print(f'Merged: {merged}')
                                PK_peptides[PK] = {'pep_df': [pep_df], 'seq':[seq], 'res_range':[res_range]}

                    else:
                        if seq not in nonPK_peptides:
                            nonPK_peptides[seq] = {'pep_df': [pep_df], 'res_range':[res_range]}
                        else:
                            nonPK_peptides[seq]['pep_df'] += [pep_df]
                            nonPK_peptides[seq]['res_range'] += [res_range]
               
                #if uniprot == 'P07395':
                #    quit()
                
                ### update the SPA output df
                spa_df['Accession'] += [uniprot]
                spa_df['N_spa'] += [spa]

                ### get the merged PK peptide ratios and pvalues
                print(f'\n### Get the merged PK peptide ratios and pvalues')
                print(f'{"#"*50}\nunique PK sites: {len(PK_peptides)}')
                for PK, PK_data in PK_peptides.items():
                    PK_df = pd.concat(PK_data['pep_df'])
                    PK_seq = ';'.join(PK_data['seq'])
                    PK_res_range = ';'.join(PK_data['res_range'])
                    print(f'\n{uniprot}\nPK: {PK}\n{PK_seq}\n{PK_res_range}\n{PK_df.to_string()}')
                    if len(PK_df) != 0:
                        median_log2_ratio, combined_pvalue = self.merge_ratios_pvalues(PK_df)
                        print(f'median_log2_ratio: {median_log2_ratio} | combined_pvalue: {combined_pvalue}')

                        if median_log2_ratio != False and combined_pvalue != False:
                            merged_data['Accession'] += [uniprot]
                            merged_data['Peptide Sequence'] += [PK_seq]
                            merged_data['Peptide Residue Range'] += [PK_res_range]
                            merged_data['proteinaseKsite'] += [PK]
                            merged_data['PeptideRatio1'] += [median_log2_ratio]
                            #merged_data['PeptidePValue1'] += [combined_pvalue]
                            merged_data['PeptidePValue1'] += [(-1)*np.log10(combined_pvalue)]


                ### get the merged non PK peptide ratios and pvalues
                print(f'### Get the merged non PK peptide ratios and pvalues')
                print(f'{"#"*50}\nunique NON PK sites: {len(nonPK_peptides)}')
                for nonPK, nonPK_data in nonPK_peptides.items():
                    nonPK_df = pd.concat(nonPK_data['pep_df'])
                    nonPK_res_range = nonPK_data['res_range'][0]

                    print(f'\n{uniprot}\nnonPK: {nonPK}\n{nonPK_res_range}\n{nonPK_df.to_string()}')

                    if len(nonPK_df) != 0:
                        median_log2_ratio, combined_pvalue = self.merge_ratios_pvalues(nonPK_df)
                        print(f'median_log2_ratio: {median_log2_ratio} | combined_pvalue: {combined_pvalue}')

                        if median_log2_ratio != False and combined_pvalue != False:
                            merged_data['Accession'] += [uniprot]
                            merged_data['Peptide Sequence'] += [nonPK]
                            merged_data['Peptide Residue Range'] += [nonPK_res_range]
                            merged_data['proteinaseKsite'] += [nonPK_res_range]
                            merged_data['PeptideRatio1'] += [median_log2_ratio]
                            #merged_data['PeptidePValue1'] += [combined_pvalue]
                            merged_data['PeptidePValue1'] += [(-1)*np.log10(combined_pvalue)]
                
                
            merged_data = pd.DataFrame(merged_data)
            spa_df = pd.DataFrame(spa_df)

            ### SAVE the per protein SPA data
            spa_df.to_csv(SPA_outfile, index=False)
            print(f'SAVED: {SPA_outfile}')

            ### SAVE the binary merged PK data
            with open(PK_binary_outfile, 'wb') as fh:
                pickle.dump(PK_peptides, fh)
            print(f'SAVED: {PK_binary_outfile}')

            ### SAVE the binary merged nonPK data
            with open(nonPK_binary_outfile, 'wb') as fh:
                pickle.dump(nonPK_peptides, fh)
            print(f'SAVED: {PK_binary_outfile}')

            ### SAVE the final FLiPPR csv file
            merged_data.to_csv(outfile, index=False)
            print(f'SAVED: {outfile}')
        else:
            print(f'Merged data exists: {outfile}')
        
    ##############################################
    def merge_ratios_pvalues(self, df):
        # Step 1: Determine if there is an unequal number of positive and negative values in 'log2(R/N)'
        log2_ratios = df['log2(R/N)']
        num_positive = sum(log2_ratios > 0)
        num_negative = sum(log2_ratios < 0)
       
        if len(log2_ratios) == 2:
            if num_positive == num_negative:
                return False, False
        
        # Keep rows in the majority
        if num_positive > num_negative:
            df = df[log2_ratios > 0]
        if num_positive < num_negative:
            df = df[log2_ratios < 0]
        
        # Step 2: Calculate the median of 'log2(R/N)' and combine p-values using Fisher's method
        median_log2_ratio = df['log2(R/N)'].median()
        combined_pvalue = combine_pvalues(df['pvalue'])[1]

        # if there is more than 3 ions mapping to the same cutsite and there is equal disagreement
        # set the pvalue to 1
        if len(log2_ratios) > 3:
            if num_positive == num_negative:
                combined_pvalue = 1
        
        # Step 3: Return the median log2(R/N) ratio and the combined p-value
        return median_log2_ratio, combined_pvalue


########################################################################################################
def main():
    """
    Main function to parse arguments and run the DataAnalysis class.
    """

    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-i", "--inpfile", type=str, required=True, help="Path to raw LiPMS file to process")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory")
    args = parser.parse_args()

    inpfile = args.inpfile
    outpath = args.outpath
    print(f'inpfile: {inpfile}')
    print(f'outpath: {outpath}')

    Analyzer = DataAnalysis(inpfile, outpath)
    print(Analyzer)

    # load and parse the text based raw LiPMS file
    # do some preprocessing to identify PK cut sites and the missing data cases (but do not fill them yet)
    Analyzer.parse_protein_file(inpfile)

    # Fill missing values and calculate abundance ratios and pvalues
    Analyzer.calc_peptide_features()

    # Merge all peptides then merge all PK cut sites
    Analyzer.merge_peptides()

    print('NORMAL TERMINATION')


if __name__ == "__main__":
    main()

