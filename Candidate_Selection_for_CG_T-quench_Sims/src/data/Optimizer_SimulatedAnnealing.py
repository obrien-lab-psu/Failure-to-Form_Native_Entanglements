import time, sys
import pickle
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
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import poisson, binom, fisher_exact, chi2, norm, kstest
import scipy.stats as st
from matplotlib.ticker import MultipleLocator
import math, random
import itertools
from tqdm import tqdm

pd.set_option('display.max_rows', 4000)
pd.options.mode.chained_assignment = None  # default='warn'

class DataAnalysis:
    """
    A class to handle the data analysis process including encoding, regression, and statistical tests.
    """

    def __init__(self, resFeat_files, outpath, ent_gene_list, nonRefold_gene_list, tag, buffer, spa, cov, reg_formula, random, n_groups, steps, C1, C2, C3, restart_path, beta, linearT):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
        - resFeat_files (str): Path to residue feature files.
        - outpath (str): Path to the output directory.
        - gene_lists (str): Path to gene lists to use.
        - tag (str): Tag for output filenames.
        - buffer (str): Buffer system to use.
        - spa (str): SPA threshold.
        - cov (str): LiPMS cov threshold.
        - reg_formula (str): Regression formula.
        """
        self.resFeat_files = resFeat_files
        self.outpath = outpath
        self.ent_gene_list = ent_gene_list
        self.nonRefold_gene_list = nonRefold_gene_list
        self.tag = tag
        self.buffer = buffer
        self.spa = spa
        self.cov = cov
        self.reg_formula = reg_formula
        self.steps = steps
        self.C1 = C1
        self.C2 = C2
        self.C3 = C3
        self.beta = beta
        self.restart_path = restart_path
        self.linearT = linearT
        self.data = {}
        #self.logger = self.setup_logging()
        #self.gene_list_files = glob.glob(self.gene_list)
        self.timepoint = 'Rall'
        if random == 'False':
            self.random = False
        if random == 'True':
            self.random = True
        print(random, type(random))
        self.n_groups = n_groups

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

    def regression(self, df, formula, genes):
        """
        Performs quasi-binomial regression analysis on the provided DataFrame.

        Parameters:
        - df (pd.DataFrame): DataFrame containing the data for regression.
        - formula (str): The formula specifying the regression model.

        Returns:
        - table_1_df (pd.DataFrame): DataFrame containing the regression results with p-values.
        """

        model = sm.GLM.from_formula(formula, family=sm.families.Binomial(), data=df)
        result = model.fit()

        ## recalculate the pvalue to add more digits as statsmodels truncates it to 0 if it is below 0.0001 for some reason.
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

        coef = table_df[table_df['var'] == 'region']['coef'].values[0]
        std = table_df[table_df['var'] == 'region']['std err'].values[0]

        # get size dist
        size_dist = self.prot_size[self.prot_size['gene'].isin(genes)]['prot_size'].values
        return table_df, float(coef), float(std), df[f'cut_{self.buffer}_Rall'].sum(), table_df[table_df['var'] == 'region'], size_dist

    def metrics(self, df, formula, genes, cutkey):

        ctable = pd.crosstab(df[cutkey], df['region'])
        res = fisher_exact(ctable)
        OR, pvalue = res.statistic, res.pvalue

        # get size dist
        size_dist = self.prot_size[self.prot_size['gene'].isin(genes)]['prot_size'].values

        return OR, pvalue, df[f'cut_{self.buffer}_Rall'].sum(), size_dist

    def load_data(self):
        """
        Loads the residue feature files and filters the data for analysis.
        """
        res_files = glob.glob(os.path.join(self.resFeat_files, '*'))
        print(f'Number of res_files: {len(res_files)}')

        self.prot_size = {'gene':[], 'prot_size':[]}
        self.data = []
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
            temp = pd.read_csv(gene_resFeat_file, sep='|')
            temp = temp[~temp['mapped_resid'].isnull()]
            self.data += [temp]

            self.prot_size['gene'] += [gene]
            self.prot_size['prot_size'] += [temp['uniprot_length'].values[0]]

        self.data = pd.concat(self.data)
        self.data = self.data[self.data['gene'].isin(self.reg_genes)]
        self.data = self.data[self.data['AA'] != 'NC']
        self.data = self.data[self.data['mapped_resid'].notna()]
        self.data = self.data[self.data['AA'].notna()]
        self.data = self.data.reset_index()
        print(f"Data loaded and filtered. Number of unique genes: {len(self.data['gene'].unique())}")
        self.prot_size = pd.DataFrame(self.prot_size)
        #print(f'self.prot_size:\n{self.prot_size}')
     

    def run(self):
        """
        Orchestrates the workflow by loading data, performing regression, and saving results.
        """
        start_time = time.time()

        #####################################################################################################
        ## Initialize all the data required for the optimizer
        # Get list of entangled genes
        ent_genes = np.loadtxt(self.ent_gene_list, dtype=str)
        self.ent_genes = ent_genes
        logging.info(f'Number of genes in entangled list: {len(ent_genes)}')

        # Get list of non-refolded genes
        nonRefold_genes = np.loadtxt(self.nonRefold_gene_list, dtype=str)
        self.nonRefold_genes = nonRefold_genes
        logging.info(f'Number of genes in nonRefold list: {len(nonRefold_genes)}')

        # Get list of ent and non-refolded
        self.reg_genes = list(set(ent_genes).intersection(nonRefold_genes))
        logging.info(f'reg_genes: {len(self.reg_genes)}')

        if len(self.reg_genes) == 0:
            logging.info(f"Gene list is empty: {self.gene_list} {len(self.reg_genes)}")
            quit()

        # Load data
        self.load_data()

        ## load reference size dist
        self.ref_sizes = self.prot_size['prot_size'].values
        print(f'ref_sizes:\n{self.ref_sizes}')

        # Define output files and get gene list
        reg_outfile = os.path.join(self.outpath, f"regression_results_{self.tag}_{self.buffer}_{self.timepoint}_spa{self.spa}_LiPMScov{self.cov}.csv")
        self.data = self.data[self.data['gene'].isin(self.reg_genes)]
        logging.info(self.data)

        # Quality check to ensure all columns are present in the df
        reg_vars = [v for v in self.reg_formula.split(' ') if v not in ['~', '+']]
        keys = []
        cut_key = ''
        for reg_var in reg_vars:
            if '*' in reg_var:
                reg_var = reg_var.split('*')
            else:
                reg_var = [reg_var]
            for v in reg_var:

                if 'cut_' in v:
                    cut_key = v

                if v not in self.data:
                    logging.info(f"Regression variable {reg_var} not in the DataFrame")
                    sys.exit()
                else:
                    keys += [v]

        logging.info(f'cut_key: {cut_key}')
        print(f'cut_key: {cut_key}')
        keys = list(set(keys))
        bin_keys = keys.copy()
        if 'AA' in keys:
            bin_keys.remove('AA')
        logging.info(f'bin_keys: {bin_keys}')

        # Encode boolean columns
        encoded_df = self.encode_boolean_columns(self.data, bin_keys)
        encoded_df = encoded_df[['gene']+reg_vars]
        logging.info(f'encoded_df:\n{encoded_df}')

        # Check if random was requested
        if self.random:
            num_cuts = sum(encoded_df[cut_key].values)
            logging.info(f'num_cuts: {num_cuts}')
            genes = encoded_df['gene'].values
            encoded_df = pd.DataFrame({
                cut_key: encoded_df[cut_key].sample(frac=1).reset_index(drop=True),
                'AA': encoded_df['AA'].sample(frac=1).reset_index(drop=True),
                'region': encoded_df['region'].sample(frac=1).reset_index(drop=True)
            })
            encoded_df.insert(0, 'gene', genes)
            logging.info(f'random encoded_df:\n{encoded_df}')
            num_cuts = sum(encoded_df[cut_key].values)
            logging.info(f'num_cuts: {num_cuts}')
           
        #####################################################################################################
        ## Get ref perresidue cut dist
        ref_perResCuts_df = get_per_res_cuts(encoded_df, cut_key)
        ref_perResCuts = ref_perResCuts_df['per_res_cuts'].values
        #print(f'ref_perResCuts: {ref_perResCuts}')


        #####################################################################################################
        ## if restart_path was not specified initialize random groups of genes
        ## do it randomly untill the groups ceof are in an assending order
        if self.restart_path != None:
            print(f'Restart path detected. Will load files from {self.restart_path} and start from there')
            groups = {}
            restart_dfs = {}
            for state in range(self.n_groups):
                reg_results_file = glob.glob(os.path.join(self.restart_path, f'State{state}_reg_results_{self.tag}_{self.buffer}_{self.timepoint}_spa{self.spa}_LiPMScov{self.cov}.csv'))
                final_gene_list_file = glob.glob(os.path.join(self.restart_path, f'State{state}_final_genelist_{self.tag}_{self.buffer}_{self.timepoint}_spa{self.spa}_LiPMScov{self.cov}.txt'))
                print(state, reg_results_file, final_gene_list_file)
                
                ## QC to ensure only 1 file is restarted from
                if len(reg_results_file) != 1:
                    raise ValueError(f'There is more than 1 reg_results file for state {state} to restart from!')
                if len(final_gene_list_file) != 1:
                    raise ValueError(f'There is more than 1 final_gene_list file for state {state} to restart from!')                   

                # load final set of genes to restart from for this state
                groups[state] = {}
                groups[state]['genes'] = list(np.loadtxt(final_gene_list_file[0], dtype=str))

                # load final df
                df = pd.read_csv(reg_results_file[0], sep='|')
                restart_dfs[state] = df
                last_row = df.tail(1)
                last_step = last_row['step'].values[0]
                print(f'Last step: {last_step}')
                print(last_row)

                # get the inital OR, number of cuts, and regression row
                _, coef, std, cuts, reg_row, size_dist = self.regression(encoded_df[encoded_df['gene'].isin(groups[state]['genes'])][reg_vars], self.reg_formula, groups[state]['genes'])
                state_size_dist_bootres = bootstrap((size_dist,) , np.mean)
                state_size_mean, state_size_lb, state_size_ub = np.mean(size_dist), state_size_dist_bootres.confidence_interval.low, state_size_dist_bootres.confidence_interval.high
                ks_stat = kstest(self.ref_sizes, size_dist).statistic
                E = -1*self.C1*coef + self.C2*(ks_stat)
                groups[state]['coef'] = [coef]
                groups[state]['std'] = [std]
                groups[state]['cuts'] = [cuts]
                groups[state]['reg_row'] = [reg_row]
                groups[state]['size_dist'] = [size_dist]
                groups[state]['psize_mean'] = [state_size_mean]
                groups[state]['psize_lb'] = [state_size_lb]
                groups[state]['psize_ub'] = [state_size_ub]
                groups[state]['step'] = [last_step]
                groups[state]['E'] = [E]
                groups[state]['ks_stat'] = [ks_stat]

        else:
            print(f'No restart detected. will use randomized starting states')
            groups = {}
            subgroups = create_unique_subgroups(np.arange(len(self.reg_genes)), self.n_groups)
            coefs = []
            last_step = -1
            for n, subgroup in enumerate(subgroups):
                #logging.info(n, len(subgroup), subgroup, type(subgroup))

                subgroup = np.array(subgroup, dtype=int)
                groups[n] = {}
                groups[n]['genes'] = [self.reg_genes[s] for s in subgroup]

                # get the inital OR, number of cuts, and regression row
                #_, coef, std, cuts, reg_row, size_dist = self.regression(encoded_df[encoded_df['gene'].isin(groups_temp[n]['genes'])][reg_vars], self.reg_formula, groups_temp[n]['genes'])
                OR, pvalue, cuts, size_dist = self.metrics(encoded_df[encoded_df['gene'].isin(groups[n]['genes'])][reg_vars], self.reg_formula, groups[n]['genes'], cut_key)

                state_size_dist_bootres = bootstrap((size_dist,) , np.mean)
                state_size_mean, state_size_lb, state_size_ub = np.mean(size_dist), state_size_dist_bootres.confidence_interval.low, state_size_dist_bootres.confidence_interval.high
                ks_stat_size = kstest(self.ref_sizes, size_dist).statistic

                perResCuts = ref_perResCuts_df[ref_perResCuts_df['gene'].isin(groups[n]['genes'])]['per_res_cuts'].values
                ks_stat_perResCuts = kstest(ref_perResCuts, perResCuts).statistic
                #E = -1*self.C1*np.log(OR) + self.C2*(ks_stat)
                E = -1*self.C1*np.log(OR) + self.C2*(ks_stat_size) + self.C3*(ks_stat_perResCuts)

                groups[n]['OR'] = [OR]
                groups[n]['pvalue'] = [pvalue]
                groups[n]['cuts'] = [cuts]
                groups[n]['size_dist'] = [size_dist]
                groups[n]['psize_mean'] = [state_size_mean]
                groups[n]['psize_lb'] = [state_size_lb]
                groups[n]['psize_ub'] = [state_size_ub]
                groups[n]['step'] = [last_step]
                groups[n]['E'] = [E]
                groups[n]['ks_stat_size'] = [ks_stat_size]
                groups[n]['ks_stat_perResCuts'] = [ks_stat_perResCuts]
                groups[n]['beta'] = [self.beta]

        ## log the starting information 
        for state in range(len(groups)):
            state_data = groups[state]
            state_OR = state_data['OR'][-1]
            state_pvalue = state_data['pvalue'][-1]
            state_cuts = state_data['cuts'][-1]
            num_genes = len(state_data['genes'])
            state_step = state_data['step'][-1]
            state_E = state_data['E'][-1]
            ks_stat_size = state_data['ks_stat_size'][-1]
            ks_stat_perResCuts = state_data['ks_stat_perResCuts'][-1]
            state_beta = state_data['beta'][-1]

            logging.info(f'STEP: {state_step} | State: {state} | OR: {state_OR} | pvalue: {state_pvalue} | cuts: {state_cuts} | num_genes: {num_genes} | E: {state_E} | ks_stat_size: {ks_stat_size} | ks_stat_perResCuts: {ks_stat_perResCuts} | beta: {state_beta}')
            print(f'STEP: {state_step} | State: {state} | OR: {state_OR} | pvalue: {state_pvalue} | cuts: {state_cuts} | num_genes: {num_genes} | E: {state_E} | ks_stat_size: {ks_stat_size} | ks_stat_perResCuts: {ks_stat_perResCuts} | beta: {state_beta}')

        #####################################################################################################
        #####################################################################################################
        ## Start optimizer and run for X steps using simulated annealing where 
        # the energy function for each state is defined as E(i) = -C1*OR(i) + C2*ks-test(ref_sizes, sizes))|
        # for each step we swap non-overlapping states. For exampline in a 5 state system
        # step 1 we swap 1-2, 3-4
        # step 2 we swap 2-3, 4-5
        # this constitutes one full montecarlo step
        # 
        # after every 1000 MC steps move to the next beta where beta ranges from (1/15 to 1000)
        #
        # The objective function is then M = exp(-beta*deltaE)
        # if deltaE is <=0 Accept the step
        # else
        #   get random float [0,1]
        #   M > random_float and M < 1 accept, else reject
        #
        # beta should start at 50 times the starting energy scale and decrease by half every 1000 steps till it reaches 1 (?)

        def swap_n(list1, list2, n):
            # Ensure n is not larger than the smaller list size
            n = min(n, len(list1), len(list2))
            list1 = np.array(list1)
            list2 = np.array(list2)
            
            # Randomly select n indices from both lists
            indices1 = random.sample(range(len(list1)), n)
            indices2 = random.sample(range(len(list2)), n)

            # Swap the elements at the selected indices
            for i in range(n):
                idx1 = indices1[i]
                idx2 = indices2[i]
                gene1 = list1[indices1[i]]
                gene2 = list2[indices2[i]]
                #print(i, idx1, idx2, gene1, gene2)
                list1[idx1] = gene2
                list2[idx2] = gene1

            return list1, list2

        ## create swapping scheme for each monte carlo step
        print(f'Making paris for MC swaps of {self.n_groups} groups in {self.steps} steps')
        pairs = [(i, i+1) for i in range(0, self.n_groups - 1 )]
        print(f'pairs: {pairs}')
      
        reps = self.steps
        beta = self.beta
        beta_i = 0
        if self.linearT == 'False':
            betas = np.linspace(beta, 1000, 75)
        elif self.linearT == 'True':
            T = 1/beta
            Ts = np.linspace(T, 0.001, 75)
            print(f'Ts: {Ts}')
            betas = 1/Ts
        print(f'betas: {betas} wiht start beta: {betas[beta_i]}')
        logging.info(f'Starting simulations with {reps} steps and beta = {self.beta}')
        
        for step in tqdm(range(last_step + 1, reps + last_step + 1)):
            #for step in range(last_step + 1, reps + last_step + 1):
            #logging.info(f'{"#"*100}\n{"#"*100}\nSTEP: {step}')
            #print(f'{"#"*100}\n{"#"*100}\nSTEP: {step}')

            ## For each pair in the pairs to test
            for pair in pairs:
                #print(f'{"#"*100}\npair: {pair}')
                
                # Get previous energy 
                Eold = groups[pair[0]]['E'][-1] + groups[pair[1]]['E'][-1]
                #print(f'Eold: {Eold}')

                # Get old state genes
                p0_genes = groups[pair[0]]['genes']
                p1_genes = groups[pair[1]]['genes']

                # Swap n genes 
                p0_genes_prime, p1_genes_prime = swap_n(p0_genes, p1_genes, 5)

                # get new regression info
                #_, p0_coef, p0_std, p0_cuts, p0_reg_row, p0_size_dist = self.regression(encoded_df[encoded_df['gene'].isin(p0_genes_prime)][reg_vars], self.reg_formula, p0_genes_prime)
                p0_OR, p0_pvalue, p0_cuts, p0_size_dist = self.metrics(encoded_df[encoded_df['gene'].isin(p0_genes_prime)][reg_vars], self.reg_formula, p0_genes_prime, cut_key)
                p0_ks_stat_size = kstest(self.ref_sizes, p0_size_dist).statistic
                p0_perResCuts = ref_perResCuts_df[ref_perResCuts_df['gene'].isin(p0_genes_prime)]['per_res_cuts'].values
                p0_ks_stat_perResCuts = kstest(ref_perResCuts, p0_perResCuts).statistic
                Ep0 = -1*self.C1*np.log(p0_OR) + self.C2*(p0_ks_stat_size) + self.C3*(p0_ks_stat_perResCuts)

                #_, p1_coef, p1_std, p1_cuts, p1_reg_row, p1_size_dist = self.regression(encoded_df[encoded_df['gene'].isin(p1_genes_prime)][reg_vars], self.reg_formula, p1_genes_prime)
                p1_OR, p1_pvalue, p1_cuts, p1_size_dist = self.metrics(encoded_df[encoded_df['gene'].isin(p1_genes_prime)][reg_vars], self.reg_formula, p1_genes_prime, cut_key)
                p1_ks_stat_size = kstest(self.ref_sizes, p1_size_dist).statistic
                p1_perResCuts = ref_perResCuts_df[ref_perResCuts_df['gene'].isin(p1_genes_prime)]['per_res_cuts'].values
                p1_ks_stat_perResCuts = kstest(ref_perResCuts, p1_perResCuts).statistic
                Ep1 = -1*self.C1*np.log(p1_OR) + self.C2*(p1_ks_stat_size) + self.C3*(p1_ks_stat_perResCuts)

                # Calculate new E and deltaE
                Enew = Ep0 + Ep1
                #print(f'Enew: {Enew}')
                deltaE = Enew - Eold
                #print(f'deltaE: {deltaE}')

                # Apply metropolis critera
                rand_float = random.uniform(0, 1)
                M = np.exp(-1*beta*deltaE)
                if deltaE <= 0:
                    accept_M = True
                else:
                    ## Apply metropolis critera
                    if M < 1 and M > rand_float:
                        accept_M = True
                    else:
                        accept_M = False
                #print(f'M: {M} with accept_M: {accept_M}')
                
                if accept_M:
                    groups[pair[0]]['genes'] = p0_genes_prime
                    groups[pair[1]]['genes'] = p1_genes_prime
            
            ## Get step summary after all pairs have been 
            #print(f'{"#"*100}\nStep summary for {step}')
            logstr = [f'{"#"*50}']
            for state, state_data in groups.items():
                state_genes = state_data['genes']
                #_, step_coef, step_std, step_cuts, step_reg_row, step_size_dist = self.regression(encoded_df[encoded_df['gene'].isin(state_genes)][reg_vars], self.reg_formula, state_genes)
                state_OR, state_pvalue, state_cuts, state_size_dist = self.metrics(encoded_df[encoded_df['gene'].isin(state_genes)][reg_vars], self.reg_formula, state_genes, cut_key)
                state_size_dist_bootres = bootstrap((state_size_dist,) , np.mean)
                state_size_mean, state_size_lb, state_size_ub = np.mean(state_size_dist), state_size_dist_bootres.confidence_interval.low, state_size_dist_bootres.confidence_interval.high

                ks_stat_size = kstest(self.ref_sizes, state_size_dist).statistic

                perResCuts = ref_perResCuts_df[ref_perResCuts_df['gene'].isin(state_genes)]['per_res_cuts'].values
                ks_stat_perResCuts = kstest(ref_perResCuts, perResCuts).statistic

                E = -1*self.C1*np.log(state_OR) + self.C2*(ks_stat_size) + self.C3*(ks_stat_perResCuts)
                #print(f'STEP: {step} | state {state} OR: {state_OR} pvalue: {state_pvalue} cuts: {state_cuts} size_mean: {state_size_mean} | ks_stat: {ks_stat} | E: {E}')
                logstr += [f'STEP: {step} | state {state} OR: {state_OR} pvalue: {state_pvalue} cuts: {state_cuts} size_mean: {state_size_mean} | ks_stat_size: {ks_stat_size} | ks_stat_perResCuts: {ks_stat_perResCuts} | E: {E} | beta: {beta}']

                groups[state]['OR'] += [state_OR]
                groups[state]['pvalue'] += [state_pvalue]
                groups[state]['cuts'] += [state_cuts]
                groups[state]['size_dist'] += [state_size_dist]
                groups[state]['psize_mean'] += [state_size_mean]
                groups[state]['psize_lb'] += [state_size_lb]
                groups[state]['psize_ub'] += [state_size_ub]
                groups[state]['step'] += [step]
                groups[state]['ks_stat_size'] += [ks_stat_size]
                groups[state]['ks_stat_perResCuts'] += [ks_stat_perResCuts]
                groups[state]['E'] += [E]
                groups[state]['beta'] += [beta]

 
            ## check ranks
            #old_ranks = sorted(range(len(old_coefs)), key=lambda i: old_coefs[i], reverse=True)
            #new_ranks = sorted(range(len(new_coefs)), key=lambda i: new_coefs[i], reverse=True)
            #rank_cond = new_ranks == old_ranks

            # logging.info status of step to log file
            if step % 100 == 0:
                logstr += [f'{"#"*50}']
                logstr = '\n'.join(logstr)
                logging.info(logstr)


            # update beta
            if step % 750 == 0 and step > 10:
                beta_i += 1
                if beta_i < len(betas):
                    beta = betas[beta_i]
                else:
                    beta = 1000
                logging.info(f'Beta update: {beta}')


        logging.info(f'{"#*50"}Simulation complete')
        logging.info(f'{"#*50"}Final state stats')
        for state in range(len(groups)):
            state_data = groups[state]
            state_OR = state_data['OR'][-1]
            state_pvalue = state_data['pvalue'][-1]
            state_cuts = state_data['cuts'][-1]
            state_size_dist = state_data['size_dist'][-1]
            state_size_mean = state_data['psize_mean'][-1]
            state_size_lb = state_data['psize_lb'][-1]
            state_size_ub = state_data['psize_ub'][-1]
            state_ks_stat_size = state_data['ks_stat_size'][-1]
            state_ks_stat_perResCuts = state_data['ks_stat_perResCuts'][-1]
            state_E = state_data['E'][-1]
            state_beta = state_data['beta'][-1]
            logging.info(f'State: {state} | OR: {state_OR} | pvalue: {state_pvalue} | cuts: {state_cuts} | state_size_mean: {state_size_mean:.0f} ({state_size_lb:.0f}, {state_size_ub:.0f}) | state_ks_stat_size: {state_ks_stat_size} | E: {state_E} | state_ks_stat_perResCuts: {state_ks_stat_perResCuts} | beta: {state_beta}')
            print(f'State: {state} | OR: {state_OR} | pvalue: {state_pvalue} | cuts: {state_cuts} | state_size_mean: {state_size_mean:.0f} ({state_size_lb:.0f}, {state_size_ub:.0f}) | state_ks_stat_size: {state_ks_stat_size} | E: {state_E} | state_ks_stat_perResCuts: {state_ks_stat_perResCuts} | beta: {state_beta}')
        #####################################################################################################

        #####################################################################################################
        ## Save final results
        dfs = []
        for state, state_data in groups.items():
            #logging.info(n, len(subgroup), subgroup, type(subgroup))

            # get the inital OR, number of cuts, and regression row
            _, coef, std, cuts, reg_row, size_dist = self.regression(encoded_df[encoded_df['gene'].isin(groups[state]['genes'])][reg_vars], self.reg_formula, groups[state]['genes'])
            state_size_dist_bootres = bootstrap((size_dist,) , np.mean)
            state_size_mean, state_size_lb, state_size_ub = np.mean(size_dist), state_size_dist_bootres.confidence_interval.low, state_size_dist_bootres.confidence_interval.high

            ks_stat_size = kstest(self.ref_sizes, size_dist).statistic

            perResCuts = ref_perResCuts_df[ref_perResCuts_df['gene'].isin(groups[state]['genes'])]['per_res_cuts'].values
            ks_stat_perResCuts = kstest(ref_perResCuts, perResCuts).statistic

            E = -1*self.C1*coef + self.C2*(ks_stat_size) + self.C3*(ks_stat_perResCuts)

            reg_row['state'] = state
            reg_row['beta'] = beta
            reg_row['cuts'] = cuts
            reg_row['ks_stat_size'] = ks_stat_size
            reg_row['ks_stat_perResCuts'] = ks_stat_perResCuts
            reg_row['E'] = E
            reg_row['psize_mean'] = state_size_mean
            reg_row['psize_lb'] = state_size_lb
            reg_row['psize_ub'] = state_size_ub
            reg_row['step'] = step
            reg_row['OR'] = np.exp(coef)
            reg_row['OR_lb'] = np.exp(reg_row['[0.025'].astype(float))
            reg_row['OR_ub'] = np.exp(reg_row['0.975]'].astype(float))
            reg_row['tag'] = self.tag
            reg_row['buff'] = self.buffer
            reg_row['timepoint'] = self.timepoint
            reg_row['spa'] = self.spa
            reg_row['cov'] = self.cov
            reg_row['n'] = len(state_data['genes'])
            dfs += [reg_row]

            ## save the final gene list for the state
            state_final_genelist_outfile = f'{self.outpath}State{state}_final_genelist_{self.tag}_{self.buffer}_{self.timepoint}_spa{self.spa}_LiPMScov{self.cov}.txt'
            logging.info(state_data['genes'])
            np.savetxt(state_final_genelist_outfile, list(state_data['genes']), fmt='%s')
            logging.info(f'SAVED: {state_final_genelist_outfile}')

            ## save the state data for this state
            state_df = {'state':[], 'step':[], 'OR':[], 'pvalue':[], 'cuts':[], 'psize_mean':[], 'psize_lb':[], 'psize_ub':[], 'ks_stat_size':[], 'E':[], 'ks_stat_perResCuts':[], 'beta':[]}
            for i in range(len(state_data['step'])):
                state_df['state'] += [state]
                state_df['step'] += [state_data['step'][i]]
                state_df['OR'] += [state_data['OR'][i]]
                state_df['pvalue'] += [state_data['pvalue'][i]]
                state_df['cuts'] += [state_data['cuts'][i]]
                state_df['psize_mean'] += [state_data['psize_mean'][i]]
                state_df['psize_lb'] += [state_data['psize_lb'][i]]
                state_df['psize_ub'] += [state_data['psize_ub'][i]]
                state_df['ks_stat_size'] += [state_data['ks_stat_size'][i]]
                state_df['ks_stat_perResCuts'] += [state_data['ks_stat_perResCuts'][i]]
                state_df['E'] += [state_data['E'][i]]
                state_df['beta'] += [state_data['beta'][i]]
            state_df = pd.DataFrame(state_df)
            #print(state_df)
            state_final_traj_outfile = f'{self.outpath}State{state}_final_traj_{self.tag}_{self.buffer}_{self.timepoint}_spa{self.spa}_LiPMScov{self.cov}.csv'
            state_df.to_csv(state_final_traj_outfile, index=False)
            logging.info(f'SAVED: {state_final_traj_outfile}')


        outdf = pd.concat(dfs)
        print(outdf)

        ## Save the final step regression data
        final_step_outfile = f'{self.outpath}Final_step_reg_{self.tag}_{self.buffer}_{self.timepoint}_spa{self.spa}_LiPMScov{self.cov}.csv'
        outdf.to_csv(final_step_outfile, index=False)
        print(f'SAVED: {final_step_outfile}')
        logging.info(f'SAVED: {final_step_outfile}')

        #####################################################################################################

        logging.info(f'NORMAL TERMINATION {time.time() - start_time}')

def create_unique_subgroups(array, m):
    # Shuffle the array
    np.random.shuffle(array)
    
    # Calculate the size of each subgroup
    # If you want each subgroup to have equal size:
    subgroup_size = len(array) // m
    
    # Initialize the list of subgroups
    subgroups = []
    
    # Split the array into subgroups
    for i in range(m):
        start_index = i * subgroup_size
        end_index = start_index + subgroup_size
        
        # Handle the case where the last subgroup might be larger
        if i == m - 1:
            subgroups.append(array[start_index:])
        else:
            subgroups.append(array[start_index:end_index])
    
    return subgroups

def is_decending(array):
    # Iterate through the array and compare adjacent elements
    for i in range(len(array) - 1):
        if array[i] < array[i + 1]:
            return False
    return True

def get_per_res_cuts(df, cutkey):
    per_res_cuts_df = {'gene':[], 'per_res_cuts':[]}
    for gene, gene_df in df.groupby('gene'):
        per_res_cuts_df['gene'] += [gene]
        per_res_cuts_df['per_res_cuts'] += [np.sum(gene_df[cutkey])/len(gene_df)]
    per_res_cuts_df = pd.DataFrame(per_res_cuts_df)
    print(per_res_cuts_df)
    return per_res_cuts_df

def main():
    """
    Main function to parse arguments and run the DataAnalysis class.
    """
    script_name = f'Optimizer_SimulatedAnnealing'
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-f", "--resFeat_files", type=str, required=True, help="Path to residue feature files")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("--ent_gene_list", type=str, required=True, help="Path to gene list to use of entangled proteins")
    parser.add_argument("--nonRefold_gene_list", type=str, required=True, help="Path to gene list to use of non-refoldable proteins")
    parser.add_argument("-t", "--tag", type=str, required=True, help="Tag for output filenames")
    parser.add_argument("-b", "--buffer", type=str, required=True, help="Buffer system to use")
    parser.add_argument("-s", "--spa", type=str, required=True, help="SPA threshold")
    parser.add_argument("-c", "--cov", type=str, required=True, help="LiPMS cov threshold")
    parser.add_argument("-n", "--n_groups", type=int, required=True, help="number of groups to optimize")
    parser.add_argument("-r", "--reg_formula", type=str, required=True, help="Regression formula")
    parser.add_argument("--random", type=str, required=True, help="Randomize dataset")
    parser.add_argument("-l", "--log", type=str, required=True, help="Path to logging file")
    parser.add_argument("--restart_path", type=str, required=False, help="Path to a folder containing files to restart from")
    parser.add_argument("--steps", type=int, required=True, help="Number of steps to run")
    parser.add_argument("-C1", type=float, required=True, help="C1 coefficient for optimization function")
    parser.add_argument("-C2", type=float, required=True, help="C2 coefficient for optimization function")
    parser.add_argument("-C3", type=float, required=True, help="C3 coefficient for optimization function")
    parser.add_argument("-beta", type=float, required=True, help="Starting beta. If >= 1000 then no temperature quenching is done")
    parser.add_argument("-linearT", type=str, required=True, help="use a linear T scale instead of a linear beta value")
    args = parser.parse_args()

    # Setup logging configuration
    logging.basicConfig(filename=args.log, level=logging.INFO, format='%(asctime)s %(message)s')
    logging.info(f'{"#"*100}\nNEW RUN {script_name}')

    analysis = DataAnalysis(
        resFeat_files=args.resFeat_files,
        outpath=args.outpath,
        ent_gene_list=args.ent_gene_list,
        nonRefold_gene_list=args.nonRefold_gene_list,
        tag=args.tag,
        buffer=args.buffer,
        spa=args.spa,
        cov=args.cov,
        reg_formula=args.reg_formula,
        random=args.random,
        n_groups=args.n_groups,
        steps=args.steps,
        C1=args.C1,
        C2=args.C2,
        C3=args.C3,
        restart_path=args.restart_path,
        beta=args.beta,
        linearT=args.linearT)
    analysis.run()
    print('NORMAL TERMINATION')
if __name__ == "__main__":
    main()

