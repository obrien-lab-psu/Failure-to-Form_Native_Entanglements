import time, sys
import multiprocessing as mp
from scipy.stats import bootstrap
import logging
import argparse
import glob
import numpy as np
import pandas as pd
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

    def __init__(self, outpath):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
        - outpath (str): Path to the output directory.
        """
        self.outpath = outpath
        self.data = {}

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


    def load_data(self, dataset):
        """
        Loads the residue feature files and filters the data for analysis.
        """
        res_files = glob.glob(f'../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gen_proteome_features_{dataset}/res_features_lib/*.csv')
        print(f'Number of {dataset} res_files: {len(res_files)}')
        for i, gene_resFeat_file in enumerate(res_files):
            #gene_resFeat = [f for f in res_files if gene in f]

            #gene_resFeat_file = gene_resFeat[0]
            #print(f'gene_resFeat_file: {gene_resFeat_file} {i}')
            if len(self.data) == 0:
                self.data = pd.read_csv(gene_resFeat_file, sep='|')
            else:
                self.data = pd.concat((self.data, pd.read_csv(gene_resFeat_file, sep='|')))

        #self.data = self.data[self.data['gene'].isin(self.reg_genes)]
        self.data = self.data[self.data['AA'] != 'NC']
        self.data = self.data[self.data['mapped_resid'].notna()]
        self.data = self.data[self.data['AA'].notna()]
        self.data = self.data.reset_index()
        print(f"Data loaded and filtered. Number of unique genes: {len(self.data['gene'].unique())}")

#########################################################################################################
    def Population_anal(self):
        """
        Calculate the different populations of proteins in the dataset
        (A) Proteins with cuts only in the entangled region
        (B) Proteins with cuts only in the non-entangled region
        (C) Proteins with cuts in both regions but a greater proportion in the entangled regions 
        (D) Proteins with cuts in both regions but a greater proportion in the non-entangled regions 
        (E) Proteins with cuts in both regions and an ~ equal proportion in both regions 
        (F) Proteins with at least 1 cut site in the entangled region = A + C
        (G) Proteins with a greater proportion in the entangled regions 
        (H) Proteins with a greater proportion in the non-entangled regions 

        Calculate the average fraction of cuts in each region for each population 
        """
        start_time = time.time()

        ### Calculate the fraction of proteins with atleast 1 cut in the entangled region
        outfile_csv = os.path.join(self.outpath, 'FractionGenesCutinEntR.csv')
        tags = ['all_genes',
                    'ent_genes',
                    'essential_ent_genes',
                    'essential_genes',
                    'nonessential_ent_genes',
                    'nonessential_genes']

        tags = ['ent_genes',
                'essential_ent_genes',
                'nonessential_ent_genes']

        if not os.path.exists(outfile_csv):
            stats_df = {'dataset':[], 'tag':[], 'buff':[], 'spa':[], 'PopClass':[], 'mean':[], 'lb':[], 'ub':[], 'n':[], 'FracEntR':[], 'FracEntR_lb':[], 'FracEntR_ub':[], 'FracNonEntR':[], 'FracNonEntR_lb':[], 'FracNonEntR_ub':[]} 
            #for dataset in ['EXP', 'AF']:
            for dataset in ['EXP']:
                self.load_data(dataset)

                for tag in tags:

                    for buff in ['C', 'CD', 'CG']:
                        loc_df = self.data[['gene', 'region', f'cut_{buff}_Rall']]

                        for spa in [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]:
                            
                            A, B, C, D, E, F, G, H = 0, 0, 0, 0, 0, 0, 0, 0
                            AFracEntR, BFracEntR, CFracEntR, DFracEntR, EFracEntR, FFracEntR, GFracEntR, HFracEntR = [], [], [], [], [], [], [], []
                            AFracNonEntR, BFracNonEntR, CFracNonEntR, DFracNonEntR, EFracNonEntR, FFracNonEntR, GFracNonEntR, HFracNonEntR = [], [], [], [], [], [], [], []
                            
                            
                            gene_list = f'../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/{dataset}/{dataset}_0.6g_{buff}_Rall_spa{spa}_LiPMScov50_{tag}.txt'
                            genes = np.loadtxt(gene_list, dtype=str)
                            print(f'gene_list: {gene_list} {len(genes)}')
                            N = len(genes)
                            # get the binary array of genes that have cuts in their entangled regions

                            for gene, gene_df in loc_df.groupby('gene'):
                                if gene in genes:
                                    
                                    entR = gene_df[gene_df['region'] == 1]  
                                    entR_n = len(entR)  
                                    entR_cut = np.sum(entR[f'cut_{buff}_Rall'].values)

                                    nonentR = gene_df[gene_df['region'] == 0]    
                                    nonentR_n = len(nonentR)  
                                    nonentR_cut = np.sum(nonentR[f'cut_{buff}_Rall'].values)

                                    prop_entR_cuts = (entR_cut/entR_n)
                                    prop_nonentR_cuts = (nonentR_cut/nonentR_n)
                                    #print(gene, entR_n, entR_cut, nonentR_n, nonentR_cut, prop_entR_cuts, prop_nonentR_cuts)
                                    
                                    #(A) Proteins with cuts only in the entangled region
                                    if entR_cut != 0 and nonentR_cut == 0:
                                        A += 1
                                        AFracEntR += [prop_entR_cuts]
                                        AFracNonEntR += [prop_nonentR_cuts]

                                    #(B) Proteins with cuts only in the non-entangled region
                                    if entR_cut == 0 and nonentR_cut != 0:
                                        B += 1
                                        BFracEntR += [prop_entR_cuts]
                                        BFracNonEntR += [prop_nonentR_cuts]

                                    #(C) Proteins with cuts in both regions but a greater proportion in the entangled regions 
                                    if entR_cut != 0 and nonentR_cut != 0:
                                        if prop_entR_cuts > prop_nonentR_cuts:
                                            C += 1      
                                            CFracEntR += [prop_entR_cuts]
                                            CFracNonEntR += [prop_nonentR_cuts]

                                    #(D) Proteins with cuts in both regions but a greater proportion in the non-entangled regions
                                    if entR_cut != 0 and nonentR_cut != 0:
                                        if prop_entR_cuts < prop_nonentR_cuts:
                                            D += 1        
                                            DFracEntR += [prop_entR_cuts]
                                            DFracNonEntR += [prop_nonentR_cuts]

                                    #(E) Proteins with cuts in both regions and an ~ equal proportion in both regions
                                    if entR_cut != 0 and nonentR_cut != 0:
                                        if prop_entR_cuts == prop_nonentR_cuts:
                                            E += 1    
                                            EFracEntR += [prop_entR_cuts]
                                            EFracNonEntR += [prop_nonentR_cuts]

                                    #(F) Proteins with at least 1 cut site in the entangled region = A + C  
                                    if entR_cut != 0:
                                        F += 1 
                                        FFracEntR += [prop_entR_cuts]
                                        FFracNonEntR += [prop_nonentR_cuts]

                                    #(G) Proteins with a greater proportion in the entangled regions 
                                    if prop_entR_cuts > prop_nonentR_cuts:
                                        G += 1
                                        GFracEntR += [prop_entR_cuts]
                                        GFracNonEntR += [prop_nonentR_cuts]
                            
                                    #(H) Proteins with a greater proportion in the non-entangled regions 
                                    if prop_entR_cuts < prop_nonentR_cuts:
                                        H += 1
                                        HFracEntR += [prop_entR_cuts]
                                        HFracNonEntR += [prop_nonentR_cuts]

                            #print(f'A: {A} | B: {B} | C: {C} | D: {D} | E: {E} | F: {F} | N: {N}')

                            # get the stats for each class and add to the dataframe
                            #  stats_df = {'dataset':[], 'tag':[], 'buff':[], 'spa':[], 'PopClass':[], 'mean':[], 'lb':[], 'ub':[], 'n':[], 
                            # 'FracEntR':[], 'FracEntR_lb':[], 'FracEntR_ub':[], 'FracNonEntR':[], 'FracNonEntR_lb':[], 'FracNonEntR_ub':[]} 
                            for count, FracEntR, FracNonEntR, label in [(A, AFracEntR, AFracNonEntR, 'A'), (B, BFracEntR, BFracNonEntR, 'B'), (C, CFracEntR, CFracNonEntR, 'C'), 
                                                                        (D, DFracEntR, DFracNonEntR, 'D'), (E, EFracEntR, EFracNonEntR, 'E'), (F, FFracEntR, FFracNonEntR, 'F'),
                                                                        (G, GFracEntR, GFracNonEntR, 'G'), (H, HFracEntR, HFracNonEntR, 'H')]:
                                arr = [1]*count + [0]*(N-count)
                                median, mean, std, lb, ub = get_stats(arr, n_resamples=10000)
                                print(f'{label} Fraction of genes cut in entR: {mean} ({lb}, {ub}) (n={len(arr)}) {buff} {tag}')  
                                stats_df['dataset'] += [dataset]
                                stats_df['tag'] += [tag]
                                stats_df['buff'] += [buff]
                                stats_df['spa'] += [spa]
                                stats_df['PopClass'] += [label]                            
                                stats_df['mean'] += [mean]
                                stats_df['lb'] += [lb]
                                stats_df['ub'] += [ub]
                                stats_df['n'] += [len(arr)]

                                median, mean, std, lb, ub = get_stats(FracEntR, n_resamples=10000)
                                print(f'{label} Fraction cuts in entR: {mean} ({lb}, {ub}) (n={len(arr)}) {buff} {tag}')                                 
                                stats_df['FracEntR'] += [mean]
                                stats_df['FracEntR_lb'] += [lb]
                                stats_df['FracEntR_ub'] += [ub]

                                median, mean, std, lb, ub = get_stats(FracNonEntR, n_resamples=10000)
                                print(f'{label} Fraction cuts in nonentR: {mean} ({lb}, {ub}) (n={len(arr)}) {buff} {tag}')                                 
                                stats_df['FracNonEntR'] += [mean]
                                stats_df['FracNonEntR_lb'] += [lb]
                                stats_df['FracNonEntR_ub'] += [ub]

            stats_df = pd.DataFrame(stats_df)
            print(stats_df)
            outfile_csv = os.path.join(self.outpath, 'FractionGenesCutinEntR.csv')
            stats_df.to_csv(outfile_csv, index=False)
            print(f'SAVED: {outfile_csv}')
        else:
            stats_df = pd.read_csv(outfile_csv)
            print(f'LOADED: {outfile_csv}')
            print(stats_df)


        #############################################################
        ### Plot type v2
        # Create a figure with 3 rows and 1 column
        tags = [('essential_ent_genes', 'Essential'),
                ('nonessential_ent_genes', 'Non-essential')]
        colors = {'A':'orange', 'B':'blue', 'C':'green', 'D':'red', 'F':'black', 'G':'purple', 'H':'grey'}
        labels = {'A':'only entR cuts', 'B':'only Non-entR cuts', 'C':'both cut (more entR)', 'D':'both cut (more Non-entR)', 'F':'at least 1 cut entR', 'G':'more entR', 'H':'more non-entR'}
        for dataset, dataset_df in stats_df.groupby('dataset'):
            print(dataset_df)
            fig1, axs1 = plt.subplots(8, 3, figsize=(10, 12)) # figure for fraction of genes in class
            outpng1 = os.path.join(self.outpath, f'PopulationAnal_{dataset}_v2.1.png')

            fig2, axs2 = plt.subplots(8, 3, figsize=(10, 12)) # figure for fraction of genes in class
            outpng2 = os.path.join(self.outpath, f'PopulationAnal_{dataset}_FractionCutEntR_v2.1.png')   

            fig3, axs3 = plt.subplots(8, 3, figsize=(10, 12)) # figure for fraction of genes in class
            outpng3 = os.path.join(self.outpath, f'PopulationAnal_{dataset}_FractionCutNonEntR_v2.1.png')     

            for tagi, (tag, ylabel) in enumerate(tags):
                tag_df = dataset_df[dataset_df['tag'] == tag]

                for i, (buff, buff_df) in enumerate(tag_df.groupby('buff')):
                    #print(buff_df)
                    for j, (PopClass, PopClass_df) in enumerate(buff_df.groupby('PopClass')):
                        PopClass_df = PopClass_df[PopClass_df['spa'] == 50]
                        print(PopClass_df)
                        if PopClass in ['F']:
                            rowi = tagi 
                        elif PopClass in ['A', 'B']:
                            rowi = tagi + 4
                        elif PopClass in ['C', 'D']:
                            rowi = tagi + 6
                        elif PopClass in ['G', 'H']:
                            rowi = tagi + 2
                        else:
                            continue
                        
                        
                        yerr_lower = PopClass_df['mean'] - PopClass_df['lb']
                        yerr_upper = PopClass_df['ub'] - PopClass_df['mean']
                        yerr = [yerr_lower, yerr_upper]

                        axs1[rowi, i].errorbar(PopClass_df['spa'], PopClass_df['mean'], yerr=yerr, fmt='o',  capsize=5,  markersize=5, label=labels[PopClass], color=colors[PopClass])
                        axs1[rowi, i].set_ylabel(f'{ylabel}')
                        axs1[rowi, i].set_xlabel('spa')
                        axs1[rowi, i].set_title(f'{buff}')

                        if rowi in [0, 1, 2, 3]:
                            axs1[rowi, i].set_yticks(np.arange(-1, 1.25, 0.25))  # Add y-axis tick marks every 0.25
                            axs1[rowi, i].set_ylim(0,1)
                        if rowi in [4, 5]:
                            axs1[rowi, i].set_yticks(np.arange(-1, 1.25, 0.25))  # Add y-axis tick marks every 0.25
                            axs1[rowi, i].set_ylim(0,0.35)
                        if rowi in [6, 7]:
                            axs1[rowi, i].set_yticks(np.arange(-1, 1.25, 0.25))  # Add y-axis tick marks every 0.25
                            axs1[rowi, i].set_ylim(0,0.6)    
                        axs1[rowi, i].grid(True)  # Add grid lines to each subplot
            

                        yerr_lower = PopClass_df['FracEntR'] - PopClass_df['FracEntR_lb']
                        yerr_upper = PopClass_df['FracEntR_ub'] - PopClass_df['FracEntR']
                        yerr = [yerr_lower, yerr_upper]

                        axs2[rowi, i].errorbar(PopClass_df['spa'], PopClass_df['FracEntR'], yerr=yerr, fmt='o',  capsize=5,  markersize=5, label=labels[PopClass], color=colors[PopClass])
                        axs2[rowi, i].set_ylim(0,0.05)
                        axs2[rowi, i].set_ylabel(f'{ylabel}')
                        axs2[rowi, i].set_xlabel('spa')
                        axs2[rowi, i].set_title(f'{buff}')


                        yerr_lower = PopClass_df['FracNonEntR'] - PopClass_df['FracNonEntR_lb']
                        yerr_upper = PopClass_df['FracNonEntR_ub'] - PopClass_df['FracNonEntR']
                        yerr = [yerr_lower, yerr_upper]

                        axs3[rowi, i].errorbar(PopClass_df['spa'], PopClass_df['FracNonEntR'], yerr=yerr, fmt='o',  capsize=5,  markersize=5, label=labels[PopClass], color=colors[PopClass])
                        axs3[rowi, i].set_ylim(0,0.05)
                        axs3[rowi, i].set_ylabel(f'{ylabel}')
                        axs3[rowi, i].set_xlabel('spa')
                        axs3[rowi, i].set_title(f'{buff}')
                                                

            axs1[0,2].legend(loc='upper left', bbox_to_anchor=(1, 1))
            axs2[0,2].legend(loc='upper left', bbox_to_anchor=(1, 1))
            axs3[0,2].legend(loc='upper left', bbox_to_anchor=(1, 1))

            axs1[2,2].legend(loc='upper left', bbox_to_anchor=(1, 1))
            axs2[2,2].legend(loc='upper left', bbox_to_anchor=(1, 1))
            axs3[2,2].legend(loc='upper left', bbox_to_anchor=(1, 1))

            axs1[4,2].legend(loc='upper left', bbox_to_anchor=(1, 1))
            axs2[4,2].legend(loc='upper left', bbox_to_anchor=(1, 1))
            axs3[4,2].legend(loc='upper left', bbox_to_anchor=(1, 1))  

            axs1[6,2].legend(loc='upper left', bbox_to_anchor=(1, 1))
            axs2[6,2].legend(loc='upper left', bbox_to_anchor=(1, 1))
            axs3[6,2].legend(loc='upper left', bbox_to_anchor=(1, 1))  

            # Adjust the layout to prevent clipping of the legend
            fig1.suptitle('Fraction of proteins in class')
            fig1.tight_layout()
            fig1.savefig(outpng1)
            print(f'SAVED: {outpng1}')

            # Adjust the layout to prevent clipping of the legend
            fig2.suptitle('Fraction entR cut')
            fig2.tight_layout()
            fig2.savefig(outpng2)
            print(f'SAVED: {outpng2}')

            # Adjust the layout to prevent clipping of the legend
            fig3.suptitle('Fraction non-entR cut')
            fig3.tight_layout()
            fig3.savefig(outpng3)
            print(f'SAVED: {outpng3}')    
        print('NORMAL TERMINATION')
        #############################################################

        #############################################################
        ### Plot type v3
        # Create a figure with 3 rows and 1 column
        tags = [('essential_ent_genes', 'Essential'),
                ('nonessential_ent_genes', 'Non-essential')]
        colors = {'A':'orange', 'B':'blue', 'C':'green', 'D':'red', 'F':'black', 'G':'purple', 'H':'grey'}
        labels = {'A':'only entR cuts', 'B':'only Non-entR cuts', 'C':'both cut (more entR)', 'D':'both cut (more Non-entR)', 'F':'at least 1 cut entR', 'G':'more entR', 'H':'more non-entR'}
        for dataset, dataset_df in stats_df.groupby('dataset'):
            print(dataset_df)
            fig1, axs1 = plt.subplots(4, 3, figsize=(10, 10)) # figure for fraction of genes in class
            outpng1 = os.path.join(self.outpath, f'PopulationAnal_{dataset}_v3.png') 

            for tagi, (tag, ylabel) in enumerate(tags):
                tag_df = dataset_df[dataset_df['tag'] == tag]

                for i, (buff, buff_df) in enumerate(tag_df.groupby('buff')):
                    #print(buff_df)
                    for j, (PopClass, PopClass_df) in enumerate(buff_df.groupby('PopClass')):
                        #print(PopClass_df)
                        if PopClass in ['F']:
                            rowi = tagi 
                        elif PopClass in ['G', 'H']:
                            rowi = tagi + 2
                        else:
                            continue
                        
                        
                        yerr_lower = PopClass_df['mean'] - PopClass_df['lb']
                        yerr_upper = PopClass_df['ub'] - PopClass_df['mean']
                        yerr = [yerr_lower, yerr_upper]

                        axs1[rowi, i].errorbar(PopClass_df['spa'], PopClass_df['mean'], yerr=yerr, fmt='o',  capsize=5,  markersize=5, label=labels[PopClass], color=colors[PopClass])
                        axs1[rowi, i].set_ylabel(f'{ylabel}')
                        axs1[rowi, i].set_xlabel('spa')
                        axs1[rowi, i].set_title(f'{buff}')
                        axs1[rowi, i].set_ylim(0,1)
                        #if rowi in [2, 3]:
                        #    axs1[rowi, i].set_ylim(0,0.35)
                        #if rowi in [4, 5]:
                        #    axs1[rowi, i].set_ylim(0,0.6)                


            axs1[0,2].legend(loc='upper left', bbox_to_anchor=(1, 1))

            axs1[2,2].legend(loc='upper left', bbox_to_anchor=(1, 1))


            # Adjust the layout to prevent clipping of the legend
            fig1.suptitle('Fraction of proteins in class')
            fig1.tight_layout()
            fig1.savefig(outpng1)
            print(f'SAVED: {outpng1}')  
        print('NORMAL TERMINATION')
#######################################################################################


#######################################################################################
def get_stats(arr, n_resamples=100, metric='mean'):
    """
    Get the <> and 95% ci for a stats array
    """
    if metric == 'mean':
        mean = np.mean(arr)
        std = np.std(arr)
        (lb, ub) = bootstrap(arr, n_resamples=100)
        median = np.median(arr)
        return (median, mean, std, lb, ub)
    
    elif metric == 'sum':
        sum = np.sum(arr)
        std = np.std(arr)
        (lb, ub) = bootstrap(arr, n_resamples=100, metric='sum')
        median = np.median(arr)
        return (median, sum, std, lb, ub)       
#######################################################################################

#######################################################################################
def bootstrap(data, n_resamples=10000, metric='mean'):
    boot_vals = []
    for b in range(n_resamples):
        boot_samp = np.random.choice(data, size=len(data))

        if metric == 'mean':
            boot_metric = np.mean(boot_samp)
        elif metric == 'sum':
            boot_metric = np.sum(boot_samp)

        #print(b, boot_metric)
        boot_vals += [boot_metric]

    lb = np.percentile(boot_vals, 2.5)
    ub = np.percentile(boot_vals, 97.5)
    return (lb, ub)
#######################################################################################

def main():
    """
    Main function to parse arguments and run the DataAnalysis class.
    """
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory")


    args = parser.parse_args()

    analysis = DataAnalysis(outpath=args.outpath)
    analysis.Population_anal()

if __name__ == "__main__":
    main()

