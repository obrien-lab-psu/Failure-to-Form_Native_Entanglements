import os
from scipy.stats import false_discovery_control
import glob
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

class Plotter:
    """
    slug_path = ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/

    Figure 1a: 
        Association between non-refoldability and the presence of a native entanglement
        C buffers at spa50 and LiPMScov 50
        EXP, AF, and both EXP & AF
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Association_Native_Entanglements_and_Misfolding/EntanglementsAndNonrefoldability/Plot_EntanglementsAndNonrefoldability/entanglement_and_nonrefoldability_plot_data_all_genes.csv

    Figure 1b:
        Association between change in proteolysis suseptibility and protein region (entangled versus non-entangled)
        C, Cpsm, CD, CG buffers at spa50 and LiPMScov 50 and has to have an entanglement
        EXP, AF, and both EXP & AF
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/AF/whole_proteome/ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/PSM/ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/AF/PSM/ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
        
    Figure 1c: 
        Association between change in proteolysis suseptibility and protein region (entangled versus non-entangled)
        C, CD, CG buffers at spa50 and LiPMScov 50 and has to have an entanglement and be essential
        C, CD, CG buffers at spa50 and LiPMScov 50 and has to have an entanglement and be non-essential   
        EXP, AF, and both EXP & AF
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/essential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/nonessential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv  
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/AF/whole_proteome/essential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/AF/whole_proteome/nonessential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv

    """
    def __init__(self, args):
        """
        Initializes the plotter with regression data.

        :param data: A DataFrame containing regression data.
        """
        self.slug_path = args.slug_path
        print(f'self.slug_path: {self.slug_path}')
        self.out_path = args.out_path
    #################################################################################################################

    #################################################################################################################
    def load_data(self,):
        print(f'Loading data from SLUG: {self.slug_path}')

        #######################################
        ## Load Figure 1a data
        #../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Association_Native_Entanglements_and_Misfolding/EntanglementsAndNonrefoldability/Plot_EntanglementsAndNonrefoldability/entanglement_and_nonrefoldability_plot_data_all_genes.cs
        inp = f'{self.slug_path}/Association_Native_Entanglements_and_Misfolding/EntanglementsAndNonrefoldability/Plot_EntanglementsAndNonrefoldability/entanglement_and_nonrefoldability_plot_data_all_genes.csv'
        print(f'inp: {inp}')
        self.Figure_1a_df = pd.read_csv(inp)
        print(f'Figure_1a_df:\n{self.Figure_1a_df}')

        #######################################
        ## Load Figure 1b data
        #../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
        #../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/AF/whole_proteome/ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
        #../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/PSM/ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
        #../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/AF/PSM/ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
        EXP_whole = f'{self.slug_path}/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv'
        print(f'EXP_whole: {EXP_whole}')
        EXP_whole_df = pd.read_csv(EXP_whole)
        EXP_whole_df['tag'] = 'EXP_Whole'
        print(f'EXP_whole_df:\n{EXP_whole_df}')

        AF_whole = f'{self.slug_path}/Modeling_Odds_of_Misfolding/Regressions/Plots/AF/whole_proteome/ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv'
        print(f'AF_whole: {AF_whole}')
        AF_whole_df = pd.read_csv(AF_whole)
        AF_whole_df['tag'] = 'AF_Whole'
        print(f'AF_whole_df:\n{AF_whole_df}')

        EXP_PSM = f'{self.slug_path}/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/PSM/ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv'
        print(f'EXP_PSM: {EXP_PSM}')
        EXP_PSM_df = pd.read_csv(EXP_PSM)
        EXP_PSM_df['tag'] = 'EXP_PSM'
        print(f'EXP_PSM_df:\n{EXP_PSM_df}')
        #AF_PSM = f'{self.slug_path}/Modeling_Odds_of_Misfolding/Regressions/Plots/AF/PSM/ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv'
        #print(f'AF_PSM: {AF_PSM}')
        #AF_PSM_df = pd.read_csv(AF_PSM)
        AF_PSM_df = EXP_PSM_df.copy()
        AF_PSM_df['tag'] = 'AF_PSM'
        print(f'AF_PSM_df:\n{AF_PSM_df}')

        self.Figure_1b_df = pd.concat((EXP_whole_df, AF_whole_df, EXP_PSM_df, AF_PSM_df))
        print(f'self.Figure_1b_df:\n{self.Figure_1b_df}')

        #######################################
        ## Load Figure 1c data
        #../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/essential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
        #../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/nonessential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv  
        #../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/AF/whole_proteome/essential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
        #./../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/AF/whole_proteome/nonessential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
        EXP_Ess = f'{self.slug_path}/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/essential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv'
        print(f'EXP_Ess: {EXP_Ess}')
        EXP_Ess_df = pd.read_csv(EXP_Ess)
        EXP_Ess_df['tag'] = 'EXP_Ess'
        print(f'EXP_Ess_df:\n{EXP_Ess_df}')
        EXP_NonEss = f'{self.slug_path}/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/nonessential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv'
        print(f'EXP_NonEss: {EXP_NonEss}')
        EXP_NonEss_df = pd.read_csv(EXP_NonEss)
        EXP_NonEss_df['tag'] = 'EXP_NonEss'
        print(f'EXP_NonEss_df:\n{EXP_NonEss_df}')

        AF_Ess = f'{self.slug_path}/Modeling_Odds_of_Misfolding/Regressions/Plots/AF/whole_proteome/essential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv'
        print(f'AF_Ess: {AF_Ess}')
        AF_Ess_df = pd.read_csv(AF_Ess)
        AF_Ess_df['tag'] = 'AF_Ess'
        print(f'AF_Ess_df:\n{AF_Ess_df}')
        AF_NonEss = f'{self.slug_path}/Modeling_Odds_of_Misfolding/Regressions/Plots/AF/whole_proteome/nonessential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv'
        print(f'AF_NonEss: {AF_NonEss}')
        AF_NonEss_df = pd.read_csv(AF_NonEss)
        AF_NonEss_df['tag'] = 'AF_NonEss'
        print(f'AF_NonEss_df:\n{AF_NonEss_df}')

        self.Figure_1c_df = pd.concat((EXP_Ess_df, EXP_NonEss_df, AF_Ess_df, AF_NonEss_df))
        print(f'self.Figure_1c_df:\n{self.Figure_1c_df}')
        
    #################################################################################################################

    #################################################################################################################
    def plot_Figure_1a(self,):
        """
        Figure 1a: 
            Association between non-refoldability and the presence of a native entanglement
            C buffers at spa50 and LiPMScov 50
            EXP, AF, and both EXP & AF
        """
        buff_tag = {'C': 'cyto-serum', 'CD': '+DnaK', 'CG': '+GroEL'}

        ## Plot Figure 1a
        Figure_1a_outfile = os.path.join(self.out_path, f'Figure_1a.png')
        Figure_1a_outfile_csv = os.path.join(self.out_path, f'Figure_1a.csv')
        print(f'Figure_1a_outfile: {Figure_1a_outfile}')
        print(f'Figure_1a_outfile_csv: {Figure_1a_outfile_csv}')

        Figure_1a_df = self.Figure_1a_df.copy()
        Figure_1a_df = Figure_1a_df[Figure_1a_df['spa'] == 50]
        Figure_1a_df['xlabel'] = [buff_tag[b] for b in Figure_1a_df['buff'].values]

        # Define the custom order for the 'buff' column
        custom_order = ['C', 'CD', 'CG']
        Figure_1a_df['buff'] = pd.Categorical(Figure_1a_df['buff'], categories=custom_order, ordered=True)
        Figure_1a_df = Figure_1a_df.sort_values('buff')

        # Define the custom order for the 'label' column
        custom_order = ['EXP', 'AF']
        Figure_1a_df['label'] = pd.Categorical(Figure_1a_df['label'], categories=custom_order, ordered=True)
        Figure_1a_df = Figure_1a_df.sort_values('label')
        print(Figure_1a_df)

        # Save the plot df
        Figure_1a_df.to_csv(Figure_1a_outfile_csv)
        print(f'SAVED: {Figure_1a_outfile_csv}')

        # Set up the plot
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(7, 3), sharey=True)
        
        # Iterate through the groups and plot each trace with error bars
        maxy = 0
        for i, (label, group) in enumerate(Figure_1a_df.groupby('label')):
            
            # Calculate error bars
            yerr_lower = group['OR'] - group['OR_lb']
            yerr_upper = group['OR_ub'] - group['OR']
            yerr = [yerr_lower, yerr_upper]
            
            # Plot the trace with error bars
            axes[i].errorbar(
                group['xlabel'], 
                group['OR'], 
                yerr=yerr, 
                fmt='o',  # Line and point markers
                label=f"{label}",
                capsize=4  # Add caps to the error bars
            )

            ## mark pvalue annotations
            # Annotate significant p-values with '*'
            for x, y, err, p in zip(group['xlabel'], group['OR'], yerr_upper, group['pvalue']):
                if p > 0.05:
                    p_annot = ""
                elif p < 0.05 and p >= 0.01:
                    p_annot = "*"
                elif p < 0.01 and p >= 0.001:
                    p_annot = "**"
                elif p < 0.001:
                    p_annot = "***" 

                axes[i].text(
                    x, 
                    y + err + 0.1,  # Place the annotation slightly above the upper error bar
                    p_annot, 
                    fontsize=12, 
                    ha='center', 
                    color='red')                  
                print(p_annot, y + err + 0.1)

                if y + err + 0.1 > maxy:
                    maxy = y + err + 0.1
                
            # Add a dashed black line at y=1    
            axes[i].axhline(y=1, color='black', linestyle='--', linewidth=1)
            axes[i].set_xlabel("Buffer")
            axes[i].set_ylabel("Odds Ratio (OR)")
            axes[i].set_xticklabels(group['xlabel'], rotation=45, ha='right')  # Rotate x-axis labels for better readability
            #axes[i].legend()
            axes[i].set_title(f"Association of Non-refoldability\nand Native Entanglements ({label})", fontsize=10)
            axes[i].grid(alpha=0.3)
            
        # Customize the plot
        print(f'Max(y) ceil: {maxy} {np.ceil(maxy)}')
        plt.ylim(0,np.ceil(maxy) + 1)
        
        # Show the plot
        plt.tight_layout()
        plt.savefig(Figure_1a_outfile)
        print(f'SAVED: Figure_1a_outfile: {Figure_1a_outfile}')
    #################################################################################################################

    #################################################################################################################
    def plot_Figure_1b(self,):
        """
        Figure 1b:
            Association between change in proteolysis suseptibility and protein region (entangled versus non-entangled)
            C, Cpsm, CD, CG buffers at spa50 and LiPMScov 50 and has to have an entanglement
            EXP, AF, and both EXP & AF
        """
        buff_tag = {'C': 'cyto-serum', 'CD': '+DnaK', 'CG': '+GroEL'}

        ## Plot Figure 1a
        Figure_1b_outfile = os.path.join(self.out_path, f'Figure_1b.png')
        Figure_1b_outfile_csv = os.path.join(self.out_path, f'Figure_1b.csv')
        print(f'Figure_1b_outfile: {Figure_1b_outfile}')
        print(f'Figure_1b_outfile_csv: {Figure_1b_outfile_csv}')
        
        Figure_1b_df = self.Figure_1b_df.copy()
        print(Figure_1b_df)
        Figure_1b_df = Figure_1b_df[Figure_1b_df['spa'] == 50]
        xlabels = []
        for ti, t in enumerate(Figure_1b_df['tag'].values):
            b = Figure_1b_df['buff'].values[ti]
            b = buff_tag[b]
            if 'PSM' in t:
                xlabels += [f"{b}(PSM)"]
            else:
                xlabels += [f"{b}"]
        Figure_1b_df['xlabel'] = xlabels
        Figure_1b_df['label'] = [t.split('_')[0] for t in Figure_1b_df['tag'].values]

        # Define the custom order for the 'buff' column
        custom_order = ['C', 'CD', 'CG']
        Figure_1b_df['buff'] = pd.Categorical(Figure_1b_df['buff'], categories=custom_order, ordered=True)
        Figure_1b_df = Figure_1b_df.sort_values('buff')

        # Define the custom order for the 'label' column
        custom_order = ['EXP', 'AF']
        Figure_1b_df['label'] = pd.Categorical(Figure_1b_df['label'], categories=custom_order, ordered=True)
        Figure_1b_df = Figure_1b_df.sort_values('label')
        print(Figure_1b_df)

        # Save the plot df
        Figure_1b_df.to_csv(Figure_1b_outfile_csv)
        print(f'SAVED: {Figure_1b_outfile_csv}')

        # Set up the plot
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(7, 3), sharey=True)
        
        # Iterate through the groups and plot each trace with error bars
        maxy = 0
        for i, (label, group) in enumerate(Figure_1b_df.groupby('label')):
            
            # Calculate error bars
            yerr_lower = group['OR'] - group['OR_lb']
            yerr_upper = group['OR_ub'] - group['OR']
            yerr = [yerr_lower, yerr_upper]
            
            # Plot the trace with error bars
            axes[i].errorbar(
                group['xlabel'], 
                group['OR'], 
                yerr=yerr, 
                fmt='o',  # Line and point markers
                label=f"{label}",
                capsize=4  # Add caps to the error bars
            )

            ## mark pvalue annotations
            # Annotate significant p-values with '*'
            for x, y, err, p in zip(group['xlabel'], group['OR'], yerr_upper, group['pvalues']):
                if p > 0.05:
                    p_annot = ""
                elif p < 0.05 and p >= 0.01:
                    p_annot = "*"
                elif p < 0.01 and p >= 0.001:
                    p_annot = "**"
                elif p < 0.001:
                    p_annot = "***" 

                axes[i].text(
                    x, 
                    y + err + 0.1,  # Place the annotation slightly above the upper error bar
                    p_annot, 
                    fontsize=12, 
                    ha='center', 
                    color='red')                  
                print(p_annot, y + err + 0.1)

                if y + err + 0.1 > maxy:
                    maxy = y + err + 0.1
                
            # Add a dashed black line at y=1    
            axes[i].axhline(y=1, color='black', linestyle='--', linewidth=1)
            axes[i].set_xlabel("Buffer")
            axes[i].set_ylabel("Odds Ratio (OR)")
            axes[i].set_xticklabels(group['xlabel'], rotation=45, ha='right')  # Rotate x-axis labels for better readability
            #axes[i].legend()
            axes[i].set_title(f"Association of change in prot. suseptibility\nand Entanglement region ({label})", fontsize=10)
            axes[i].grid(alpha=0.3)
            
        # Customize the plot
        print(f'Max(y) ceil: {maxy} {np.ceil(maxy)}')
        plt.ylim(0,np.ceil(maxy) + 0.5)
        
        # Show the plot
        plt.tight_layout()
        plt.savefig(Figure_1b_outfile)
        print(f'SAVED: Figure_1b_outfile: {Figure_1b_outfile}')
    #################################################################################################################

    #################################################################################################################
    def plot_Figure_1c(self,):
        """
        Figure 1c: 
            Association between change in proteolysis suseptibility and protein region (entangled versus non-entangled)
            C, CD, CG buffers at spa50 and LiPMScov 50 and has to have an entanglement and be essential
            C, CD, CG buffers at spa50 and LiPMScov 50 and has to have an entanglement and be non-essential   
            EXP, AF, and both EXP & AF
        """
        buff_tag = {'C': 'cyto-serum', 'CD': '+DnaK', 'CG': '+GroEL'}

        ## Plot Figure 1c
        Figure_1c_outfile = os.path.join(self.out_path, f'Figure_1c.png')
        Figure_1c_outfile_csv = os.path.join(self.out_path, f'Figure_1c.csv')
        print(f'Figure_1c_outfile: {Figure_1c_outfile}')
        print(f'Figure_1c_outfile_csv: {Figure_1c_outfile_csv}')
        
        Figure_1c_df = self.Figure_1c_df.copy()
        Figure_1c_df = Figure_1c_df[Figure_1c_df['spa'] == 50]
        print(Figure_1c_df)

        xlabels = []
        for ti, t in enumerate(Figure_1c_df['tag'].values):
            b = Figure_1c_df['buff'].values[ti]
            b = buff_tag[b]
            xlabels += [f"{b}"]
        Figure_1c_df['xlabel'] = xlabels
        Figure_1c_df['label'] = [t.split('_')[0] for t in Figure_1c_df['tag'].values]
        tags = []
        for t in Figure_1c_df['tag'].values:
            if 'NonEss' in t:
                tags += ['NonEssential']
            else:
                tags += ['Essential']
        Figure_1c_df['tag'] = tags

        # Define the custom order for the 'label' column
        custom_order = ['EXP', 'AF']
        Figure_1c_df['label'] = pd.Categorical(Figure_1c_df['label'], categories=custom_order, ordered=True)
        Figure_1c_df = Figure_1c_df.sort_values('label')
        print(Figure_1c_df)

        # Save the plot df
        Figure_1c_df.to_csv(Figure_1c_outfile_csv)
        print(f'SAVED: {Figure_1c_outfile_csv}')

        # Set up the plot
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(7, 7), sharey=True)
        
        # Iterate through the groups and plot each trace with error bars
        maxy = 0
        for i, (label, group) in enumerate(Figure_1c_df.groupby('label')):
            for j, (tag, tag_group) in enumerate(group.groupby('tag')):
                print(tag_group)
            
                # Calculate error bars
                yerr_lower = tag_group['OR'] - tag_group['OR_lb']
                yerr_upper = tag_group['OR_ub'] - tag_group['OR']
                yerr = [yerr_lower, yerr_upper]
                
                # Plot the trace with error bars
                axes[i, j].errorbar(
                    tag_group['xlabel'], 
                    tag_group['OR'], 
                    yerr=yerr, 
                    fmt='o',  # Line and point markers
                    label=f"{tag}",
                    capsize=4  # Add caps to the error bars
                )

                ## mark pvalue annotations
                # Annotate significant p-values with '*'
                for x, y, err, p in zip(tag_group['xlabel'], tag_group['OR'], yerr_upper, tag_group['pvalues']):
                    if p > 0.05:
                        p_annot = ""
                    elif p < 0.05 and p >= 0.01:
                        p_annot = "*"
                    elif p < 0.01 and p >= 0.001:
                        p_annot = "**"
                    elif p < 0.001:
                        p_annot = "***" 

                    axes[i, j].text(
                        x, 
                        y + err + 0.1,  # Place the annotation slightly above the upper error bar
                        p_annot, 
                        fontsize=12, 
                        ha='center', 
                        color='red')                  
                    print(p_annot, y + err + 0.1)

                    if y + err + 0.1 > maxy:
                        maxy = y + err + 0.1
                    
                # Add a dashed black line at y=1    
                axes[i, j].axhline(y=1, color='black', linestyle='--', linewidth=1)
                axes[i, j].set_xlabel("Buffer")
                axes[i, j].set_ylabel("Odds Ratio (OR)")
                axes[i, j].set_xticklabels(tag_group['xlabel'], rotation=45, ha='right')  # Rotate x-axis labels for better readability
                axes[i, j].set_title(f"Association of change in prot. suseptibility\nand Entanglement region ({label}, {tag})", fontsize=10)
                axes[i, j].grid(alpha=0.3)
                #axes[i, j].legend()
            
        # Customize the plot
        print(f'Max(y) ceil: {maxy} {np.ceil(maxy)}')
        plt.ylim(0,np.ceil(maxy))
        
        # Show the plot
        plt.tight_layout()
        plt.savefig(Figure_1c_outfile)
        print(f'SAVED: Figure_1c_outfile: {Figure_1c_outfile}')
        #################################################################################################################


def main():
    """
    slug_path = ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/

    Figure 2a: 
        Entanglement complexity as a discriminator
        C, CD, CG buffers at spa50 and LiPMScov 50
        EXP, AF
    
    Figure 2b:
        DnaK binding motif its
        C, CD, CG buffers at spa50 and LiPMScov 50 and has to have an entanglement
        EXP, AF
        Dall, Dthread, Dloop
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Chaperone_Client_Associations/Dnak_simple_motif_scan/Plots/consolidated_Dnak_binding_data_EXP.csv

    """
    parser = argparse.ArgumentParser(description="Process regression data and generate plots.")
    parser.add_argument("-s", "--slug_path", type=str, required=True, help="Path to the slug containing all the raw data for this paper.")
    parser.add_argument("-o", "--out_path", type=str, required=True, help="Path to output directory.")
    args = parser.parse_args()
    print(args)

    if not os.path.exists(args.out_path):
        os.makedirs(args.out_path)
        print(f'MADE: {args.out_path}')

    plotter = Plotter(args)
    plotter.load_data()
    plotter.plot_Figure_1a()
    plotter.plot_Figure_1b()
    plotter.plot_Figure_1c()

    print('NORMAL TERMINATION')

if __name__ == "__main__":
    main()

