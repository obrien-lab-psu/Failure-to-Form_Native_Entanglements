import os
from scipy.stats import false_discovery_control
import glob
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib import font_manager as fm
import matplotlib as mpl

class Plotter:
    """
    -------------------------------------------------------------------------------------------------
    Figure formating requirements for Nature
    https://research-figure-guide.nature.com/figures/preparing-figures-our-specifications/?utm_source=chatgpt.com#we-require

    Figure Sizing and Positioning
    Dimensions: Figures should fit within a single page, ideally leaving space for the legend below. 
                The maximum page dimensions are 180 mm wide by 170 mm tall (170 mm to accommodate a legend underneath).
    
    Placement: Position each figure centrally on a new page. Avoid placing multiple figures on the same page. 

    Line Weights: Set lines and strokes between 0.25 and 1 point to ensure clarity.

    Text and Fonts: All text within figures should be legible and editable. 
                    Use standard sans-serif fonts like Arial or Helvetica. 
                    Avoid outlining text and ensure fonts are embedded (True Type 2 or 42). 
                    Text size should range between 5-point (minimum) and 7-point (maximum)
                    Present amino-acid sequences in Courier (or other monospaced) font using the one-letter code in lines of 50 or 100 characters
                    Separate panels in multi-panelled figures should be labelled with 8-pt bold, upright (not italic) and lowercase a, b, c, etc.
                    If you are using Python please use the following setting: Matplotlib.rcParams['pdf.fonttype']=42

    Color and Accessibility: Use accessible color palettes to accommodate readers with color vision deficiencies. 
                             Avoid red/green combinations and rainbow scales. Ensure high-contrast text (>4.5 contrast ratio) for readability.
                             Use the RGB color space never the CMYK color space  

    Panel Arrangement: Arrange multi-panel figures neatly, minimizing white space and ordering panels alphabetically. 
                       Separate panels in multi-panelled figures should be labelled with 8-pt bold, upright (not italic) and lowercase a, b, c, etc.

    Axis labels and tickmarks: All tick marks should be included for any number on the axis. 
                               All axis should have a label with units in parentheses

    Avoid: Background grid lines
           superfluous icons and other decorative elements
           drop shadows
           text place on top of a busy image and hard-to-read background
           overlaping text
           coloured text
           pattern filling of bars, pies, ect...

    legends: should use color boxes not colored text
        
    Exporting: please export figure panels as vector artwork 
                .pdf or .eps preferred
                For images, minimum 450 dpi
    -------------------------------------------------------------------------------------------------
    Raw data paths required and summary of each panel

    slug_path = ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/

    Figure is 1 row by 3 column 
    Figure 4a (row 1 column 1): 
        Association between change in proteolysis suseptibility and protein region (entangled versus non-entangled)
        Essential versus non-essentail for Cyto-serum at spa50 and LiPMScov 50 in EXP dataset
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/essential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/nonessential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv  
                
    Figure 4b (row 1 column 2):
        Association between change in proteolysis suseptibility and protein region (entangled versus non-entangled)
        Essential versus non-essentail for +DnaK at spa50 and LiPMScov 50 in EXP dataset
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/essential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/nonessential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv  

    Figure 4c: 
        Association between change in proteolysis suseptibility and protein region (entangled versus non-entangled)
        Essential versus non-essentail for +GroEL at spa50 and LiPMScov 50 in EXP dataset
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/essential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/nonessential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv  
                
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
    def plot_Figure_4a(self, ax):
        """
        Figure 4a (row 1 column 1): 
            Association between change in proteolysis suseptibility and protein region (entangled versus non-entangled)
            Essential versus non-essentail for Cyto-serum at spa50 and LiPMScov 50 in EXP dataset
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/essential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/nonessential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv       
        """
        #######################################
        ## Load Figure 1a data
        inp = f'{self.slug_path}/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/essential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv'
        print(f'inp: {inp}')
        self.Figure_4a_ess_df = pd.read_csv(inp)
        self.Figure_4a_ess_df = self.Figure_4a_ess_df[self.Figure_4a_ess_df['spa'] == 50]
        self.Figure_4a_ess_df = self.Figure_4a_ess_df[self.Figure_4a_ess_df['cov'] == 50]
        self.Figure_4a_ess_df = self.Figure_4a_ess_df[self.Figure_4a_ess_df['buff'] == 'C']
        self.Figure_4a_ess_df['xlabel'] = 'Essential'
        self.Figure_4a_ess_df['color'] = 'blue'
        print(f'Figure_4a_ess_df:\n{self.Figure_4a_ess_df}')

        inp = f'{self.slug_path}/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/nonessential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv'
        print(f'inp: {inp}')
        self.Figure_4a_noness_df = pd.read_csv(inp)
        self.Figure_4a_noness_df = self.Figure_4a_noness_df[self.Figure_4a_noness_df['spa'] == 50]
        self.Figure_4a_noness_df = self.Figure_4a_noness_df[self.Figure_4a_noness_df['cov'] == 50]
        self.Figure_4a_noness_df = self.Figure_4a_noness_df[self.Figure_4a_noness_df['buff'] == 'C']
        self.Figure_4a_noness_df['xlabel'] = 'Non-Essential\nproteins'
        self.Figure_4a_noness_df['color'] = 'red'
        print(f'Figure_4a_noness_df:\n{self.Figure_4a_noness_df}')
        #######################################
      
        # Save the figure 1a raw plot df
        Figure_4a_df = pd.concat([self.Figure_4a_ess_df, self.Figure_4a_noness_df])
        print(Figure_4a_df)
        Figure_4a_outfile_csv = os.path.join(self.out_path, f'Figure_4a.csv')
        Figure_4a_df.to_csv(Figure_4a_outfile_csv)
        print(f'SAVED: {Figure_4a_outfile_csv}')
                
        # Adjust the linewidth of the axis spines
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)  # Set the linewidth for all spines

        # Adjust the linewidth of the ticks
        ax.tick_params(width=0.5)  # Both major and minor ticks

        maxy = 0
        # Calculate error bars
        yerr_lower = Figure_4a_df['OR'] - Figure_4a_df['OR_lb']
        yerr_upper = Figure_4a_df['OR_ub'] - Figure_4a_df['OR']
        yerr = [yerr_lower, yerr_upper]
        
        # Loop through each row of the dataframe and plot each point
        for index, row in Figure_4a_df.iterrows():
            ax.errorbar(
                row['xlabel'], 
                row['OR'], 
                yerr=[[row['OR'] - row['OR_lb']], [row['OR_ub'] - row['OR']]],  # Use the individual yerr for each point
                fmt='o', 
                capsize=3,
                color=row['color'],  # Use the color from the dataframe
                markersize=3,
                markerfacecolor=row['color'],
                elinewidth=0.5,
                markeredgewidth=0.5)
    
        ## mark pvalue annotations
        # Annotate significant p-values with '*'
        for x, y, err, p in zip(Figure_4a_df['xlabel'], Figure_4a_df['OR'], yerr_upper, Figure_4a_df['pvalues']):
            if p > 0.05:
                p_annot = ""
            elif p < 0.05 and p >= 0.01:
                p_annot = "*"
            elif p < 0.01 and p >= 0.001:
                p_annot = "**"
            elif p < 0.001:
                p_annot = "***" 

            ax.text(
                x, 
                y + err + 0.1,  # Place the annotation slightly above the upper error bar
                p_annot, 
                fontsize=6, 
                ha='center', 
                color='black')                  
            #print(p_annot, y + err + 0.1)

            if y + err + 0.1 > maxy:
                maxy = y + err + 0.1
            
        # add buffer label to subplot
        ax.text(0.95, 0.98, "cyto-serum", transform=ax.transAxes, va='top', ha='right', fontsize=6) 

        # Add a dashed black line at y=1    
        ax.axhline(y=1, color='#7E8D85', linestyle='--', linewidth=0.5)
        ax.set_ylabel("Odds ratio between\nentangled protein regions and misfolding")
        ax.set_xticklabels(Figure_4a_df['xlabel'], rotation=45, ha='right')  # Rotate x-axis labels for better readability
            
        # Customize the plot
        print(f'Max(y) ceil: {maxy} {np.ceil(maxy)}')
        ax.set_ylim(0.5,np.ceil(maxy) + 0.1)
        ax.set_xlim(-0.5, 1.5)

        ax.tick_params(axis='y', labelsize=6)
        ax.tick_params(axis='x', labelsize=6)
        # Remove the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        return ax   
    #################################################################################################################

    #################################################################################################################
    def plot_Figure_4b(self, ax):
        """
        Figure 4b (row 1 column 1): 
            Association between change in proteolysis suseptibility and protein region (entangled versus non-entangled)
            Essential versus non-essentail for +DnaK at spa50 and LiPMScov 50 in EXP dataset
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/essential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/nonessential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv       
        """
        #######################################
        ## Load Figure 1a data
        inp = f'{self.slug_path}/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/essential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv'
        print(f'inp: {inp}')
        self.Figure_4b_ess_df = pd.read_csv(inp)
        self.Figure_4b_ess_df = self.Figure_4b_ess_df[self.Figure_4b_ess_df['spa'] == 50]
        self.Figure_4b_ess_df = self.Figure_4b_ess_df[self.Figure_4b_ess_df['cov'] == 50]
        self.Figure_4b_ess_df = self.Figure_4b_ess_df[self.Figure_4b_ess_df['buff'] == 'CD']
        self.Figure_4b_ess_df['xlabel'] = 'Essential'
        self.Figure_4b_ess_df['color'] = 'blue'
        print(f'Figure_4b_ess_df:\n{self.Figure_4b_ess_df}')

        inp = f'{self.slug_path}/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/nonessential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv'
        print(f'inp: {inp}')
        self.Figure_4b_noness_df = pd.read_csv(inp)
        self.Figure_4b_noness_df = self.Figure_4b_noness_df[self.Figure_4b_noness_df['spa'] == 50]
        self.Figure_4b_noness_df = self.Figure_4b_noness_df[self.Figure_4b_noness_df['cov'] == 50]
        self.Figure_4b_noness_df = self.Figure_4b_noness_df[self.Figure_4b_noness_df['buff'] == 'CD']
        self.Figure_4b_noness_df['xlabel'] = 'Non-Essential\nproteins'
        self.Figure_4b_noness_df['color'] = 'red'
        print(f'Figure_4b_noness_df:\n{self.Figure_4b_noness_df}')
        #######################################
      
        # Save the figure 1a raw plot df
        Figure_4b_df = pd.concat([self.Figure_4b_ess_df, self.Figure_4b_noness_df])
        print(Figure_4b_df)
        Figure_4b_outfile_csv = os.path.join(self.out_path, f'Figure_4b.csv')
        Figure_4b_df.to_csv(Figure_4b_outfile_csv)
        print(f'SAVED: {Figure_4b_outfile_csv}')
                
        # Adjust the linewidth of the axis spines
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)  # Set the linewidth for all spines

        # Adjust the linewidth of the ticks
        ax.tick_params(width=0.5)  # Both major and minor ticks

        maxy = 0
        # Calculate error bars
        yerr_lower = Figure_4b_df['OR'] - Figure_4b_df['OR_lb']
        yerr_upper = Figure_4b_df['OR_ub'] - Figure_4b_df['OR']
        yerr = [yerr_lower, yerr_upper]
        
        # Loop through each row of the dataframe and plot each point
        for index, row in Figure_4b_df.iterrows():
            ax.errorbar(
                row['xlabel'], 
                row['OR'], 
                yerr=[[row['OR'] - row['OR_lb']], [row['OR_ub'] - row['OR']]],  # Use the individual yerr for each point
                fmt='o', 
                capsize=3,
                color=row['color'],  # Use the color from the dataframe
                markersize=3,
                markerfacecolor=row['color'],
                elinewidth=0.5,
                markeredgewidth=0.5)

        ## mark pvalue annotations
        # Annotate significant p-values with '*'
        for x, y, err, p in zip(Figure_4b_df['xlabel'], Figure_4b_df['OR'], yerr_upper, Figure_4b_df['pvalues']):
            if p > 0.05:
                p_annot = ""
            elif p < 0.05 and p >= 0.01:
                p_annot = "*"
            elif p < 0.01 and p >= 0.001:
                p_annot = "**"
            elif p < 0.001:
                p_annot = "***" 

            ax.text(
                x, 
                y + err + 0.1,  # Place the annotation slightly above the upper error bar
                p_annot, 
                fontsize=6, 
                ha='center', 
                color='black')                  
            #print(p_annot, y + err + 0.1)

            if y + err + 0.1 > maxy:
                maxy = y + err + 0.1
            
        # add buffer label to subplot
        ax.text(0.95, 0.98, "+DnaK", transform=ax.transAxes, va='top', ha='right', fontsize=6) 

        # Add a dashed black line at y=1    
        ax.axhline(y=1, color='#7E8D85', linestyle='--', linewidth=0.5)
        #ax.set_ylabel("Odds ratio for association of\nmisfolding and entanglement region")
        ax.set_xticklabels(Figure_4b_df['xlabel'], rotation=45, ha='right')  # Rotate x-axis labels for better readability
            
        # Customize the plot
        print(f'Max(y) ceil: {maxy} {np.ceil(maxy)}')
        ax.set_ylim(0.5,np.ceil(maxy) + 0.1)
        ax.set_xlim(-0.5, 1.5)

        ax.tick_params(axis='y', labelsize=6)
        ax.tick_params(axis='x', labelsize=6)
        # Remove the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_yticklabels([])
        return ax   
    #################################################################################################################

    #################################################################################################################
    def plot_Figure_4c(self, ax):
        """
        Figure 4c (row 1 column 1): 
            Association between change in proteolysis suseptibility and protein region (entangled versus non-entangled)
            Essential versus non-essentail for +GroEL at spa50 and LiPMScov 50 in EXP dataset
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/essential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/nonessential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv       
        """
        #######################################
        ## Load Figure 4c data
        inp = f'{self.slug_path}/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/essential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv'
        print(f'inp: {inp}')
        self.Figure_4c_ess_df = pd.read_csv(inp)
        self.Figure_4c_ess_df = self.Figure_4c_ess_df[self.Figure_4c_ess_df['spa'] == 50]
        self.Figure_4c_ess_df = self.Figure_4c_ess_df[self.Figure_4c_ess_df['cov'] == 50]
        self.Figure_4c_ess_df = self.Figure_4c_ess_df[self.Figure_4c_ess_df['buff'] == 'CG']
        self.Figure_4c_ess_df['xlabel'] = 'Essential'
        self.Figure_4c_ess_df['color'] = 'blue'
        print(f'Figure_4c_ess_df:\n{self.Figure_4c_ess_df}')

        inp = f'{self.slug_path}/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/nonessential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv'
        print(f'inp: {inp}')
        self.Figure_4c_noness_df = pd.read_csv(inp)
        self.Figure_4c_noness_df = self.Figure_4c_noness_df[self.Figure_4c_noness_df['spa'] == 50]
        self.Figure_4c_noness_df = self.Figure_4c_noness_df[self.Figure_4c_noness_df['cov'] == 50]
        self.Figure_4c_noness_df = self.Figure_4c_noness_df[self.Figure_4c_noness_df['buff'] == 'CG']
        self.Figure_4c_noness_df['xlabel'] = 'Non-Essential\nproteins'
        self.Figure_4c_noness_df['color'] = 'red'
        print(f'Figure_4c_noness_df:\n{self.Figure_4c_noness_df}')
        #######################################
      
        # Save the figure 1a raw plot df
        Figure_4c_df = pd.concat([self.Figure_4c_ess_df, self.Figure_4c_noness_df])
        print(Figure_4c_df)
        Figure_4c_outfile_csv = os.path.join(self.out_path, f'Figure_4c.csv')
        Figure_4c_df.to_csv(Figure_4c_outfile_csv)
        print(f'SAVED: {Figure_4c_outfile_csv}')
                
        # Adjust the linewidth of the axis spines
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)  # Set the linewidth for all spines

        # Adjust the linewidth of the ticks
        ax.tick_params(width=0.5)  # Both major and minor ticks

        maxy = 0
        # Calculate error bars
        yerr_lower = Figure_4c_df['OR'] - Figure_4c_df['OR_lb']
        yerr_upper = Figure_4c_df['OR_ub'] - Figure_4c_df['OR']
        yerr = [yerr_lower, yerr_upper]
        
        # Loop through each row of the dataframe and plot each point
        for index, row in Figure_4c_df.iterrows():
            ax.errorbar(
                row['xlabel'], 
                row['OR'], 
                yerr=[[row['OR'] - row['OR_lb']], [row['OR_ub'] - row['OR']]],  # Use the individual yerr for each point
                fmt='o', 
                capsize=3,
                color=row['color'],  # Use the color from the dataframe
                markersize=3,
                markerfacecolor=row['color'],
                elinewidth=0.5,
                markeredgewidth=0.5)

        ## mark pvalue annotations
        # Annotate significant p-values with '*'
        for x, y, err, p in zip(Figure_4c_df['xlabel'], Figure_4c_df['OR'], yerr_upper, Figure_4c_df['pvalues']):
            if p > 0.05:
                p_annot = ""
            elif p < 0.05 and p >= 0.01:
                p_annot = "*"
            elif p < 0.01 and p >= 0.001:
                p_annot = "**"
            elif p < 0.001:
                p_annot = "***" 

            ax.text(
                x, 
                y + err + 0.1,  # Place the annotation slightly above the upper error bar
                p_annot, 
                fontsize=6, 
                ha='center', 
                color='black')                  
            #print(p_annot, y + err + 0.1)

            if y + err + 0.1 > maxy:
                maxy = y + err + 0.1
            
        # add buffer label to subplot
        ax.text(0.95, 0.98, "+GroEL", transform=ax.transAxes, va='top', ha='right', fontsize=6) 

        # Add a dashed black line at y=1    
        ax.axhline(y=1, color='#7E8D85', linestyle='--', linewidth=0.5)
        #ax.set_ylabel("Odds ratio for association of\nmisfolding and entanglement region")
        ax.set_xticklabels(Figure_4c_df['xlabel'], rotation=45, ha='right')  # Rotate x-axis labels for better readability
            
        # Customize the plot
        print(f'Max(y) ceil: {maxy} {np.ceil(maxy)}')
        ax.set_ylim(0.5,np.ceil(maxy) + 0.1)
        ax.set_xlim(-0.5, 1.5)

        ax.tick_params(axis='y', labelsize=6)
        ax.tick_params(axis='x', labelsize=6)
        # Remove the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_yticklabels([])
        return ax   
    #################################################################################################################

    #################################################################################################################
    def plot_arrows(self, ax):
        # Arrow pointing up (Positive Association)
        ax.annotate(
            '',  # No text for the arrow itself
            xy=(0.5, 1),  # End point of the arrow (top center)
            xytext=(0.5, 0.32),  # Start point of the arrow (middle)
            color='#7E8D85',
            arrowprops=dict(facecolor='#7E8D85', edgecolor='#7E8D85', arrowstyle='->', lw=0.5)
        )
        # Label for the upward arrow
        ax.text(0.8, 0.65, 'Positive\nAssociation', rotation=90, va='center', ha='center', color='#7E8D85')

        # Arrow pointing down (Negative Association)
        ax.annotate(
            '',  # No text for the arrow itself
            xy=(0.5, 0),  # End point of the arrow (bottom center)
            xytext=(0.5, 0.32),  # Start point of the arrow (middle)
            color='#7E8D85',
            arrowprops=dict(facecolor='#7E8D85', edgecolor='#7E8D85', arrowstyle='->', lw=0.5)
        )
        # Label for the downward arrow
        ax.text(0.8, 0.175, 'Negative\nAssociation', rotation=90, va='center', ha='center', color='#7E8D85')

        # Hide the axes
        #ax.axis('tight')
        ax.axis('off')
    #################################################################################################################

def mm_to_inches(mm):
    return mm / 25.4

# Function to compute the CDF
def compute_cdf(data):
    sorted_data = np.sort(data)
    cdf = np.arange(1, len(sorted_data) + 1) / len(sorted_data)
    return sorted_data, cdf

def main():
    """
    Creates Figure 4
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

    # Create a master figure and axes
    fig_width_mm = 89  # Width in mm
    fig_height_mm = 70  # Height in mm
    #custom_font_path = "/storage/group/epo2/default/ims86/miniconda3/envs/FtoF/fonts/Arial.ttf" # Path to your custom font
    #arial_font = fm.FontProperties(fname=custom_font_path) # Create a FontProperties object
    plt.rcParams['font.family'] = 'Arial'  # Change to your desired font, e.g., 'Times New Roman', 'DejaVu Sans', etc.
    plt.rcParams['font.size'] = 6  # Default font size
    plt.rcParams['pdf.fonttype'] = 42
    fig, axs = plt.subplots(1, 4, 
                            figsize=(mm_to_inches(fig_width_mm), mm_to_inches(fig_height_mm)), 
                            dpi=600, constrained_layout=True, 
                            gridspec_kw={'width_ratios': [1, 1, 1, 0.5]})  

    ## make subplot figures 
    plotter.plot_Figure_4a(axs[0])
    plotter.plot_Figure_4b(axs[1])
    plotter.plot_Figure_4c(axs[2])
    plotter.plot_arrows(axs[3])


    #########################################################
    # Get the current position of the subplots and adjust their positions manually
    axs0_position = axs[0].get_position()
    width0, height0 = axs0_position.extents[2] - axs0_position.extents[0], axs0_position.extents[3] - axs0_position.extents[1]
    print(axs0_position, width0, height0)
    axs[0].set_position([0.13, 0.2, 0.15, 0.7])  # [left, bottom, width, height]

    bbox_in_fig_coords = axs[0].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bbox_in_fig_coords.x0, 0.975, 'a', fontsize=8, fontweight='bold', va='top', ha='left')


    axs1_position = axs[1].get_position()
    width1, height1 = axs1_position.extents[2] - axs1_position.extents[0], axs1_position.extents[3] - axs1_position.extents[1]
    print(axs1_position, width1, height1)
    axs[1].set_position([0.4375, 0.2, 0.15, 0.7])

    bbox_in_fig_coords = axs[1].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bbox_in_fig_coords.x0, 0.975, 'b', fontsize=8, fontweight='bold', va='top', ha='left')


    axs2_position = axs[2].get_position()
    width2, height2 = axs2_position.extents[2] - axs2_position.extents[0], axs2_position.extents[3] - axs2_position.extents[1]
    print(axs2_position, width2, height2)
    axs[2].set_position([0.75, 0.2, 0.15, 0.7])

    bbox_in_fig_coords = axs[2].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bbox_in_fig_coords.x0, 0.975, 'c', fontsize=8, fontweight='bold', va='top', ha='left')


    axs3_position = axs[3].get_position()
    width3, height3 = axs3_position.extents[2] - axs3_position.extents[0], axs3_position.extents[3] - axs3_position.extents[1]
    print(axs3_position, width3, height3)
    axs[3].set_position([0.9, 0.2, width3, 0.7])
    #########################################################

    # final formating and output
    # Automatically adjust layout
    #fig.tight_layout(pad=2.0)  # 'pad' is the overall padding
    figure_outpath = os.path.join(args.out_path, 'Figure4.pdf')
    plt.savefig(figure_outpath)
    print(f'SAVED: {figure_outpath}')

    figure_outpath = os.path.join(args.out_path, 'Figure4.svg')
    plt.savefig(figure_outpath)
    print(f'SAVED: {figure_outpath}')

    figure_outpath = os.path.join(args.out_path, 'Figure4.png')
    plt.savefig(figure_outpath)
    print(f'SAVED: {figure_outpath}')

    print('NORMAL TERMINATION')

if __name__ == "__main__":
    main()

