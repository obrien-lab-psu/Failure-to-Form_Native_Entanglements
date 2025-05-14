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
    Figure 3a (row 1 column 1): 
        Association between change in proteolysis suseptibility and protein region (entangled versus non-entangled)
        Essential versus non-essentail for Cyto-serum at spa50 and LiPMScov 50 in EXP dataset
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/essential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/nonessential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv  
                
    Figure 3b (row 1 column 2):
        Association between change in proteolysis suseptibility and protein region (entangled versus non-entangled)
        Essential versus non-essentail for +DnaK at spa50 and LiPMScov 50 in EXP dataset
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/essential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/nonessential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv  

    Figure 3c: 
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
    def plot_MainTextFigure3_Supp(self, essential=True):
        """
        Figure 3a (row 1 column 1): 
            Association between change in proteolysis suseptibility and protein region (entangled versus non-entangled)
            Essential versus non-essentail for Cyto-serum at spa50 and LiPMScov 50 in EXP dataset
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/essential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/nonessential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv       
        """

        #######################################
        if essential == True:
            inp = f'{self.slug_path}/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/essential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv'
        else:
            inp = f'{self.slug_path}/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/nonessential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv'
        print(f'inp: {inp}')
        Figure_3_EXP_df = pd.read_csv(inp)
        Figure_3_EXP_df['label'] = 'EXP'

        if essential == True:
            inp = f'{self.slug_path}/Modeling_Odds_of_Misfolding/Regressions/Plots/AF/whole_proteome/essential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv'
        else:
            inp = f'{self.slug_path}/Modeling_Odds_of_Misfolding/Regressions/Plots/AF/whole_proteome/nonessential_ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv'
        print(f'inp: {inp}')
        Figure_3_AF_df = pd.read_csv(inp)
        Figure_3_AF_df['label'] = 'AF'     
        
        Figure_3_df = pd.concat([Figure_3_EXP_df, Figure_3_AF_df])
        print(f'Figure_3_df:\n{Figure_3_df}')

        # Save the figure 1a raw plot df
        Figure_3_outfile_csv = os.path.join(self.out_path, f'MainTextFigure3_essential-{essential}_Supp.csv')
        Figure_3_df.to_csv(Figure_3_outfile_csv)
        print(f'SAVED: {Figure_3_outfile_csv}')
        ######################################

        #######################################
        ## set up figure(s)
        fig_width_mm = 180  # Width in mm
        fig_height_mm = 150  # Height in mm
        fig, axes = plt.subplots(3, 4, figsize=(mm_to_inches(fig_width_mm), mm_to_inches(fig_height_mm)), gridspec_kw={'width_ratios': [1, 1, 1, 0.175]})
        #######################################

        #######################################
        ## plot data
        lines = []
        for buffi, buff in enumerate(['C', 'CD', 'CG']):
            buff_df = Figure_3_df[Figure_3_df['buff'] == buff]
            #print(buff_df)

            for label in ['EXP', 'AF']:
                label_df = buff_df[buff_df['label'] == label]
                #print(label)
                x = label_df['spa'].values

                ## Plot the OR
                y = label_df['OR'].values
                ylb = y - label_df['OR_lb'].values
                yub = label_df['OR_ub'].values - y
                yerr = [ylb, yub]
                line1 = axes[0, buffi].errorbar(x, y, yerr=yerr, fmt='o', capsize=3, markersize=3,  elinewidth=0.5,  markeredgewidth=0.5, label=label)
                #axes[0, buffi].set_yscale('log')
                axes[0, buffi].axhline(y=1.0, color='black', linestyle='--', linewidth=0.5)
                if buffi == 0:
                    lines += [line1]

                ## Plot the pvalues
                y = label_df['pvalues'].values
                axes[1, buffi].plot(x, y, marker='o', ls='none', markersize=3,  markeredgewidth=0.5, label=label)
                axes[1, buffi].set_yscale('log')
                axes[1, buffi].axhline(y=0.05, color='black', linestyle='--', linewidth=0.5)

                ## Plot the n values
                y = label_df['n'].values
                axes[2, buffi].plot(x, y, marker='o', ls='none', markersize=3,  markeredgewidth=0.5, label=label)
      
        #######################################

        ######################################
        ## adjust axes locations and labels
        panel_labels = np.asarray([['a', 'b', 'c'], ['d', 'e', 'f'], ['g', 'h', 'i']])
        panel_label_ys = {0:1, 1:0.66, 2:0.325}
        x0s = {0:0.075, 1:0.365, 2:2/3}
        y0s = {0:0.74, 1:0.4, 2:0.07}
        for coli, colname in enumerate(['cyto-serum', '+DnaK', '+GroEL']):
            for rowi, rowname in enumerate(['Odds ratio for association of\nmisfolding and entanglement regeion', 'p-value', 'Number of proteins in dataset']):

                # Adjust the linewidth of the axis spines
                for spine in axes[rowi, coli].spines.values():
                    spine.set_linewidth(0.5)  # Set the linewidth for all spines

                # Adjust the linewidth of the ticks
                axes[rowi, coli].tick_params(width=0.5)  # Both major and minor ticks

                axes[rowi, coli].set_ylabel(rowname)
                axes[rowi, coli].set_xlabel('Sum of peptide abundances')

                axes[rowi, coli].tick_params(axis='y', labelsize=6)
                axes[rowi, coli].tick_params(axis='x', labelsize=6)

                # Remove the right and top spines
                axes[rowi, coli].spines['right'].set_visible(False)
                axes[rowi, coli].spines['top'].set_visible(False)

                if rowi == 1:
                    axes[rowi, coli].set_ylim(top=1)
                if rowi == 0:
                    axes[rowi, coli].set_ylim(0.5, 2)

                axs_position = axes[rowi, coli].get_position()
                x0 = x0s[coli]
                y0 = y0s[rowi]
                width0, height0 = axs_position.extents[2] - axs_position.extents[0], axs_position.extents[3] - axs_position.extents[1]
                #print(axs_position, width0, height0)
                axes[rowi, coli].set_position([x0, y0, width0, height0])  # [left, bottom, width, height]

                bbox_in_fig_coords = axes[rowi, coli].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
                fig.text(bbox_in_fig_coords.x0, panel_label_ys[rowi], panel_labels[rowi, coli], fontsize=8, fontweight='bold', va='top', ha='left')
        #######################################

        axes[0, 3].set_axis_off()
        axes[1, 3].set_axis_off()
        axes[2, 3].set_axis_off()

        axes[0, 0].set_title('Cyto-serum', y=0.92)
        axes[0, 1].set_title('+Dnak', y=0.92)
        axes[0, 2].set_title('+GroEL', y=0.92)

        # Adjust layout for clarity
        labels = ['Crystal Strutures', 'AlphaFold']
        if essential:
            label = 'Essential'
        else:
            label = 'Non-essential'
        fig.legend(lines, labels, loc='upper right', bbox_to_anchor=(1.009, 1), title=label, frameon=False)

        #######################################
        figure_outpath = os.path.join(self.out_path, f'MainTextFigure3_essential-{essential}_Supp.pdf')
        plt.savefig(figure_outpath)
        print(f'SAVED: {figure_outpath}')

        figure_outpath = os.path.join(self.out_path, f'MainTextFigure3_essential-{essential}_Supp.png')
        plt.savefig(figure_outpath)
        print(f'SAVED: {figure_outpath}')

        figure_outpath = os.path.join(self.out_path, f'MainTextFigure3_essential-{essential}_Supp.svg')
        plt.savefig(figure_outpath)
        print(f'SAVED: {figure_outpath}')
        ####################################### 
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
    Creates Figure 3
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
    fig_width_mm = 180  # Width in mm
    fig_height_mm = 150  # Height in mm
    #custom_font_path = "/storage/group/epo2/default/ims86/miniconda3/envs/FtoF/fonts/Arial.ttf" # Path to your custom font
    #arial_font = fm.FontProperties(fname=custom_font_path) # Create a FontProperties object
    plt.rcParams['font.family'] = 'Arial'  # Change to your desired font, e.g., 'Times New Roman', 'DejaVu Sans', etc.
    plt.rcParams['font.size'] = 6  # Default font size
    plt.rcParams['pdf.fonttype'] = 42


    ## make subplot figures 
    plotter.plot_MainTextFigure3_Supp(essential=True)
    plotter.plot_MainTextFigure3_Supp(essential=False)
    print('NORMAL TERMINATION')

if __name__ == "__main__":
    main()

