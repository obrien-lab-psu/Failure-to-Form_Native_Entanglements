import os
import glob
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib import font_manager as fm
import matplotlib as mpl
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import PathPatch
from matplotlib.path import Path
import matplotlib.image as mpimg
from PIL import Image
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
import matplotlib.image as image
#os.environ["PYMOL_LICENSE_FILE"] = "/storage/home/ims86/pymolLicenseFile.lic"
#import pymol

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

    Figure is 1 row by 4 column 
    Figure 6a (row 1 column 1)
        Schematic of CG -> temperature quenching simulations
        (All atom) ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Rebuild_AllAtom_structures/PDBs/P0AD61_4YNG_C_rebuilt.pdb
        (CG cor) ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/P0AD61_4YNG_C/setup/P0AD61_4YNG_C_rebuilt_clean_ca.cor
        (CG psf) ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/P0AD61_4YNG_C/setup/P0AD61_4YNG_C_rebuilt_clean_ca.psf
        (CG unfolded) ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/P0AD61_4YNG_C/Unfolding/P0AD61_4YNG_C_t48_unfolding_finalframe1433.pdb
        (CG misfolded) 

    Figure 6b (row 1 column 2): 
        Figure 6b (row 1 column 1): 
            Misfolding propensity 
            (dataset 1, All)../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingPropensity/All/Plots/Fract_misfolded_set1_NativeByRef_Stats_Scope-full_plot.csv        
            (dataset 2, All)../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingPropensity/All/Plots/Fract_misfolded_set2_NativeByRef_Stats_Scope-full_plot.csv
            (comp 1 and 2 stats, All)../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingPropensity/All/DATA/Fract_misfolded_permutation_TwoSampleStats_Scope-full_NativeByRef.csv
            
    Figure 6c (row 1 column 3): 
        Figure 6c (row 1 column 1): 
            Misfolding propensity 
            (dataset 1, Reduced)../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingPropensity/Reduced/Plots/Fract_misfolded_set1_NativeByRef_Stats_Scope-full_plot.csv        
            (dataset 2, Reduced)../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingPropensity/Reduced/Plots/Fract_misfolded_set2_NativeByRef_Stats_Scope-full_plot.csv
            (comp 1 and 2 stats, Reduced)../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingPropensity/Reduced/DATA/Fract_misfolded_permutation_TwoSampleStats_Scope-full_NativeByRef.csv
            
    Figure 6d (row 2 column 1):
        Examples of misfolded structures
        Loss example # 1
            (Native) ../Simulations_of_Native_Entanglement_Misfolding/Rebuild_AllAtom_structures/data/post_rebuilt/P0AES0_2IO9_B_rebuilt.pdb
            (Misfolded) 
            (EntInfo) ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/P0AES0_2IO9_B/Cluster_ChangesInEnt/P0AES0_2IO9_B_t29_clustered.EntInfo
            
    
    Figure 6e: 
        Figure 6c (row 2 column 1): 
            Misfolding mechanism
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingMechanism/setID3/DATA/Mechanism_stats_summary_setID3.csv 

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
    def plot_Figure_6a(self, ax):
        """
        Figure 6a (row 0 column 0)
            Temperature quenching simulation and misfolding mechanism schematics
        """
        ax.axis('tight')
        ax.axis('off')
    #################################################################################################################

    #################################################################################################################
    def plot_Figure_6b(self, ax):
        """
        Figure 6b (row 1 column 0): 
            Loss example
        """
        ax.axis('tight')
        ax.axis('off')
    #################################################################################################################

    #################################################################################################################
    def plot_Figure_6b(self, ax):
        """
        Figure 6b (row 0 column 1): 
            Misfolding propensity 
            (dataset 1, Reduced)../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingPropensity/Reduced/DATA/Fract_misfolded_setID1_OneSampleStats_full_NativeByRef.csv        
            (dataset 2, Reduced)../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingPropensity/Reduced/DATA/Fract_misfolded_setID2_OneSampleStats_full_NativeByRef.csv
            (comp 1 and 2 stats, Reduced)../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingPropensity/Reduced/DATA/Fract_misfolded_permutation_TwoSampleStats_Scope-full_NativeByRef.csv

            
        """
        buff_tag = {'C': 'cyto-serum', 'CD': '+DnaK', 'CG': '+GroEL'}
        #######################################
        ## Load Figure 6b data
        #../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingMechanism/setID3/DATA/Mechanism_stats_summary_setID3.csv
        inp = f'{self.slug_path}/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingPropensity/Reduced/DATA/Fract_misfolded_setID1_OneSampleStats_full_NativeByRef.csv'
        print(f'inp: {inp}')
        Figure_6b_Reduced_set1_df = pd.read_csv(inp)
        Figure_6b_Reduced_set1_df = Figure_6b_Reduced_set1_df[Figure_6b_Reduced_set1_df['Label'] == '90-100%']
        Figure_6b_Reduced_set1_df['xlabel'] = 'Low'
        print(f'Figure_6b_Reduced_set1_df:\n{Figure_6b_Reduced_set1_df}')

        inp = f'{self.slug_path}/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingPropensity/Reduced/DATA/Fract_misfolded_setID2_OneSampleStats_full_NativeByRef.csv'
        print(f'inp: {inp}')
        Figure_6b_Reduced_set2_df = pd.read_csv(inp)
        Figure_6b_Reduced_set2_df = Figure_6b_Reduced_set2_df[Figure_6b_Reduced_set2_df['Label'] == '90-100%']
        Figure_6b_Reduced_set2_df['xlabel'] = 'High'
        print(f'Figure_6b_Reduced_set2_df:\n{Figure_6b_Reduced_set2_df}')

        inp = f'{self.slug_path}/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingPropensity/Reduced/DATA/Fract_misfolded_permutation_TwoSampleStats_Scope-full_NativeByRef.csv'
        print(f'inp: {inp}')
        Figure_6b_Reduced_comp_df = pd.read_csv(inp)
        Figure_6b_Reduced_comp_df = Figure_6b_Reduced_comp_df[Figure_6b_Reduced_comp_df['Label'] == '90-100%']
        Figure_6b_Reduced_comp_df['Full'] = False
        print(f'Figure_6b_Reduced_comp_df:\n{Figure_6b_Reduced_comp_df}')     

        Figure_6b_Reduced_set1_df['perm_stat'] = Figure_6b_Reduced_comp_df['statistic']
        Figure_6b_Reduced_set1_df['pvalue'] = Figure_6b_Reduced_comp_df['pvalue']
        Figure_6b_Reduced_set2_df['perm_stat'] = Figure_6b_Reduced_comp_df['statistic']
        Figure_6b_Reduced_set2_df['pvalue'] = Figure_6b_Reduced_comp_df['pvalue']

        Figure_6b_df = pd.concat([Figure_6b_Reduced_set1_df, Figure_6b_Reduced_set2_df])  
        Figure_6b_df = Figure_6b_df.drop('ttest', axis=1)  
        print(f'Figure_6b_df:\n{Figure_6b_df}')

        Figure_6b_outfile_csv = os.path.join(self.out_path, f'Figure_6b.csv')
        Figure_6b_df.to_csv(Figure_6b_outfile_csv)
        print(f'SAVED: {Figure_6b_outfile_csv}')
        #######################################


        # Adjust the linewidth of the axis spines
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)  # Set the linewidth for Reduced spines

        # Adjust the linewidth of the ticks
        ax.tick_params(width=0.5)  # Both major and minor ticks


        # Calculate error bars
        yerr_lower = Figure_6b_df['avg'] - Figure_6b_df['lb']
        yerr_upper = Figure_6b_df['ub'] - Figure_6b_df['avg']
        yerr = [yerr_lower, yerr_upper]
        
        # Adjust to be fraction of misfolded frames
        Figure_6b_df['avg'] = Figure_6b_df['avg']
        #print(f'Figure_6b_df:\n{Figure_6b_df}')

        # Plot the trace with error bars
        ax.errorbar(
            Figure_6b_df['xlabel'], 
            Figure_6b_df['avg'], 
            yerr=yerr, 
            fmt='o',  # Line and point markers
            capsize=3,  # Add caps to the error bars
            color='black',
            markersize=3,
            markerfacecolor='black',  
            elinewidth=0.5,  
            markeredgewidth=0.5)

                
        # Add a dashed black line at y=1    
        ax.set_ylabel('Simulation observed\nmisfolding propensity')
        ax.set_xticklabels(Figure_6b_df['xlabel'])  # Rotate x-axis labels for better readability
        ax.set_xlabel('Experimentaly observed\nmisfolding propensity')

        # Customize the plot
        ax.set_ylim(0.0, 1.0)
        ax.set_xlim(-0.5, 1.5)

        ax.tick_params(axis='y', labelsize=6)
        ax.tick_params(axis='x', labelsize=6)
        # Remove the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        # add the curly brackets
        ax.annotate("***", xy=(0.5, 0.775), xytext=(0.5, 0.825), xycoords='axes fraction', arrowprops=dict(arrowstyle='-[, widthB=3.25, lengthB=0.7', lw=0.5, color='k'), va='center', ha='center', fontweight='bold')
    #################################################################################################################

    #################################################################################################################
    def plot_Figure_6c(self, ax):
        """
        Figure 6c (row 1 column 0): 
            Loss example
        """
        ax.axis('tight')
        ax.axis('off')
    #################################################################################################################

    #################################################################################################################
    def plot_Figure_6d(self, ax):
        """
        Figure 6d (row 1 column 0): 
            Loss example
        """
        ax.axis('tight')
        ax.axis('off')
    #################################################################################################################


    #################################################################################################################
    def plot_Figure_6e(self, ax):
        """
        Figure 6e: 
            Misfolding mechanism (fraction frames with only loss, only gain, both)
            ../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingMechanism/setID3/DATA/Mechanism_stats_summary_setID3.csv

                    """

        #######################################
        ## Load Figure 6c data
        #../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingMechanism/setID3/DATA/Mechanism_stats_summary_setID3.csv
        inp = f'{self.slug_path}/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingMechanism/setID3/DATA/Mechanism_stats_summary_setID3.csv'
        print(f'inp: {inp}')
        Figure_6e_df = pd.read_csv(inp)
        print(Figure_6e_df)
        Figure_6e_df = Figure_6e_df[Figure_6e_df['type'].isin(['LossOnly', 'GainOnly', 'BothLossGain'])]
        Figure_6e_df['xlabel'] = ['Short\nlived']*3 + ['Long\nlived']*3
        #Figure_6e_df['xlabel'] = ['Short,\n(<100 ns)']*3 + ['Long,\n(>=100 ns)']*3
        Figure_6e_df['color'] = ['red', 'blue', 'purple']*2
        Figure_6e_df['markers'] = ['o', '^', 's']*2
        Figure_6e_df['type'] = ['Only loss', 'Only gain', 'Both']*2
        print(f'Figure_6e_df:\n{Figure_6e_df}')

        Figure_6e_outfile_csv = os.path.join(self.out_path, f'Figure_6e.csv')
        Figure_6e_df.to_csv(Figure_6e_outfile_csv)
        print(f'SAVED: {Figure_6e_outfile_csv}')
        #######################################

        # Adjust the linewidth of the axis spines
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)  # Set the linewidth for Reduced spines

        # Adjust the linewidth of the ticks
        ax.tick_params(width=0.5)  # Both major and minor ticks

        for changetype, changetype_df in Figure_6e_df.groupby('type'):
            print(changetype_df)
            # Calculate error bars
            yerr_lower = changetype_df['mean'] - changetype_df['lb']
            yerr_upper = changetype_df['ub'] - changetype_df['mean']
            yerr = [yerr_lower, yerr_upper]

            # Plot the trace with error bars
            ax.errorbar(
                changetype_df['xlabel'], 
                changetype_df['mean'], 
                yerr=yerr, 
                fmt='o',  # Line and point markers
                capsize=3,  # Add caps to the error bars
                color=changetype_df['color'].values[0],
                marker=changetype_df['markers'].values[0],
                markersize=3,
                elinewidth=0.5,  
                markeredgewidth=0.5, label=changetype)

                
        # Add a dashed black line at y=1    
        ax.set_ylabel('Probability of structure\nmisfolded by mechanism')
        #ax.set_xlabel('Misfolded state lifetime, (ns)')

        # Customize the plot
        ax.set_ylim(0.0, 1.0)
        ax.set_xlim(-0.5, 1.5)

        ax.tick_params(axis='y', labelsize=6)
        ax.tick_params(axis='x', labelsize=6)
        # Remove the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        legend = ax.legend(frameon=True, loc='upper right', bbox_to_anchor=(1, 1), fontsize=5)
        legend.get_frame().set_linewidth(0.5)  # Set frame linewidth
        legend.get_frame().set_edgecolor('black')  # Change frame background color        
        # add the curly brackets
        #ax.annotate("***", xy=(0.5, 0.75), xytext=(0.5, 0.77), xycoords='axes fraction', arrowprops=dict(arrowstyle='-[, widthB=2.25, lengthB=0.75', lw=0.5, color='k'), va='center', ha='center', fontweight='bold')
    #################################################################################################################

    #################################################################################################################
    def plot_Figure_6f(self, ax):
        """
        Figure 6f: 
            Misfolding mechanism (fraction frames with only loss, only gain, both)
            ../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingMechanism/setID3/DATA/Mechanism_stats_summary_setID3.csv
        """

        #######################################
        ## Load Figure 6f data
        #../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingMechanism/setID3/DATA/Mechanism_stats_summary_setID3.csv
        inp = f'{self.slug_path}/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingMechanism/setID3/DATA/Mechanism_stats_summary_setID3.csv'
        print(f'inp: {inp}')
        Figure_6f_df = pd.read_csv(inp)
        print(Figure_6f_df)
        Figure_6f_df = Figure_6f_df[Figure_6f_df['type'].isin(['TrueLossOverlap', 'TrueLossOverlap_minusOne'])]
        Figure_6f_df['xlabel'] = ['Short\nlived']*2 + ['Long\nlived']*2
        #Figure_6f_df['xlabel'] = ['Short,\n(<100 ns)']*2 + ['Long,\n(>=100 ns)']*2
        Figure_6f_df['color'] = ['red', 'black']*2
        Figure_6f_df['markers'] = ['o', 's']*2
        Figure_6f_df['type'] = ['Paired', 'Unpaired']*2
        print(f'Figure_6f_df:\n{Figure_6f_df}')

        Figure_6f_outfile_csv = os.path.join(self.out_path, f'Figure_6f.csv')
        Figure_6f_df.to_csv(Figure_6f_outfile_csv)
        print(f'SAVED: {Figure_6f_outfile_csv}')
        #######################################

        # Adjust the linewidth of the axis spines
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)  # Set the linewidth for Reduced spines

        # Adjust the linewidth of the ticks
        ax.tick_params(width=0.5)  # Both major and minor ticks

        for changetype, changetype_df in Figure_6f_df.groupby('type'):
            print(changetype_df)
            # Calculate error bars
            yerr_lower = changetype_df['mean'] - changetype_df['lb']
            yerr_upper = changetype_df['ub'] - changetype_df['mean']
            yerr = [yerr_lower, yerr_upper]

            # Plot the trace with error bars
            ax.errorbar(
                changetype_df['xlabel'], 
                changetype_df['mean'], 
                yerr=yerr, 
                fmt='o',  # Line and point markers
                capsize=3,  # Add caps to the error bars
                color=changetype_df['color'].values[0],
                marker=changetype_df['markers'].values[0],
                markersize=3,
                elinewidth=0.5,  
                markeredgewidth=0.5, label=changetype)

                
        # Add a dashed black line at y=1    
        ax.set_ylabel('Probability of entanglement\nloss paired with a gain')
        #ax.set_xlabel('Misfolded state lifetime, (ns)')

        # Customize the plot
        ax.set_ylim(0.0, 1.0)
        ax.set_xlim(-0.5, 1.5)

        ax.tick_params(axis='y', labelsize=6)
        ax.tick_params(axis='x', labelsize=6)
        # Remove the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        legend = ax.legend(frameon=True, loc='upper right', bbox_to_anchor=(1, 1), fontsize=5)
        legend.get_frame().set_linewidth(0.5)  # Set frame linewidth
        legend.get_frame().set_edgecolor('black')  # Change frame background color
        # add the curly brackets
        #ax.annotate("***", xy=(0.5, 0.75), xytext=(0.5, 0.77), xycoords='axes fraction', arrowprops=dict(arrowstyle='-[, widthB=2.25, lengthB=0.75', lw=0.5, color='k'), va='center', ha='center', fontweight='bold')
    #################################################################################################################

def mm_to_inches(mm):
    return mm / 25.4

# Function to compute the CDF
def compute_cdf(data):
    sorted_data = np.sort(data)
    cdf = np.arange(1, len(sorted_data) + 1) / len(sorted_data)
    return sorted_data, cdf

def format_scientific(number, precision=2):
    # Split the number into mantissa and exponent
    formatted = f"{number:.{precision}e}"
    mantissa, exponent = formatted.split("e")
    exponent = int(exponent)
    #print(mantissa, exponent)

    # Convert the exponent into superscript
    superscript_map = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
    exponent_superscript = str(exponent).translate(superscript_map)
    # Return the formatted string
    formatted_str = f"{mantissa} × $10^{{{exponent}}}$"
    print(formatted_str)
    return formatted_str

def main():
    """
    Creates Figure 6
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
    fig_width_mm = 183  # Width in mm
    fig_height_mm = 100  # Height in mm
    #custom_font_path = "/storage/group/epo2/default/ims86/miniconda3/envs/FtoF/fonts/Arial.ttf" # Path to your custom font
    #arial_font = fm.FontProperties(fname=custom_font_path) # Create a FontProperties object
    plt.rcParams['font.family'] = 'Arial'  # Change to your desired font, e.g., 'Times New Roman', 'DejaVu Sans', etc.
    plt.rcParams['font.size'] = 6  # Default font size
    plt.rcParams['pdf.fonttype'] = 42
    #plt.rc('text', usetex=True)

    ### use gridspec. give lower quality images for some reason?
    
    fig = plt.figure(figsize=(mm_to_inches(fig_width_mm), mm_to_inches(fig_height_mm)), dpi=600, constrained_layout=True)
    gs = gridspec.GridSpec(2, 4, figure=fig)
    # Add subplots
    ax1 = fig.add_subplot(gs[0, 0:3])  # 5a
    ax2 = fig.add_subplot(gs[0, 3])  # 5b
    ax3 = fig.add_subplot(gs[1, 0])  # 5c
    ax4 = fig.add_subplot(gs[1, 1])  # 5d
    ax5 = fig.add_subplot(gs[1, 2])  # 5e
    ax6 = fig.add_subplot(gs[1, 3])  # 5f

    axs = [ax1, ax2, ax3, ax4, ax5, ax6]

    # make subplot figures 
    plotter.plot_Figure_6a(axs[0])
    plotter.plot_Figure_6b(axs[1])
    plotter.plot_Figure_6c(axs[2])
    plotter.plot_Figure_6d(axs[3])
    plotter.plot_Figure_6e(axs[4])
    plotter.plot_Figure_6f(axs[5])
    

    ########################################################
    # Get the current position of the subplots and adjust their positions manually
    ## panel a
    axs0_position = axs[0].get_position()
    width0, height0 = axs0_position.extents[2] - axs0_position.extents[0], axs0_position.extents[3] - axs0_position.extents[1]
    print('a', axs0_position, width0, height0)
    axs[0].set_position([0.001, 0.575, 0.49, 0.305])  # [left, bottom, width, height]
    axs[0].axis('off')

    bbox_in_fig_coords = axs[0].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bbox_in_fig_coords.x0, 0.995, 'a', fontsize=8, fontweight='bold', va='top', ha='left')

    #####################################################
    ## panel b
    axs1_position = axs[1].get_position()
    width1, height1 = axs1_position.extents[2] - axs1_position.extents[0], axs1_position.extents[3] - axs1_position.extents[1]
    print('b', axs1_position, width1, height1)
    axs[1].set_position([.815, 0.625, 0.15, 0.305])

    bbox_in_fig_coords = axs[1].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bbox_in_fig_coords.x0, 0.995, 'b', fontsize=8, fontweight='bold', va='top', ha='left')

    #####################################################
    ## panel c (row 1 column 3)
    axs2_position = axs[2].get_position()
    width2, height2 = axs2_position.extents[2] - axs2_position.extents[0], axs2_position.extents[3] - axs2_position.extents[1]
    print('c', axs2_position, width2, height2)
    axs[2].set_position([0.001, 0.15, 0.1, 0.305])
    axs[2].axis('off')

    bbox_in_fig_coords = axs[2].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bbox_in_fig_coords.x0, 0.995, 'c', fontsize=8, fontweight='bold', va='top', ha='left')

    #####################################################
    ##panel d
    axs3_position = axs[3].get_position()
    width3, height3 = axs3_position.extents[2] - axs3_position.extents[0], axs3_position.extents[3] - axs3_position.extents[1]
    print('d', axs3_position, width3, height3)
    axs[3].set_position([0.25, 0.1, 0.175, 0.305])
    #axs[3].axis('off')

    bbox_in_fig_coords = axs[3].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bbox_in_fig_coords.x0, 0.5, 'd', fontsize=8, fontweight='bold', va='top', ha='left')

    #####################################################
    ##panel e
    axs4_position = axs[4].get_position()
    width4, height4 = axs4_position.extents[2] - axs4_position.extents[0], axs4_position.extents[3] - axs4_position.extents[1]
    print('e', axs4_position, width4, height4)
    axs[4].set_position([0.565, 0.1, 0.15, 0.375])

    bbox_in_fig_coords = axs[4].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bbox_in_fig_coords.x0, 0.5, 'e', fontsize=8, fontweight='bold', va='top', ha='left')

    #####################################################
    ##panel f
    axs5_position = axs[5].get_position()
    width5, height5 = axs5_position.extents[2] - axs5_position.extents[0], axs5_position.extents[3] - axs5_position.extents[1]
    print('f', axs5_position, width5, height5)
    axs[5].set_position([0.815, 0.1, 0.15, 0.375])

    bbox_in_fig_coords = axs[5].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bbox_in_fig_coords.x0, 0.5, 'f', fontsize=8, fontweight='bold', va='top', ha='left')


    # final formating and output
    # Automatically adjust layout
    #fig.tight_layout()  # 'pad' is the overall padding
    figure_outpath = os.path.join(args.out_path, 'Figure6.pdf')
    plt.savefig(figure_outpath)
    print(f'SAVED: {figure_outpath}')

    figure_outpath = os.path.join(args.out_path, 'Figure6.png')
    plt.savefig(figure_outpath)
    print(f'SAVED: {figure_outpath}')

    figure_outpath = os.path.join(args.out_path, 'Figure6.svg')
    plt.savefig(figure_outpath)
    print(f'SAVED: {figure_outpath}')

    print('NORMAL TERMINATION')

if __name__ == "__main__":
    main()

