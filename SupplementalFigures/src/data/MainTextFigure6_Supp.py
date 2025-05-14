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
    def plot_MisfoldedTrajSelection_Supp(self,):
        """
        Figure 6e: 
            Misfolding mechanism (fraction frames with only loss, only gain, both)
        """

        #######################################
        ## Load Figure 6c data
        #../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingMechanism/setID3/DATA/Mechanism_stats_summary_setID3.csv
        inp = f'{self.slug_path}/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingMechanism/setID3/DATA/Combined_and_Processed_threshold_metrics_setID3.csv'
        print(f'inp: {inp}')
        Figure_MisfoldedTrajSelection_df = pd.read_csv(inp)
        print(Figure_MisfoldedTrajSelection_df)

        noNative_MisfoldedTrajSelection_df = Figure_MisfoldedTrajSelection_df[Figure_MisfoldedTrajSelection_df['NativeFolded'] == 0]
        print(f'noNative_MisfoldedTrajSelection_df:\n{noNative_MisfoldedTrajSelection_df}')

        filteredQmode_MisfoldedTrajSelection_df = noNative_MisfoldedTrajSelection_df[noNative_MisfoldedTrajSelection_df['Qmode'] >= 0.6]
        print(f'filteredQmode_MisfoldedTrajSelection_df:\n{filteredQmode_MisfoldedTrajSelection_df}')

        filteredQmodeAndMisfProp_MisfoldedTrajSelection_df = filteredQmode_MisfoldedTrajSelection_df[filteredQmode_MisfoldedTrajSelection_df['MisfoldingProp'] >= 0.8]
        print(f'filteredQmodeAndMisfProp_MisfoldedTrajSelection_df:\n{filteredQmodeAndMisfProp_MisfoldedTrajSelection_df}')

        Figure_MisfoldedTrajSelection_outfile_csv = os.path.join(self.out_path, f'Figure_MisfoldedTrajSelection.csv')
        Figure_MisfoldedTrajSelection_df.to_csv(Figure_MisfoldedTrajSelection_outfile_csv)
        print(f'SAVED: {Figure_MisfoldedTrajSelection_outfile_csv}')
        #######################################

        #######################################
        ## set up figure(s)
        fig_width_mm = 180  # Width in mm
        fig_height_mm = 90  # Height in mm
        fig, axes = plt.subplots(2, 4, figsize=(mm_to_inches(fig_width_mm), mm_to_inches(fig_height_mm)), gridspec_kw={'width_ratios': [1, 1, 1, 1]})
        #######################################

        #######################################
        ## plot row 1 GQ data
        labels = {0:'All', 1:'Non-native', 2:'Non-native\nand folded', 3:'Non-native\nfolded and\npersistent'}
        for idx, df in {0:Figure_MisfoldedTrajSelection_df, 1:noNative_MisfoldedTrajSelection_df, 2:filteredQmode_MisfoldedTrajSelection_df, 3:filteredQmodeAndMisfProp_MisfoldedTrajSelection_df}.items():
            print(idx, df)
            G, Q, MisfoldingProp = df['Gmode'].values, df['Qmode'].values, df['MisfoldingProp'].values
            axes[0, idx].plot(Q, G, marker='o', ls='none', markersize=3,  markeredgewidth=0.5)
            axes[1, idx].plot(Q, MisfoldingProp, marker='o', ls='none', markersize=3,  markeredgewidth=0.5)

            n = len(df['Gmode'].values)
            axes[1, idx].text(0.05, 0.05, f"{labels[idx]}\n(n={n})", transform=axes[1, idx].transAxes, va='bottom', ha='left', fontsize=7) 
            
        #######################################

        ######################################
        ## adjust axes locations and labels
        panel_labels = np.asarray([['a', 'b', 'c', 'd'], ['e', 'f', 'g', 'h']])
        panel_label_ys = {0:1, 1:0.5}
        x0s = {0:0.075, 1:0.325, 2:0.575, 3:0.815}
        y0s = {0:0.595, 1:0.11}
        for coli, colname in {0:'All Trajectories', 1:'Non-nativly folded trajectories', 2:'Compact and Non-natively folded trajectories', 3:'Compact, persistent, and Non-natively folded trajectories'}.items():
            for rowi, rowname in enumerate([(r'$G_{mode}$', r'$Q_{mode}$'), (r'$P_{misfolded}$',r'$Q_{mode}$')]):

                # Adjust the linewidth of the axis spines
                for spine in axes[rowi, coli].spines.values():
                    spine.set_linewidth(0.5)  # Set the linewidth for all spines

                # Adjust the linewidth of the ticks
                axes[rowi, coli].tick_params(width=0.5)  # Both major and minor ticks

                axes[rowi, coli].set_ylabel(rowname[0])
                axes[rowi, coli].set_xlabel(rowname[1])

                axes[rowi, coli].tick_params(axis='y', labelsize=6)
                axes[rowi, coli].tick_params(axis='x', labelsize=6)

                # Remove the right and top spines
                axes[rowi, coli].spines['right'].set_visible(False)
                axes[rowi, coli].spines['top'].set_visible(False)

                if rowi == 1:
                    axes[rowi, coli].set_ylim(0.0,1)
                    #axes[rowi, coli].set_yscale('log')
                    axes[rowi, coli].set_xlim(0.0,1)
                
                if rowi == 0:
                    axes[rowi, coli].set_xlim(0.0,1)

                axs_position = axes[rowi, coli].get_position()
                x0 = x0s[coli]
                y0 = y0s[rowi]
                width0, height0 = axs_position.extents[2] - axs_position.extents[0], axs_position.extents[3] - axs_position.extents[1]
                print(axs_position, width0, height0)
                width0, height0 = 0.15, 0.35
                axes[rowi, coli].set_position([x0, y0, width0, height0])  # [left, bottom, width, height]

                bbox_in_fig_coords = axes[rowi, coli].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
                fig.text(bbox_in_fig_coords.x0, panel_label_ys[rowi], panel_labels[rowi, coli], fontsize=8, fontweight='bold', va='top', ha='left')
        #######################################


        # Adjust layout for clarity
        #labels = ['AlphaFold', 'Crystal Strutures']
        #fig.legend(lines, labels, loc='upper right', bbox_to_anchor=(1.009, 1), title="Dataset", frameon=False)

        #######################################
        figure_outpath = os.path.join(self.out_path, 'MainTextFigure6_MisfoldedTrajSelection_Supp.pdf')
        plt.savefig(figure_outpath)
        print(f'SAVED: {figure_outpath}')

        figure_outpath = os.path.join(self.out_path, 'MainTextFigure6_MisfoldedTrajSelection_Supp.png')
        plt.savefig(figure_outpath)
        print(f'SAVED: {figure_outpath}')

        figure_outpath = os.path.join(self.out_path, 'MainTextFigure6_MisfoldedTrajSelection_Supp.svg')
        plt.savefig(figure_outpath)
        print(f'SAVED: {figure_outpath}')
        ####################################### 
        # quit()       
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
    fig_width_mm = 183  # Width in mm
    fig_height_mm = 100  # Height in mm
    #custom_font_path = "/storage/group/epo2/default/ims86/miniconda3/envs/FtoF/fonts/Arial.ttf" # Path to your custom font
    #arial_font = fm.FontProperties(fname=custom_font_path) # Create a FontProperties object
    plt.rcParams['font.family'] = 'Arial'  # Change to your desired font, e.g., 'Times New Roman', 'DejaVu Sans', etc.
    plt.rcParams['font.size'] = 7  # Default font size
    plt.rcParams['pdf.fonttype'] = 42
    #plt.rc('text', usetex=True)

    plotter.plot_MisfoldedTrajSelection_Supp()

    print('NORMAL TERMINATION')

if __name__ == "__main__":
    main()

