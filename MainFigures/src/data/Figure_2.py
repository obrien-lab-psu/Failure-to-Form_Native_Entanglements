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

    Figure is 2 row by 4 column 
    Figure 2a (row 1 column 1)


    Figure 2b (row 1 column 2): 
        Pie charts showing counts of proteins that passed our filters in this study

    
    Figure 2c (row 1 column 3): 

        
    
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
    def plot_Figure_2a(self, ax):
        """
        Figure 5a

        """

    #################################################################################################################

    #################################################################################################################
    def plot_Figure_2b(self, ax):
        """
        Figure 5b (row 1 column 1): 
            ../Make_Protein_Feature_Files/data/Combined_genelist_counts_EXP.csv

        """
        buff_tag = {'C': 'cyto-serum', 'CD': '+DnaK', 'CG': '+GroEL'}
        #######################################
        ## Load Figure 2b data
        inp = f'../Make_Protein_Feature_Files/data/Combined_genelist_counts_EXP.csv'
        print(f'inp: {inp}')
        Figure_2b_df = pd.read_csv(inp)
        Figure_2b_df = Figure_2b_df[(Figure_2b_df['timepoint'] == 'Rall') & (Figure_2b_df['spa_threshold'] == 50) & (Figure_2b_df['LiPMScov_threshold'] == 50)]
        Figure_2b_df = Figure_2b_df[['buff', 'essential_ent_genes_n', 'essential_nonent_genes_n', 'nonessential_ent_genes_n', 'nonessential_nonent_genes_n']]
        Figure_2b_df = Figure_2b_df.rename(columns={'essential_ent_genes_n':'Ess + Ent', 'essential_nonent_genes_n':'Ess - Ent', 'nonessential_ent_genes_n':'NonEss + Ent', 'nonessential_nonent_genes_n':'NonEss - Ent'})
        custom_order = ["C", "CD", "CG"]
        Figure_2b_df["buff"] = pd.Categorical(Figure_2b_df["buff"], categories=custom_order, ordered=True)
        Figure_2b_df = Figure_2b_df.sort_values("buff")
        Figure_2b_df['buff'] = [buff_tag[buff] for buff in Figure_2b_df['buff']]
        print(f'Figure_2b_df:\n{Figure_2b_df}')
        print(Figure_2b_df.keys())
        #######################################


        # Save the Figure 2b raw plot df
        Figure_2b_outfile_csv = os.path.join(self.out_path, f'Figure_2b.csv')
        Figure_2b_df.to_csv(Figure_2b_outfile_csv)
        print(f'SAVED: {Figure_2b_outfile_csv}')

        # Loop through each row and create a pie chart
        colors = ["#37B2E6", "#B7E3F6", "#EC5B69", "#FAD4D8"]
        #colors = ["blue", "purple", "red", "pink"]
        labels = ["Ess + Ent", "Ess - Ent", "NonEss + Ent", "NonEss - Ent"]
        x_offset = 0  # Initial x offset for positioning
        for _, row in Figure_2b_df.iterrows():
            # Exclude the "buff" column for the pie chart values
            pie_values = row.drop("buff")
            print(pie_values)
            labels = pie_values.index
            print(labels)
            
            # Create a pie chart and position it using offsets
            wedges, texts, autotexts = ax.pie(pie_values, colors=colors, autopct=lambda pct: f"{int(round(pct/100.0 * sum(pie_values)))}", radius=1.0, center=(x_offset, 0))

            # Add label from "buff" below the corresponding pie chart
            ax.text(x_offset, -1.2, row["buff"], ha="center", fontsize=7)

            x_offset += 2.1  # Increase x offset for next pie chart

        # Remove axes for a clean look
        # Add a single legend for all pie charts
        ax.legend(labels, loc="upper left", bbox_to_anchor=(1, 1))
        ax.set_xlim(-1, x_offset - 1)
        ax.axis('off')
    #################################################################################################################

    #################################################################################################################
    def plot_Figure_2c(self, ax):
        """
        Figure 5c (row 1 column 1): 

        """

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
    Creates Figure 2
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
    fig = plt.figure(figsize=(mm_to_inches(fig_width_mm), mm_to_inches(fig_height_mm)))
    gs = gridspec.GridSpec(2, 4, figure=fig)
    # Add subplots
    ax1 = fig.add_subplot(gs[0, :])  # First row, first column
    ax2 = fig.add_subplot(gs[1, :2])  # First row, second column
    ax3 = fig.add_subplot(gs[1, 2:])  # First row, second column
    axs = [ax1, ax2, ax3]

    ## make subplot figures 
    plotter.plot_Figure_2a(axs[0])
    plotter.plot_Figure_2b(axs[1])
    plotter.plot_Figure_2c(axs[2])

    #########################################################
    # Get the current position of the subplots and adjust their positions manually
    ## panel a
    axs0_position = axs[0].get_position()
    width0, height0 = axs0_position.extents[2] - axs0_position.extents[0], axs0_position.extents[3] - axs0_position.extents[1]
    print('a', axs0_position, width0, height0)
    #axs[0].set_position([0.005, 0.625, 0.5, 0.325])  # [left, bottom, width, height]
    axs[0].axis('off')

    bbox_in_fig_coords = axs[0].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bbox_in_fig_coords.x0, 0.995, 'a', fontsize=8, fontweight='bold', va='top', ha='left')

    ## panel b
    axs1_position = axs[1].get_position()
    width1, height1 = axs1_position.extents[2] - axs1_position.extents[0], axs1_position.extents[3] - axs1_position.extents[1]
    print('b', axs1_position, width1, height1)
    axs[1].set_position([-0.075, 0.05, 0.65, 0.325])

    bbox_in_fig_coords = axs[1].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bbox_in_fig_coords.x0, 0.5, 'b', fontsize=8, fontweight='bold', va='top', ha='left')

    ## panel c
    axs2_position = axs[2].get_position()
    width2, height2 = axs2_position.extents[2] - axs2_position.extents[0], axs2_position.extents[3] - axs2_position.extents[1]
    print('c', axs2_position, width2, height2)
    axs[2].set_position([0.75, 0.1, 0.135, 0.325])

    bbox_in_fig_coords = axs[2].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bbox_in_fig_coords.x0, 0.5, 'c', fontsize=8, fontweight='bold', va='top', ha='left')
    axs[2].axis('off')

    #########################################################

    # final formating and output
    # Automatically adjust layout
    #fig.tight_layout()  # 'pad' is the overall padding
    figure_outpath = os.path.join(args.out_path, 'Figure2.pdf')
    plt.savefig(figure_outpath)
    print(f'SAVED: {figure_outpath}')

    figure_outpath = os.path.join(args.out_path, 'Figure2.png')
    plt.savefig(figure_outpath)
    print(f'SAVED: {figure_outpath}')

    figure_outpath = os.path.join(args.out_path, 'Figure2.svg')
    plt.savefig(figure_outpath)
    print(f'SAVED: {figure_outpath}')


    print('NORMAL TERMINATION')

if __name__ == "__main__":
    main()

