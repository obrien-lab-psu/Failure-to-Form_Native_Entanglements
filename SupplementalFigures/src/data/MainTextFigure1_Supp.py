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
import matplotlib.gridspec as gridspec
import pickle
from matplotlib import cm
from matplotlib_venn import venn2, venn3
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
    def plot_LiPMS_Supp(self,):
        """

        """
        buff_tag = {'C': 'cyto-serum', 'CD': '+DnaK', 'CG': '+GroEL'}
        #######################################
        ## Load Figure 1a data
        inp = f'{self.slug_path}/Processing_LiP-MS_data/SPA_thresholds/Data/spa_threshold_master.pkl'
        print(f'inp: {inp}')
        with open(inp, 'rb') as fh:
            SPA_data = pickle.load(fh)
        print(f'LOADED: {inp}')

        inp = f'{self.slug_path}/Processing_LiP-MS_data/COV_thresholds/Data/cov_threshold_master.pkl'
        print(f'inp: {inp}')
        with open(inp, 'rb') as fh:
            COV_data = pickle.load(fh)
        print(f'LOADED: {inp}')

        ## calcualte the CDF
        plot_dfs = []
        for key in SPA_data.keys():
            #print(key)
            buff, timepoint = key
            SPA_key = [spa_key for spa_key in SPA_data[key].keys() if spa_key[0] == 0][0]
            COV_key = [spa_key for spa_key in COV_data[key].keys() if spa_key[0] == 0][0]
            #print(SPA_key, COV_key)
            all_SPA_data = SPA_data[key][SPA_key]
            all_COV_data = COV_data[key][COV_key] 
            #print(key, buff, timepoint, len(all_SPA_data), len(all_COV_data))
            for cov in [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]:
                thresh_COV_data = all_COV_data[all_COV_data['coverage'] >= cov]
                cov_SPA_data = all_SPA_data[all_SPA_data['Accession'].isin(thresh_COV_data['Accession'])]
                sorted_data = np.sort(cov_SPA_data['N_spa'].values)
                x = np.unique(sorted_data)  # Unique sorted values
                y = np.searchsorted(sorted_data, x, side='right') / len(sorted_data)  # Cumulative probabilities
                
                cov_df = pd.DataFrame({'SPA':x, 'CDF':y})
                cov_df['buff'] = buff
                cov_df['timepoint'] = timepoint
                cov_df['cov'] = cov
                #print(cov_df)

                plot_dfs += [cov_df]

        plot_dfs = pd.concat(plot_dfs)
        print(f'plot_dfs:\n{plot_dfs}')
        # Save the figure 1a raw plot df
        plot_dfs_outfile_csv = os.path.join(self.out_path, f'MainTextFigure1_LiPMS_CDFs_Supp.csv')
        plot_dfs.to_csv(plot_dfs_outfile_csv)
        print(f'SAVED: {plot_dfs_outfile_csv}')
        #######################################
        
        #######################################
        ## set up figure(s)
        fig_width_mm = 180  # Width in mm
        fig_height_mm = 170  # Height in mm
        fig, axes = plt.subplots(4, 3, figsize=(mm_to_inches(fig_width_mm), mm_to_inches(fig_height_mm)))
        #######################################

        #######################################
        ## plot data
        rownames = []
        colormap = cm.get_cmap('viridis')
        values = np.linspace(0, 1, 10)
        label_colors = [colormap(value) for value in values]
        label_colors[5] = (1, 0, 0, 1)
        print(label_colors)
        for buffi, buff in enumerate(['C', 'CD', 'CG']):
            buff_df = plot_dfs[plot_dfs['buff'] == buff]
            print(buff_df)

            for timei, time in enumerate(['R1min', 'R5min', 'R2hr', 'Rall']):
                time_df = buff_df[buff_df['timepoint'] == time]
                print(time_df)
                if buffi == 0:
                    rownames += [time]
                for labeli, (label, label_df) in enumerate(time_df.groupby('cov')):
                    #print(timei, buffi)
                    x = label_df['SPA'].values

                    ## Plot the OR
                    y = label_df['CDF'].values

                    axes[timei, buffi].plot(x, y, label=label, color=label_colors[labeli])
                    axes[timei, buffi].set_xscale('log')
                    axes[timei, buffi].set_xlabel('LiP-MS coverage, (% primary structure)')
                    axes[timei, buffi].set_ylabel('Cumulative probability (CDF)')
                    axes[timei, buffi].set_title(f'{buff_tag[buff]}, {time}')

        axes[0,2].legend(bbox_to_anchor=(1.1, 1), title="Coverage", frameon=False)
        plt.tight_layout(pad=3)

        ######################################
        ## adjust axes locations and labels
        panel_labels = np.asarray([['a', 'b', 'c'], ['d', 'e', 'f'], ['g', 'h', 'i'], ['j', 'k', 'l']])
        #panel_label_ys = {0:1, 1:2/3, 2:1/3}
        #x0s = {0:0.075, 1:0.365, 2:2/3}
        #y0s = {0:0.7, 1:0.4, 2:0.1}
        panel_label_ys = {0:0.975, 1:0.74, 2:0.5, 3:0.26}
        x0s = {0:0.075, 1:0.415, 2:0.735}
        y0s = {0:0.83, 1:0.56, 2:0.32, 3:0.07}
        for coli, colname in enumerate(['cyto-serum', '+DnaK', '+GroEL']):
            for rowi, rowname in enumerate(rownames):

                # Adjust the linewidth of the axis spines
                for spine in axes[rowi, coli].spines.values():
                    spine.set_linewidth(0.5)  # Set the linewidth for all spines

                # Adjust the linewidth of the ticks
                axes[rowi, coli].tick_params(width=0.5)  # Both major and minor ticks

                axes[rowi, coli].set_ylabel('Cumulative probability')
                axes[rowi, coli].set_xlabel('Sum of Peptide Abundances, (SPA)')

                axes[rowi, coli].tick_params(axis='y', labelsize=6)
                axes[rowi, coli].tick_params(axis='x', labelsize=6)

                # Remove the right and top spines
                axes[rowi, coli].spines['right'].set_visible(False)
                axes[rowi, coli].spines['top'].set_visible(False)

                #if rowi == 1:
                #    axes[rowi, coli].set_ylim(top=1)

                axs_position = axes[rowi, coli].get_position()
                #x0 = x0s[coli]
                #y0 = y0s[rowi]
                #width0, height0 = axs_position.extents[2] - axs_position.extents[0], axs_position.extents[3] - axs_position.extents[1]
                #print(axs_position, width0, height0)
                #axes[rowi, coli].set_position([x0, y0, width0, 0.15])  # [left, bottom, width, height]

                bbox_in_fig_coords = axes[rowi, coli].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
                fig.text(bbox_in_fig_coords.x0, panel_label_ys[rowi], panel_labels[rowi, coli], fontsize=8, fontweight='bold', va='top', ha='left')
        #######################################

        #axes[0, 3].set_axis_off()
        #axes[1, 3].set_axis_off()
        #axes[2, 3].set_axis_off()

        #######################################
        figure_outpath = os.path.join(self.out_path, 'MainTextFigure1_Supp.pdf')
        plt.savefig(figure_outpath)
        print(f'SAVED: {figure_outpath}')

        figure_outpath = os.path.join(self.out_path, 'MainTextFigure1_Supp.png')
        plt.savefig(figure_outpath)
        print(f'SAVED: {figure_outpath}')

        figure_outpath = os.path.join(self.out_path, 'MainTextFigure1_Supp.svg')
        plt.savefig(figure_outpath)
        print(f'SAVED: {figure_outpath}')
        #######################################  
        quit()      
    #################################################################################################################

    #################################################################################################################
    def plot_overlap_Supp(self,):
        """
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Processing_LiP-MS_data/Overlap/EXP/Data/GeneOverlapTable_spa50_LiPMScov50_ent_genes_C.csv
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Processing_LiP-MS_data/Overlap/EXP/Data/PKsiteOverlap2Table_spa50_LiPMScov50_ent_genes_C.csv
        """
        buff_tag = {'C': 'cyto-serum', 'CD': '+DnaK', 'CG': '+GroEL'}
        #######################################
        dfs = []
        for dataset in ['ent_genes', 'all_genes']:
            for num in [2, 3]:
                for typeoverlap in ['GeneOverlapTableVenn', 'PKsiteOverlap']:
                    for buff in ['C', 'CD', 'CG']:
                        inp = glob.glob(f'{self.slug_path}/Processing_LiP-MS_data/Overlap/EXP/Data/{typeoverlap}{num}*_spa50_LiPMScov50_ent_genes_{buff}.csv')
                        print(f'inp: {inp}')
                        df = pd.read_csv(inp[0])
                        df['num'] = num
                        df['typeoverlap'] = typeoverlap
                        df['dataset'] = dataset
                        if typeoverlap == 'GeneOverlapTableVenn':
                            df = df.rename(columns={'count':'Count'})
        
                        dfs += [df]
        Figure_overlap_df = pd.concat(dfs)
        print(Figure_overlap_df)

        # Save the figure 1a raw plot df
        Figure_overlap_outfile = os.path.join(self.out_path, f'cutsiteOverlap_Supp.csv')
        Figure_overlap_df.to_csv(Figure_overlap_outfile)
        print(f'SAVED: {Figure_overlap_outfile}')
        ######################################

        for dataset in ['ent_genes', 'all_genes']:
            for typeoverlap in ['GeneOverlapTableVenn', 'PKsiteOverlap']:
                typeoverlap_df = Figure_overlap_df[(Figure_overlap_df['dataset'] == dataset) & (Figure_overlap_df['typeoverlap'] == typeoverlap)]
                #######################################
                ## set up figure(s)
                fig_width_mm = 180  # Width in mm
                fig_height_mm = 170  # Height in mm
                fig, axes = plt.subplots(3, 2, figsize=(mm_to_inches(fig_width_mm), mm_to_inches(fig_height_mm)))
                #######################################

                #######################################
                ## plot data
                for buffi, buff in enumerate(['C', 'CD', 'CG']):
                    buff_df = typeoverlap_df[typeoverlap_df['buff'] == buff]
                    print(buff_df)

                    for num_i, num in enumerate([3, 2]):
                        num_df = buff_df[buff_df['num'] == num]
                        print(num_df)
                        plot_venn(num_df, ax=axes[buffi, num_i], title=f'{buff_tag[buff]}')
                     
                plt.suptitle(f'{dataset}, {typeoverlap}')
                #######################################
                figure_outpath = os.path.join(self.out_path, f'cutsiteOverlap_{dataset}_{typeoverlap}_Supp.pdf')
                plt.savefig(figure_outpath)
                print(f'SAVED: {figure_outpath}')

                figure_outpath = os.path.join(self.out_path, f'cutsiteOverlap_{dataset}_{typeoverlap}_Supp.png')
                plt.savefig(figure_outpath)
                print(f'SAVED: {figure_outpath}')

                figure_outpath = os.path.join(self.out_path, f'cutsiteOverlap_{dataset}_{typeoverlap}_Supp.svg')
                plt.savefig(figure_outpath)
                print(f'SAVED: {figure_outpath}')
                #######################################  
                plt.close()
    #################################################################################################################

#################################################################################################################
def plot_venn(dataframe, ax=None, title=''):
    """
    Creates a Venn diagram on the provided axis object using the input dataframe.

    Parameters:
        dataframe (pd.DataFrame): Input dataframe containing 'OverlapClass' and 'count' columns.
        ax (matplotlib.axes.Axes): Axis object to plot the Venn diagram. If None, creates a new axis.

    Returns:
        matplotlib.axes.Axes: The axis with the Venn diagram plotted.
    """
    # Extract counts and percentages
    counts = dataframe.set_index('OverlapClass')['Count']
    percentages = dataframe.set_index('OverlapClass')['Overlap'] * 100  # Convert fractions to percentages

    # Determine the type of Venn diagram based on the unique OverlapClass labels
    overlap_classes = counts.index.tolist()

    # Create an axis if one is not provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))

    if len(dataframe) == 7:  # 3D Venn Diagram
        subsets = {
            '100': counts.get('R1min', 0),
            '010': counts.get('R5min', 0),
            '001': counts.get('R2hr', 0),
            '110': counts.get('R1min-R5min', 0),
            '101': counts.get('R1min-R2hr', 0),
            '011': counts.get('R5min-R2hr', 0),
            '111': counts.get('R1min-R5min-R2hr', 0),
        }
        venn = venn3(subsets=subsets, set_labels=('R1min', 'R5min', 'R2hr'), ax=ax)

    elif len(overlap_classes) == 3:  # 2D Venn Diagram
        subsets = {
            '10': counts.get('R5min', 0),
            '01': counts.get('R2hr', 0),
            '11': counts.get('R5min-R2hr', 0),
        }
        venn = venn2(subsets=subsets, set_labels=('R5min', 'R2hr'), ax=ax)

    # Add title to the axis
    ax.set_title(title)

    return ax
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


    ## make subplot figures 
    plotter.plot_LiPMS_Supp()
    plotter.plot_overlap_Supp()
    quit()

    # final formating and output
    # Automatically adjust layout
    #fig.tight_layout(pad=2.0)  # 'pad' is the overall padding
    figure_outpath = os.path.join(args.out_path, 'Figure2.1.pdf')
    plt.savefig(figure_outpath)
    print(f'SAVED: {figure_outpath}')

    print('NORMAL TERMINATION')

if __name__ == "__main__":
    main()

