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
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg

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
    Figure Xa (row 1 column 1): 
        MM kinetics of PK
        data/MM_kinetics/MM_Native_agg_processed.csv
        data/MM_kinetics/MM_Refolded_agg_processed.csv

    Figure Xb (row 1 column 2):
        SEC data for PK
        data/SEC/PK_E*
        
    Figure Xc: 
        Native page image
              data/NativePAGE_raw.png

    Figure Xd (row 2 column 2):
        ANS results
        data/ANS/

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
    def plot_Figure_Xa(self, ax):
        """
        Figure Xa (row 1 column 1): 
            MM kinetics of PK
            data/MM_kinetics/MM_Native_agg_processed.csv
            data/MM_kinetics/MM_Refolded_agg_processed.csv
        """
        buff_tag = {'C': 'cyto-serum', 'CD': '+DnaK', 'CG': '+GroEL'}
        #######################################
        ## Load Figure 1a data
        inp = f'data/MM_kinetics/MM_Native_agg_processed.csv'
        print(f'inp: {inp}')
        Figure_Xa_Native_df = pd.read_csv(inp)
        Figure_Xa_Native_df.rename(columns={'N_avg': 'avg', 'N_ci_lower_delta':'lb', 'N_ci_upper_delta':'ub'}, inplace=True)
        print(f'Figure_Xa_Native_df:\n{Figure_Xa_Native_df}')

        inp = f'data/MM_kinetics/MM_Refolded_agg_processed.csv'
        print(f'inp: {inp}')
        Figure_Xa_Refolded_df = pd.read_csv(inp)
        Figure_Xa_Refolded_df.rename(columns={'R_avg': 'avg', 'R_ci_lower_delta':'lb', 'R_ci_upper_delta':'ub'}, inplace=True)
        print(f'Figure_Xa_Refolded_df:\n{Figure_Xa_Refolded_df}')
        #######################################

        # Define the custom function
        def MM_fit(X, a, b, n):
            return (a * X**n) / ((b**n) + X**n)
        

        # Adjust the linewidth of the axis spines
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)  # Set the linewidth for all spines

        # Adjust the linewidth of the ticks
        ax.tick_params(width=0.5)  # Both major and minor ticks

        for label, df, color, popt in [('Untreated', Figure_Xa_Native_df, 'black', [0.0058317, 7.1, 4.6]), ('Treated', Figure_Xa_Refolded_df, 'red', [0.0042221, 7.57, 4.42])]:

            # Fit the data
            X = df['pep'].values[:7]
            Y = df['avg'].values[:7]
            #popt, pcov = curve_fit(MM_fit, X, Y, p0=[1, 1, 1])  # Provide initial guesses for a, b, and n
            print(f'fit coef: {popt}')

            # Generate fitted values
            X_fit = np.linspace(1, 20, 500)  # Generate smooth X values for the fit curve

            #Y_fit = MM_fit(df['pep'].values, *popt)
            Y_fit = MM_fit(X_fit, *popt)

            # Calculate error bars
            yerr = [df['lb'], df['ub']]
            
            # Plot the trace with error bars
            ax.errorbar(
                df['pep'], 
                df['avg'], 
                yerr=yerr, 
                fmt='o',  # Line and point markers
                capsize=3,  # Add caps to the error bars
                color=color,
                markersize=3,
                markerfacecolor=color,  
                elinewidth=0.5,  
                markeredgewidth=0.5)

            ax.plot(X_fit, Y_fit, label=label, color=color, linewidth=0.75)  # Plot the fit curve

        # Add a dashed black line at y=1    
        ax.set_ylabel('Reaction rate (μmols NADH $\cdot$s$^{-1}$)')
        ax.set_xlabel('[PEP] (mM)')

        # Customize the plot
        ax.set_ylim(0.0, 0.007)
        ax.set_xlim(-0.5, 20.5)

        ax.set_xticks([0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20], labels=[0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20])  # Rotate x-axis labels for better readability

        ax.tick_params(axis='y', labelsize=6)
        ax.tick_params(axis='x', labelsize=6, rotation=45)
        # Remove the right and top spines

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        #ax.legend(loc='upper left')
        ax.legend(fontsize='small', handlelength=1, handletextpad=0.5, borderpad=0.3, frameon=False)

        return ax
    #################################################################################################################

    #################################################################################################################
    def plot_Figure_Xb(self, axN, axR):
        """
        Figure Xb (row 1 column 2):
            SEC data for PK
            data/SEC/PK_E*
        """

        #######################################
        ## Load data
        df = pd.read_csv('data/SEC/combined_SEC.csv')
        print(df)
        X_keys = [key for key in df.keys() if '_ml' in key]
        Y_keys = [key for key in df.keys() if '_mAU' in key]
        print(f'X_keys: {X_keys}')
        print(f'Y_keys: {Y_keys}')

        ## get offset to correct y values so they are all on the same baseline
        # use values for 0.25 < x < 0.5 from the first trace to get the baseline
        for xkey, ykey in zip(X_keys, Y_keys):
            baseline = df[[xkey, ykey]]
            baseline = baseline[(baseline[xkey] >= 0.25) & (baseline[xkey] <= 0.5)]
            #print(baseline)
            baseline = np.mean(baseline[ykey])
            print(xkey, ykey, baseline)
            df[ykey] -= baseline - 0.05
        print(df)

        # Define one set of common bins across both x columns
        all_x_values = np.concatenate([df[key].dropna() for key in X_keys])
        bins = np.linspace(all_x_values.min(), all_x_values.max(), num=2000)  # Adjust number of bins as needed
        print(bins)

        # Use np.digitize to assign each x value to a bin
        x = [] # x values (bin centers)
        Ny = [] # average mAU in the bin for the Native samples
        Nci = [] # ci 95%
        Ry = [] # average mAU in the bin for the Refolded samples
        Rci = [] # ci 95%
        for i, _ in enumerate(bins):
            if i == len(bins) - 1:
                break
            start, end, center = bins[i], bins[i+1], bins[i] + (bins[i+1] - bins[i])/2
            #print(start, end, center)

            Nyvals = []
            Ryvals = []
            for xkey, ykey in zip(X_keys, Y_keys):
                #print(xkey, ykey)
                loc_df = df[[xkey, ykey]]
                loc_df = loc_df[(loc_df[xkey] >= start) & (loc_df[xkey] < end)]
                #print(loc_df)

                if len(loc_df) == 0:
                    continue
                
                if '_N_' in xkey:
                    Nyvals += [loc_df[ykey].values]
                if '_R_' in xkey:
                    Ryvals += [loc_df[ykey].values]

            # average the native sample y values and get the confidence interval
            if len(Nyvals) != 0:
                Nyvals = np.hstack(Nyvals)
                Nyavg = np.mean(Nyvals)
                Nystd = np.std(Nyvals)
                Nlen = len(Nyvals)
                Nyci = (Nystd/np.sqrt(Nlen))*1.96
            else:
                Nyavg = 0
                Nyci = 0

            if len(Ryvals) != 0:
                Ryvals = np.hstack(Ryvals)
                Ryavg = np.mean(Ryvals)
                Rystd = np.std(Ryvals)
                Rlen = len(Ryvals)
                Ryci = (Rystd/np.sqrt(Rlen))*1.96
            else:
                Ryavg = 0
                Ryci = 0

            x += [center] # x values (bin centers)
            Ny += [Nyavg] # average mAU in the bin for the Native samples
            Nci += [Nyci] # ci 95%
            Ry += [Ryavg] # average mAU in the bin for the Refolded samples
            Rci += [Ryci] # ci 95%

        plot_df = pd.DataFrame({'x':x, 'Ny':Ny, 'Nci':Nci, 'Ry':Ry, 'Rci':Rci})
        plot_df = plot_df[(plot_df['x'] >= 0) & (plot_df['x']<=2.5)]
        print(plot_df)
      
        # Save the figure 1a raw plot df
        Figure_Xb_outfile_csv = os.path.join(self.out_path, f'Figure_Xb.csv')
        plot_df.to_csv(Figure_Xb_outfile_csv)
        print(f'SAVED: {Figure_Xb_outfile_csv}')
        #######################################
       
        # Adjust the linewidth of the axis spines
        for spine in axN.spines.values():
            spine.set_linewidth(0.5)  # Set the linewidth for all spines
        # Adjust the linewidth of the axis spines
        for spine in axR.spines.values():
            spine.set_linewidth(0.5)  # Set the linewidth for all spines

        # Adjust the linewidth of the ticks
        axN.tick_params(width=0.5)  # Both major and minor ticks
        axR.tick_params(width=0.5)  # Both major and minor ticks

        # Plot the trace with error bar
        axN.plot(plot_df['x'], plot_df['Ny'], label='Untreated', ls='-', linewidth=0.75, color='black')
        axN.fill_between(plot_df['x'], 
                 plot_df['Ny'] - plot_df['Nci'], 
                 plot_df['Ny'] + plot_df['Nci'], 
                 alpha=0.2, label='Untreated', color='black')

        axR.plot(plot_df['x'], plot_df['Ry'], label='Treated', ls='-', linewidth=0.75, color='red')
        axR.fill_between(plot_df['x'], 
                 plot_df['Ry'] - plot_df['Rci'], 
                 plot_df['Ry'] + plot_df['Rci'], 
                 alpha=0.2, label='Treated', color='red')
            
        
        # Add a dashed black line at y=1    
        axN.set_ylabel('UV abs @ 280nm')
        axR.set_ylabel('UV abs @ 280nm')
        axN.set_xlabel('Volume since injection (mL)')
        axR.set_xlabel('Volume since injection (mL)')
            
        # Customize the plot
        #print(f'Max(y) ceil: {maxy} {np.ceil(maxy)}')
        axN.set_ylim(0.0, 0.6)
        axR.set_ylim(0.0, 0.6)
        axN.set_xlim(0, 2.5)
        axR.set_xlim(0, 2.5)

        axN.tick_params(axis='y', labelsize=6)
        axN.tick_params(axis='x', labelsize=6)
        axR.tick_params(axis='y', labelsize=6)
        axR.tick_params(axis='x', labelsize=6)

        # Remove the right and top spines
        axN.spines['right'].set_visible(False)
        axN.spines['top'].set_visible(False)
        axR.spines['right'].set_visible(False)
        axR.spines['top'].set_visible(False)

        #axN.legend(frameon=False, loc='upper left')
        #axR.legend(frameon=False, loc='upper left')

        axN.text(0.05, 1, "Untreated", transform=axN.transAxes, va='top', ha='left', fontsize=6) 
        axR.text(0.05, 1, "Treated", transform=axR.transAxes, va='top', ha='left', fontsize=6) 
    #################################################################################################################

    #################################################################################################################
    def plot_Figure_Xc(self, ax):
        """
        Figure Xc: 
            SEC image
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Chaperone_Client_Associations/Dnak_simple_motif_scan/Plots/consolidated_Dnak_binding_data_EXP.csv     
        """
        # Load the image from a file
        image_path = 'data/NativePAGE_raw.png'  # Replace with your image file path
        image = mpimg.imread(image_path)

        # Display the image in the axis
        ax.imshow(image, extent=[-0.5, 661.5, 638.5, -0.5])
        im = ax.imshow(image)
        print(im.get_extent())
        # Set the axis limits to center the image
        #ax.set_xlim(-2, 2)  # Adjust these limits to center and scale
        #ax.set_ylim(-2, 2)

        # Optionally remove axis ticks/labels for a cleaner look
        ax.axis('off')

        ax.text(1, 0.475, "Tetramer", transform=ax.transAxes, va='top', ha='left', fontsize=6) 
        ax.text(-0.01, 1.075, "kDa", transform=ax.transAxes, va='top', ha='right', fontsize=6)
        ax.text(-0.01, 0.985, "1236", transform=ax.transAxes, va='top', ha='right', fontsize=6) 
        ax.text(-0.01, 0.925, "1048", transform=ax.transAxes, va='top', ha='right', fontsize=6) 
        ax.text(-0.01, 0.705, "720", transform=ax.transAxes, va='top', ha='right', fontsize=6) 
        ax.text(-0.01, 0.575, "480", transform=ax.transAxes, va='top', ha='right', fontsize=6) 
        ax.text(-0.01, 0.425, "242", transform=ax.transAxes, va='top', ha='right', fontsize=6) 
        ax.text(-0.01, 0.275, "146", transform=ax.transAxes, va='top', ha='right', fontsize=6) 
        ax.text(-0.01, 0.175, "66", transform=ax.transAxes, va='top', ha='right', fontsize=6) 
        ax.text(-0.01, 0.015, "20", transform=ax.transAxes, va='bottom', ha='right', fontsize=6) 

        ax.text(0.1, 0, "ladder", transform=ax.transAxes, va='top', ha='right', fontsize=6, rotation=45)
        ax.text(0.25, 0, "Untreated", transform=ax.transAxes, va='top', ha='right', fontsize=6, rotation=45)
        ax.text(0.4, 0, "Treated", transform=ax.transAxes, va='top', ha='right', fontsize=6, rotation=45)
        ax.text(0.7, 0, "Misfolded\nPK peak", transform=ax.transAxes, va='top', ha='right', fontsize=6, rotation=45)
        ax.text(1, 0, "Folded\nPK peak", transform=ax.transAxes, va='top', ha='right', fontsize=6, rotation=45)
    #################################################################################################################

    #################################################################################################################
    def plot_Figure_Xd(self, ax):
        """
        Figure Xd (row 2 column 2):
            ANS results
            data/ANS/
        """

        #######################################
        ## Load 
        ANS_files = glob.glob(f'data/ANS/*')
        x = []
        Ny = []
        Ry = []
        for ANS_f in ANS_files:
            #print(ANS_f)
            df = pd.read_csv(ANS_f, sep=' ')
            #print(df)

            Xkey = [key for key in df.keys() if '_wavelength' in key][0]
            Ykey = [key for key in df.keys() if '_emi_intensity' in key][0]
            #print(Xkey, Ykey)
            if 'N1_' in Xkey:
                x = df[Xkey].values

            if 'N' in Ykey:
                Ny += [df[Ykey].values]
            elif 'R' in Ykey:
                Ry += [df[Ykey].values]

        x = np.asarray(x)
        #print(f'x: {x} {x.shape}')
        Ny = np.stack(Ny)
        #print(f'Ny: {Ny} {Ny.shape}')
        Ny_avg = np.mean(Ny, axis=0)
        #print(f'Ny_avg: {Ny_avg} {Ny_avg.shape}')
        Ny_std = np.std(Ny, axis=0)
        Ny_ci = (Ny_std/np.sqrt(8))*1.96
        #print(f'Ny_ci: {Ny_ci} {Ny_ci.shape}')

        Ry = np.stack(Ry)
        #print(f'Ry: {Ry} {Ry.shape}')
        Ry_avg = np.mean(Ry, axis=0)
        #print(f'Ry_avg: {Ry_avg} {Ry_avg.shape}')
        Ry_std = np.std(Ry, axis=0)
        Ry_ci = (Ry_std/np.sqrt(8))*1.96
        #print(f'Ry_ci: {Ry_ci} {Ry_ci.shape}')

        plot_df = pd.DataFrame({'x':x, 'Ny_avg':Ny_avg, 'Ny_ci':Ny_ci, 'Ry_avg':Ry_avg, 'Ry_ci':Ry_ci})
        print(plot_df)

        # Save the figure 1a raw plot df
        Figure_Xd_outfile_csv = os.path.join(self.out_path, f'Figure_Xd.csv')
        plot_df.to_csv(Figure_Xd_outfile_csv)
        print(f'SAVED: {Figure_Xd_outfile_csv}')
        #######################################

        # Adjust the linewidth of the axis spines
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)  # Set the linewidth for all spines

        # Adjust the linewidth of the ticks
        ax.tick_params(width=0.5)  # Both major and minor ticks

        for label, ykey, cikey, color in [('Untreated', 'Ny_avg', 'Ny_ci', 'black'), ('Treated', 'Ry_avg', 'Ry_ci', 'red')]:
            # Plot the trace with error bars
            ax.errorbar(
                plot_df['x'], 
                plot_df[ykey], 
                yerr=plot_df[cikey], 
                fmt='o',  # Line and point markers
                capsize=3,  # Add caps to the error bars
                color=color,
                markersize=3,
                markerfacecolor=color,  
                elinewidth=0.5,  
                markeredgewidth=0.5,
                label=label)

        #ax.set_ylim(0.5,np.ceil(maxy))
        #ax.set_xlim(-0.5, 2.5)

        ax.tick_params(axis='y', labelsize=6)
        ax.tick_params(axis='x', labelsize=6)

        # Remove the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        ax.legend(fontsize='small', handlelength=1, handletextpad=0.5, borderpad=0.3, frameon=False)
        ax.set_ylabel('Fluorescence intensity')
        ax.set_xlabel('Emission Wavelength (nm)')

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
    fig_height_mm = 60  # Height in mm
    #custom_font_path = "/storage/group/epo2/default/ims86/miniconda3/envs/FtoF/fonts/Arial.ttf" # Path to your custom font
    #arial_font = fm.FontProperties(fname=custom_font_path) # Create a FontProperties object
    plt.rcParams['font.family'] = 'Arial'  # Change to your desired font, e.g., 'Times New Roman', 'DejaVu Sans', etc.
    plt.rcParams['font.size'] = 6  # Default font size
    plt.rcParams['pdf.fonttype'] = 42
    #fig, axs = plt.subplots(1, 5, figsize=(mm_to_inches(fig_width_mm), mm_to_inches(fig_height_mm)), dpi=600, constrained_layout=False, gridspec_kw={'width_ratios': [1, 1, 1, 3]})  
    
    fig = plt.figure(figsize=(mm_to_inches(fig_width_mm), mm_to_inches(fig_height_mm)), dpi=600)
    gs = gridspec.GridSpec(2, 4, figure=fig)
    # Add subplots
    ax1 = fig.add_subplot(gs[:, 0])  # First row, first column
    ax2 = fig.add_subplot(gs[0, 1])  # First row, second column
    ax3 = fig.add_subplot(gs[1, 1])  # First row, second column
    ax4 = fig.add_subplot(gs[:, 2])  # Second row, first column
    ax5 = fig.add_subplot(gs[:, 3])  # Second row, second column
    axs = [ax1, ax2, ax3, ax4, ax5]

    plotter.plot_Figure_Xa(axs[0])
    plotter.plot_Figure_Xb(axs[1], axs[2])
    plotter.plot_Figure_Xc(axs[3])
    plotter.plot_Figure_Xd(axs[4])

    #########################################################
    # Get the current position of the subplots and adjust their positions manually
    ## MM kinetics
    axs0_position = axs[0].get_position()
    width0, height0 = axs0_position.extents[2] - axs0_position.extents[0], axs0_position.extents[3] - axs0_position.extents[1]
    print(axs0_position, width0, height0)
    axs[0].set_position([0.07, 0.17, width0, height0])  # [left, bottom, width, height]

    bbox_in_fig_coords = axs[0].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bbox_in_fig_coords.x0, 0.995, 'a', fontsize=8, fontweight='bold', va='top', ha='left')


    ## SEC
    axs1_position = axs[1].get_position()
    width1, height1 = axs1_position.extents[2] - axs1_position.extents[0], axs1_position.extents[3] - axs1_position.extents[1]
    print(axs1_position, width1, height1)
    axs[1].set_position([0.31, axs1_position.y0 + 0.1, 0.12, 0.3])

    bbox_in_fig_coords = axs[1].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bbox_in_fig_coords.x0, 0.995, 'b', fontsize=8, fontweight='bold', va='top', ha='left')

    axs2_position = axs[2].get_position()
    width2, height2 = axs2_position.extents[2] - axs2_position.extents[0], axs2_position.extents[3] - axs2_position.extents[1]
    print(axs2_position, width2, height2)
    axs[2].set_position([0.31, 0.1575, 0.12, 0.3])


    ## Native page
    axs3_position = axs[3].get_position()
    width3, height3 = axs3_position.extents[2] - axs3_position.extents[0], axs3_position.extents[3] - axs3_position.extents[1]
    print(axs3_position, width3, height3)
    axs[3].set_position([0.5, 0.1575, 0.2, 0.77])

    bbox_in_fig_coords = axs[3].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bbox_in_fig_coords.x0, 0.995, 'c', fontsize=8, fontweight='bold', va='top', ha='left')


    ## ANS
    axs4_position = axs[4].get_position()
    width4, height4 = axs4_position.extents[2] - axs4_position.extents[0], axs4_position.extents[3] - axs4_position.extents[1]
    print(axs4_position, width4, height4)
    axs[4].set_position([0.825, 0.1575, width4, height4])

    bbox_in_fig_coords = axs[4].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bbox_in_fig_coords.x0, 0.995, 'd', fontsize=8, fontweight='bold', va='top', ha='left')


    #########################################################

    # final formating and output
    # Automatically adjust layout
    #fig.tight_layout()  # 'pad' is the overall padding
    figure_outpath = os.path.join(args.out_path, 'Figure_experiments.pdf')
    plt.savefig(figure_outpath)
    print(f'SAVED: {figure_outpath}')

    figure_outpath = os.path.join(args.out_path, 'Figure_experiments.svg')
    plt.savefig(figure_outpath)
    print(f'SAVED: {figure_outpath}')

    print('NORMAL TERMINATION')

if __name__ == "__main__":
    main()

