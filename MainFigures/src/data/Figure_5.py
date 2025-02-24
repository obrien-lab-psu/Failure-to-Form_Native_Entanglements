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

    Figure is 2 row by 2 column 
    Figure 5a (row 1 column 1): 
        Association between essentiality and being a client of DnaK or GroEL
        C buffers at spa50 and LiPMScov 50 in EXP dataset
        make a contingency table with the values and put the OR and pvalue below
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Chaperone_Client_Associations/Assoc_Client_n_Essential/Plots/Assoc_Client_n_Essential_plot_data_Knockout_clientType-DnaKonly.csv
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Chaperone_Client_Associations/Assoc_Client_n_Essential/Plots/Assoc_Client_n_Essential_plot_data_Knockout_clientType-GroELonly.csv

    Figure 5b (row 1 column 2):
        Entanglement complexity as a discrimination tool
        C, CD, CG buffers at spa50 and LiPMScov 50 and has to have an entanglement
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Entanglement_Topological_Complexity_and_Discrimination/EXP/Lasso_plot_results_C_50.csv
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Entanglement_Topological_Complexity_and_Discrimination/EXP/Lasso_plot_results_CD_50.csv
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Entanglement_Topological_Complexity_and_Discrimination/EXP/Lasso_plot_results_CG_50.csv
        
    Figure 5c: 
        DnaK binding motif analysis
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Chaperone_Client_Associations/Dnak_simple_motif_scan/Plots/consolidated_Dnak_binding_data_EXP.csv       

    Figure 5d (row 2 column 2):
        Loop closing contacts enrichement matrix
        C buffers at spa50 and LiPMScov 50 in EXP dataset
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Sequence_Complexity_and_Discrimination/Compare_OR_with_permutation/EXP/C_50_p100000/FrequencyGeneratorOutput/OR_FDR_pvalues.csv
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Sequence_Complexity_and_Discrimination/Compare_OR_with_permutation/EXP/C_50_p100000/FrequencyGeneratorOutput/OR_GT.csv

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
    def plot_Figure_5a(self, ax):
        """
        Figure 5a (row 1 column 1): 
            Association between essentiality and being a client of DnaK or GroEL
            C buffers at spa50 and LiPMScov 50 in EXP dataset
            make a contingency table with the values and put the OR and pvalue below
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Chaperone_Client_Associations/Assoc_Client_n_Essential/Plots/Assoc_Client_n_Essential_plot_data_Knockout_clientType-DnaKonly.csv
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Chaperone_Client_Associations/Assoc_Client_n_Essential/Plots/Assoc_Client_n_Essential_plot_data_Knockout_clientType-GroELonly.csv
        """
        buff_tag = {'C': 'cyto-serum', 'CD': '+DnaK', 'CG': '+GroEL'}
        #######################################
        ## Load Figure 1a data
        inp = f'{self.slug_path}/Chaperone_Client_Associations/Assoc_Client_n_Essential/Plots/Assoc_Client_n_Essential_plot_data_Knockout_clientType-DnaKonly.csv'
        print(f'inp: {inp}')
        Figure_5a_DnaK_df = pd.read_csv(inp)
        Figure_5a_DnaK_df = Figure_5a_DnaK_df[(Figure_5a_DnaK_df['buff'] == 'C') & (Figure_5a_DnaK_df['spa'] == 50) & (Figure_5a_DnaK_df['LiPMScov'] == 50) & (Figure_5a_DnaK_df['label'] == 'EXP')]
        Figure_5a_DnaK_df['Chap'] = '+DnaK'
        #print(f'Figure_5a_DnaK_df:\n{Figure_5a_DnaK_df}')

        inp = f'{self.slug_path}/Chaperone_Client_Associations/Assoc_Client_n_Essential/Plots/Assoc_Client_n_Essential_plot_data_Knockout_clientType-GroELonly.csv'
        print(f'inp: {inp}')
        Figure_5a_GroEL_df = pd.read_csv(inp)
        Figure_5a_GroEL_df = Figure_5a_GroEL_df[(Figure_5a_GroEL_df['buff'] == 'C') & (Figure_5a_GroEL_df['spa'] == 50) & (Figure_5a_GroEL_df['LiPMScov'] == 50) & (Figure_5a_DnaK_df['label'] == 'EXP')]
        Figure_5a_GroEL_df['Chap'] = '+GroEL'
        #print(f'Figure_5a_GroEL_df:\n{Figure_5a_GroEL_df}')
        #######################################

        # Save the figure 1a raw plot df
        Figure_5a_df = pd.concat([Figure_5a_DnaK_df, Figure_5a_GroEL_df])
        print(Figure_5a_df)
        Figure_5a_outfile_csv = os.path.join(self.out_path, f'Figure_5a.csv')
        Figure_5a_df.to_csv(Figure_5a_outfile_csv)
        print(f'SAVED: {Figure_5a_outfile_csv}')


        for df, label, bbox, textxy in [(Figure_5a_DnaK_df, 'DnaK', [0.2, 0.675, 0.6, 0.3], [0.5, 0.49]), (Figure_5a_GroEL_df, 'GroEL', [0.2, 0.0575, 0.6, 0.3], [0.5, -0.125])]:
            # Extract the required values from the DataFrame
            values = np.asarray([
                [df['clientY_essentialY'].values[0], df['clientY_essentialN'].values[0]],
                [df['clientN_essentialY'].values[0], df['clientN_essentialN'].values[0]]])
            print(values)

            # make odds ratio and pvalue strings
            OR = df['OR'].values[0]
            OR_lb, OR_ub = df['OR_lb'].values[0], df['OR_ub'].values[0]
            OR_str = f'Odds Ratio = {OR:.2f} ({OR_lb:.2f}, {OR_ub:.2f})'
            print(OR_str)

            pvalue = df['pvalue'].values[0]
            pvalue_str = f'p-value = {format_scientific(pvalue)}'
            print(pvalue_str)
        
            # Define row and column labels
            row_labels = ['Yes ', 'No ']
            col_labels = ['Yes', 'No']

            # Hide the axes
            ax.axis('tight')
            ax.axis('off')

            # Create the table
            table = ax.table(
                cellText=values,
                rowLabels=row_labels,
                colLabels=col_labels,
                bbox=bbox,
                cellLoc='center')

            # Adjust the font size for better readability
            table.auto_set_font_size(False)
            table.set_fontsize(6)

            # Get the table bounding box in display coordinates
            table_bbox = table.get_window_extent()
            print(table_bbox)

            # Transform display coordinates to axes coordinates
            table_axes_coords = table_bbox.transformed(ax.transAxes.inverted())

            # Print the bounding box in axes coordinates
            print("Axes coordinates of the table:", table_axes_coords)

            # Remove borders around row and column label cells
            for key, cell in table.get_celld().items():
                print(key, cell, cell.get_text().get_text())
                if key[0] == 0 or key[1] == -1:  # Row or column label cells
                    cell.set_linewidth(0)  # Remove border
                    cell.set_edgecolor("none")  # Ensure no border color is shown
                else:
                    cell.set_linewidth(0.5)

            # Add rotated row header text to the left of the table
            ax.text(-0.15, bbox[1] + 0.1, f"{label} Client", transform=ax.transAxes, va='center', ha='left', rotation=90, fontsize=6, fontweight='bold')  

            # Add column header text to the top of the table
            ax.text(0.5, bbox[1] + 0.3, "Essential", transform=ax.transAxes, va='center', ha='center', fontsize=6, fontweight='bold') 
        
            # Add OR and pvalue strings below table
            ax.text(0.5, bbox[1] - 0.075, OR_str, transform=ax.transAxes, va='center', ha='center', fontsize=6)  
            ax.text(0.5, bbox[1] - 0.125, pvalue_str, transform=ax.transAxes, va='center', ha='center', fontsize=6) 
            
            # Add (insignificant) text
            ax.text(textxy[0], textxy[1], '(insignificant)', transform=ax.transAxes, va='center', ha='center', fontsize=6, fontweight='bold', color='#7E8D85')

        return ax
    #################################################################################################################

    #################################################################################################################
    def plot_Figure_5b(self, ax):
        """
        Figure 5b (row 1 column 2):
                Entanglement complexity as a discrimination tool
                C, CD, CG buffers at spa50 and LiPMScov 50 and has to have an entanglement
                ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Entanglement_Topological_Complexity_and_Discrimination/EXP/Lasso_plot_results_C_50.csv
                ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Entanglement_Topological_Complexity_and_Discrimination/EXP/Lasso_plot_results_CD_50.csv
                ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Entanglement_Topological_Complexity_and_Discrimination/EXP/Lasso_plot_results_CG_50.csv
        """
        buff_tag = {'C': 'cyto-serum', 'CD': '+DnaK', 'CG': '+GroEL'}
        #######################################
        ## Load Figure 1a data
        inp = f'{self.slug_path}/Entanglement_Topological_Complexity_and_Discrimination/EXP_022025/Lasso_plot_results_C_50.csv'
        print(f'inp: {inp}')
        Figure_5b_C_df = pd.read_csv(inp, sep='|')
        Figure_5b_C_df = Figure_5b_C_df[np.isclose(Figure_5b_C_df['C'], 1.5)]
        Figure_5b_C_df['buff'] = 'C'
        print(f'Figure_5b_C_df:\n{Figure_5b_C_df}')


        inp = f'{self.slug_path}/Entanglement_Topological_Complexity_and_Discrimination/EXP_022025/Lasso_plot_results_CD_50.csv'
        print(f'inp: {inp}')
        Figure_5b_CD_df = pd.read_csv(inp, sep='|')
        Figure_5b_CD_df = Figure_5b_CD_df[np.isclose(Figure_5b_CD_df['C'], 2.5)]
        Figure_5b_CD_df['buff'] = 'CD'
        print(f'Figure_5b_CD_df:\n{Figure_5b_CD_df}')

        inp = f'{self.slug_path}/Entanglement_Topological_Complexity_and_Discrimination/EXP_022025/Lasso_plot_results_CG_50.csv'
        print(f'inp: {inp}')
        Figure_5b_CG_df = pd.read_csv(inp, sep='|')
        Figure_5b_CG_df = Figure_5b_CG_df[np.isclose(Figure_5b_CG_df['C'], 2.5)]
        Figure_5b_CG_df['buff'] = 'CG'
        print(f'Figure_5b_CG_df:\n{Figure_5b_CG_df}')
     
        Figure_5b_df = pd.concat([Figure_5b_C_df, Figure_5b_CD_df, Figure_5b_CG_df])
        Figure_5b_df['xlabel'] = [buff_tag[b] for b in Figure_5b_df['buff'].values]
        print(Figure_5b_df)
      
        # Save the figure 1a raw plot df
        Figure_5b_outfile_csv = os.path.join(self.out_path, f'Figure_5b.csv')
        Figure_5b_df.to_csv(Figure_5b_outfile_csv)
        print(f'SAVED: {Figure_5b_outfile_csv}')
        #######################################
       
        # Adjust the linewidth of the axis spines
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)  # Set the linewidth for all spines

        # Adjust the linewidth of the ticks
        ax.tick_params(width=0.5)  # Both major and minor ticks

        maxy = 0
        # Calculate error bars
        yerr_lower = Figure_5b_df['<BA>'] - Figure_5b_df['BA_lb']
        yerr_upper = Figure_5b_df['BA_ub'] - Figure_5b_df['<BA>']
        yerr = [yerr_lower, yerr_upper]
        
        # Plot the trace with error bars
        ax.errorbar(
            Figure_5b_df['xlabel'], 
            Figure_5b_df['<BA>'], 
            yerr=yerr, 
            fmt='o',  # Line and point markers
            capsize=3,  # Add caps to the error bars
            color='black',
            markersize=3,
            markerfacecolor='black',  
            elinewidth=0.5,  
            markeredgewidth=0.5)

        ## mark pvalue annotations
        # Annotate significant p-values with '*'
        for x, y, err, p in zip(Figure_5b_df['xlabel'], Figure_5b_df['<BA>'], yerr_upper, Figure_5b_df['pvalue']):
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
        
   
        # Add a dashed black line at y=1    
        ax.axhline(y=0.5, color='#7E8D85', linestyle='--', linewidth=0.5)
        ax.set_ylabel('Balanced Accuracy')
        ax.set_xticklabels(Figure_5b_df['xlabel'], rotation=45, ha='right')  # Rotate x-axis labels for better readability
            
        # Customize the plot
        print(f'Max(y) ceil: {maxy} {np.ceil(maxy)}')
        ax.set_ylim(0.0, 1.0)
        ax.set_xlim(-0.5, 2.5)

        ax.tick_params(axis='y', labelsize=6)
        ax.tick_params(axis='x', labelsize=6)
        # Remove the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        # add the curly brackets
        ax.annotate("~7%", xy=(1.15, 0.55), xytext=(1.3, 0.55), xycoords='axes fraction', arrowprops=dict(arrowstyle='-[, widthB=1.0, lengthB=0.75', lw=0.5, color='#7E8D85'), va='center', ha='center', rotation=90, fontweight='bold', color='#7E8D85')
        #ax.annotate("(insignificant)", xy=(1.1, 0.2), xytext=(1.2, 0.2), xycoords='axes fraction', arrowprops=dict(arrowstyle='-[, widthB=4.0, lengthB=0.75', lw=0.5, color='k'), va='center', ha='center', rotation=90)

        return ax   
    #################################################################################################################

    #################################################################################################################
    def plot_Figure_5c(self, ax):
        """
        Figure 5c: 
            DnaK binding motif analysis
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Chaperone_Client_Associations/Dnak_simple_motif_scan/Plots/consolidated_Dnak_binding_data_EXP.csv       
        """
        buff_tag = {'C': 'cyto-serum', 'CD': '+DnaK', 'CG': '+GroEL'}
        #######################################
        ## Load Figure 2c data
        # Load post matched res_sasa
        inp = f'{self.slug_path}/Chaperone_Client_Associations/Dnak_simple_motif_scan/Plots/consolidated_Dnak_binding_data_EXP.csv'
        print(f'inp: {inp}')
        Figure_5c_df = pd.read_csv(inp)
        Figure_5c_df = Figure_5c_df[(Figure_5c_df['buff'] == 'C') & (Figure_5c_df['spa'] == 50) & (Figure_5c_df['OnlyEnt'] == True) & (Figure_5c_df['motif'] == 'Schymkowitz')]
        Figure_5c_df['xlabel'] = [ l.replace('D_', '') for l in Figure_5c_df['D_type'].values]
        print(Figure_5c_df)
        
        # Save the figure 1a raw plot df
        Figure_5c_outfile_csv = os.path.join(self.out_path, f'Figure_5c.csv')
        Figure_5c_df.to_csv(Figure_5c_outfile_csv)
        print(f'SAVED: {Figure_5c_outfile_csv}')
        ######################################
  
        # Adjust the linewidth of the axis spines
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)  # Set the linewidth for all spines

        # Adjust the linewidth of the ticks
        ax.tick_params(width=0.5)  # Both major and minor ticks
        maxy = 0
        for label, mean_label, lb_label, ub_label, color in [('Essential', 'Ess_mean', 'Ess_lower_ci', 'Ess_upper_ci', 'blue'), ('Non-Essential', 'NonEss_mean', 'NonEss_lower_ci', 'NonEss_upper_ci', 'red')]:
            print(label, mean_label, ub_label, lb_label)
            # Calculate error bars
            yerr_lower = Figure_5c_df[mean_label] - Figure_5c_df[lb_label]
            yerr_upper = Figure_5c_df[ub_label] - Figure_5c_df[mean_label]
            yerr = [yerr_lower.values * 100, yerr_upper.values * 100]
            print(f'{label} yerr: {yerr}')

            # Plot the trace with error bars
            ax.errorbar(
                Figure_5c_df['xlabel'], 
                Figure_5c_df[mean_label].values * 100, 
                yerr=yerr, 
                fmt='o',  # Line and point markers
                capsize=3,  # Add caps to the error bars
                color=color,
                markersize=3,
                markerfacecolor=color,  
                elinewidth=0.5,  
                markeredgewidth=0.5,
                label=label)

            ## mark pvalue annotations
            # Annotate significant p-values with '*'
            for x, y, err, p in zip(Figure_5c_df['xlabel'], Figure_5c_df[mean_label], yerr_upper, Figure_5c_df['q_value']):
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
                
            # Add a dashed black line at y=1    
            #ax.axhline(y=1, color='black', linestyle='--', linewidth=0.5)
            ax.set_ylabel('Number of DnaK binding motifs per 100 residues')
            #ax.set_ylabel(r'D ($\frac{\text{hits}}{\text{100 residues}}$)')
            ax.set_xticklabels(Figure_5c_df['xlabel'], rotation=45, ha='right')  # Rotate x-axis labels for better readability
                
        # Customize the plot
        print(f'Max(y) ceil: {maxy} {np.ceil(maxy)}')
        #ax.set_ylim(0.5,np.ceil(maxy))
        ax.set_xlim(-0.5, 2.5)

        ax.tick_params(axis='y', labelsize=6)
        ax.tick_params(axis='x', labelsize=6)
        # Remove the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.legend(fontsize='small', handlelength=1, handletextpad=0.5, borderpad=0.3, frameon=False)

        # add the curly brackets
        ax.annotate("(insignificant)", xy=(1.1, 0.2), xytext=(1.2, 0.2), xycoords='axes fraction', arrowprops=dict(arrowstyle='-[, widthB=4.0, lengthB=0.75', lw=0.5, color='#7E8D85'), va='center', ha='center', rotation=90, fontweight='bold', color='#7E8D85')
        return ax
    #################################################################################################################

    #################################################################################################################
    def plot_Figure_5d(self, ax):
        """
        Figure 5d (row 2 column 2):
            Loop closing contacts enrichement matrix
            C buffers at spa50 and LiPMScov 50 in EXP dataset
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Sequence_Complexity_and_Discrimination/Compare_OR_with_permutation/EXP/C_50_p100000/FrequencyGeneratorOutput/OR_FDR_pvalues.csv
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Sequence_Complexity_and_Discrimination/Compare_OR_with_permutation/EXP/C_50_p100000/FrequencyGeneratorOutput/OR_GT.csv
        """
        buff_tag = {'C': r'$\text{cyto-serum}_{\text{PSM}}$'}
        #######################################
        ## Load Figure 2c data
        inp = f'{self.slug_path}/Sequence_Complexity_and_Discrimination/Compare_OR_with_permutation/EXP/C_50_p100000/FrequencyGeneratorOutput/OR_FDR_pvalues.csv'
        print(f'inp: {inp}')
        Figure_5d_FDR_df = pd.read_csv(inp)
        Figure_5d_FDR_df.set_index('AA', inplace=True)
        print(f'Figure_5d_FDR_df:\n{Figure_5d_FDR_df}')
    
        #a = Figure_5d_FDR_df.values
        #symmetric_a = np.maximum(a, a.T)
        #Figure_5d_FDR_df = pd.DataFrame(symmetric_a, index=Figure_5d_FDR_df.index, columns=Figure_5d_FDR_df.columns)
        ## take care of pvalue elements that are not robust across all three buffers
        Figure_5d_FDR_df.at['C', 'C'] = 1 
        Figure_5d_FDR_df.at['F', 'K'] = 1
        Figure_5d_FDR_df.at['K', 'F'] = 1
        Figure_5d_FDR_df.at['M', 'R'] = 1
        Figure_5d_FDR_df.at['R', 'M'] = 1
        mask = np.triu(np.ones(Figure_5d_FDR_df.shape), k=1)  # k=1 excludes the diagonal
        Figure_5d_FDR_df = Figure_5d_FDR_df.mask(mask == 1)
        print(f'Figure_5d_FDR_df:\n{Figure_5d_FDR_df}')


        inp = f'{self.slug_path}/Sequence_Complexity_and_Discrimination/Compare_OR_with_permutation/EXP/C_50_p100000/FrequencyGeneratorOutput/OR_GT.csv'
        print(f'inp: {inp}')
        Figure_5d_OR_df = pd.read_csv(inp)
        #Figure_5d_OR_df = Figure_5d_OR_df.replace([np.inf, -np.inf], np.nan).fillna(0)
        Figure_5d_OR_df.set_index('AA', inplace=True)
        Figure_5d_OR_df.at['C', 'C'] = 0
        print(f'Figure_5d_OR_df:\n{Figure_5d_OR_df}')
        #######################################
      
        min_OR = min(Figure_5d_OR_df.values.flatten())
        max_OR = max(Figure_5d_OR_df.values.flatten())
        max_OR = 5
        
        # Create a TwoSlopeNorm for custom normalization
        norm = TwoSlopeNorm(vmin=min_OR, vcenter=1, vmax=max_OR)  # vcenter sets the whitest point

        # Replace `np.inf` with `max_OR` for visualization
        # Create a mask for the upper triangle
        visualization_df = Figure_5d_OR_df.replace(np.inf, max_OR).replace(0.0, np.nan)
        mask = np.triu(np.ones(visualization_df.shape), k=1)  # k=1 excludes the diagonal
        visualization_df = visualization_df.mask(mask == 1)
        print(visualization_df)

        # Create a DataFrame for annotations, replacing `np.inf` with 'inf'
        annotations = Figure_5d_OR_df.round(decimals=1).replace(np.inf, "inf").astype(str)
        annotations.at['C', 'C'] = 'NA'
        print(annotations)

        #sns.heatmap(Figure_5d_OR_df, annot=Figure_5d_OR_df.round(decimals=1).astype(str), fmt="", cmap="coolwarm", norm=norm, vmin=min_OR, vmax=max_OR, linewidths=.5, cbar_kws={"shrink": .8, "aspect": 30, "label": "Odds Ratio"}, annot_kws={"fontsize": 5}, ax=ax)
        sns.heatmap(visualization_df, annot=annotations, fmt="", cmap="coolwarm", norm=norm, vmin=min_OR, vmax=max_OR, linewidths=.5, cbar_kws={"shrink": .8, "aspect": 30, "label": "Odds Ratio"}, annot_kws={"fontsize": 5}, ax=ax)

        ## annotate the CC cell with N.A. in black
        # Get the location of the cell in the heatmap
        row_index = visualization_df.index.get_loc('C')
        col_index = visualization_df.columns.get_loc('C')
        cell_value = annotations.loc['C', 'C']
        print(row_index, col_index, cell_value)
        ax.text(col_index + 0.5, row_index + 0.5, f"{cell_value}", color="black", ha="center", va="center", fontsize=5)
    

        # if there is a valid pvalue dataframe color those with values below 0.05 in a yellow highlight
        if isinstance(Figure_5d_FDR_df, pd.DataFrame):
            print('Pvalues found')
            print(Figure_5d_FDR_df)
            # Highlight cells based on df2 values
            for i in range(Figure_5d_FDR_df.shape[0]):
                for j in range(Figure_5d_FDR_df.shape[1]):
                    if i == 9 and j == 2:
                        ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=False, edgecolor='black', lw=0.75, ls='--')) # make T-Y dashed
                    if Figure_5d_FDR_df.iloc[i, j] < 0.05:
                        ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=False, edgecolor='black', lw=0.75))

        # Set labels
        print(Figure_5d_OR_df.index)
        ax.set_yticklabels(Figure_5d_OR_df.index, rotation=0)
        print(Figure_5d_OR_df.columns)
        ax.set_xticks(np.arange(len(Figure_5d_OR_df.columns)) + 0.5)
        ax.set_xticklabels(Figure_5d_OR_df.columns, rotation=0)

        ax.set_xlabel('Amino Acid')
        ax.set_ylabel('Amino Acid')

        ax.tick_params(axis='y', labelsize=6, width=0.5)
        ax.tick_params(axis='x', labelsize=6, width=0.5)
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
    Creates Figure 5
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
    fig_height_mm = 70  # Height in mm
    #custom_font_path = "/storage/group/epo2/default/ims86/miniconda3/envs/FtoF/fonts/Arial.ttf" # Path to your custom font
    #arial_font = fm.FontProperties(fname=custom_font_path) # Create a FontProperties object
    plt.rcParams['font.family'] = 'Arial'  # Change to your desired font, e.g., 'Times New Roman', 'DejaVu Sans', etc.
    plt.rcParams['font.size'] = 6  # Default font size
    plt.rcParams['pdf.fonttype'] = 42
    fig, axs = plt.subplots(1, 4, 
                            figsize=(mm_to_inches(fig_width_mm), mm_to_inches(fig_height_mm)), 
                            dpi=600, constrained_layout=False, gridspec_kw={'width_ratios': [1, 1, 1, 3]})  

    plotter.plot_Figure_5a(axs[0])
    plotter.plot_Figure_5b(axs[1])
    plotter.plot_Figure_5c(axs[2])
    plotter.plot_Figure_5d(axs[3])

    #########################################################
    # Get the current position of the subplots and adjust their positions manually
    axs0_position = axs[0].get_position()
    width0, height0 = axs0_position.extents[2] - axs0_position.extents[0], axs0_position.extents[3] - axs0_position.extents[1]
    print(axs0_position, width0, height0)
    axs[0].set_position([0.025, 0.1575, width0, height0])  # [left, bottom, width, height]

    bbox_in_fig_coords = axs[0].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bbox_in_fig_coords.x0, 0.995, 'a', fontsize=8, fontweight='bold', va='top', ha='left')


    axs1_position = axs[1].get_position()
    width1, height1 = axs1_position.extents[2] - axs1_position.extents[0], axs1_position.extents[3] - axs1_position.extents[1]
    print(axs1_position, width1, height1)
    axs[1].set_position([0.235, 0.1575, 0.075, height0])

    bbox_in_fig_coords = axs[1].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bbox_in_fig_coords.x0, 0.995, 'b', fontsize=8, fontweight='bold', va='top', ha='left')


    axs2_position = axs[2].get_position()
    width2, height2 = axs2_position.extents[2] - axs2_position.extents[0], axs2_position.extents[3] - axs2_position.extents[1]
    print(axs2_position, width2, height2)
    axs[2].set_position([0.425, 0.1575, width0, height0])

    bbox_in_fig_coords = axs[2].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bbox_in_fig_coords.x0, 0.995, 'c', fontsize=8, fontweight='bold', va='top', ha='left')


    axs3_position = axs[3].get_position()
    width3, height3 = axs3_position.extents[2] - axs3_position.extents[0], axs3_position.extents[3] - axs3_position.extents[1]
    print(axs3_position, width3, height3)
    axs[3].set_position([0.625, 0.1575, 0.365, height0])

    bbox_in_fig_coords = axs[3].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bbox_in_fig_coords.x0, 0.995, 'd', fontsize=8, fontweight='bold', va='top', ha='left')

    # Adjust the colorbar position
    colorbar = axs[3].collections[0].colorbar
    colorbar.ax.set_position([0.935, 0.5, 0.02, 0.4])  # [left, bottom, width, height]
    colorbar.ax.tick_params(width=0.5)
    # Get the existing tick positions
    tick_positions = colorbar.get_ticks()
    #tick_positions = np.linspace(0, 5, 6)  # Equal spacing for ticks
    custom_labels = [' 0', ' 1', ' 2', ' 3', ' 4', '≥5']
    colorbar.set_ticks(tick_positions)
    print(tick_positions, custom_labels)
    colorbar.set_ticklabels(custom_labels)

    #########################################################

    # final formating and output
    # Automatically adjust layout
    #fig.tight_layout()  # 'pad' is the overall padding
    figure_outpath = os.path.join(args.out_path, 'Figure5.pdf')
    plt.savefig(figure_outpath)
    print(f'SAVED: {figure_outpath}')

    figure_outpath = os.path.join(args.out_path, 'Figure5.png')
    plt.savefig(figure_outpath)
    print(f'SAVED: {figure_outpath}')

    figure_outpath = os.path.join(args.out_path, 'Figure5.svg')
    plt.savefig(figure_outpath)
    print(f'SAVED: {figure_outpath}')
    print('NORMAL TERMINATION')

if __name__ == "__main__":
    main()

