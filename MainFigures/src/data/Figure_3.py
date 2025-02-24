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
    Figure 3a (row 1 column 1): 
        Association between non-refoldability and the presence of a native entanglement
        C buffers at spa50 and LiPMScov 50 in EXP dataset
        make a contingency table with the values and put the OR and pvalue below
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Association_Native_Entanglements_and_Misfolding/EntanglementsAndNonrefoldability/Plot_EntanglementsAndNonrefoldability/entanglement_and_nonrefoldability_plot_data_all_genes.csv

    Figure 3b (row 1 column 2):
        Association between change in proteolysis suseptibility and protein region (entangled versus non-entangled)
        C, CD, CG buffers at spa50 and LiPMScov 50 and has to have an entanglement
        C buffers at spa50 and LiPMScov 50 in EXP dataset
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
        
    Figure 3c: 
        PSM matched distributions of SASA
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Propensity_Score_Matching/EXP/matched_residue_features_res_sasa_n1.csv (Post matched data)
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gen_proteome_features_EXP/res_features_lib/ (Pre matched feature files)

    Figure 3d:
        Association between change in proteolysis suseptibility and protein region (entangled versus non-entangled)
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/PSM/ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv


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
    def plot_Figure_3a(self, ax):
        """
        Figure 3a (row 1 column 1): 
            Association between non-refoldability and the presence of a native entanglement
            C buffers at spa50 and LiPMScov 50 in EXP dataset
            make a contingency table with the values and put the OR and pvalue below
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Association_Native_Entanglements_and_Misfolding/EntanglementsAndNonrefoldability/Plot_EntanglementsAndNonrefoldability/entanglement_and_nonrefoldability_plot_data_all_genes.csv

        """
        buff_tag = {'C': 'cyto-serum', 'CD': '+DnaK', 'CG': '+GroEL'}
        #######################################
        ## Load Figure 1a data
        #../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Association_Native_Entanglements_and_Misfolding/EntanglementsAndNonrefoldability/Plot_EntanglementsAndNonrefoldability/entanglement_and_nonrefoldability_plot_data_all_genes.cs
        inp = f'{self.slug_path}/Association_Native_Entanglements_and_Misfolding/EntanglementsAndNonrefoldability/Plot_EntanglementsAndNonrefoldability/entanglement_and_nonrefoldability_plot_data_all_genes.csv'
        print(f'inp: {inp}')
        self.Figure_3a_df = pd.read_csv(inp)
        #print(f'Figure_3a_df:\n{self.Figure_3a_df}')
        #######################################

        Figure_3a_df = self.Figure_3a_df.copy()
        Figure_3a_df = Figure_3a_df[Figure_3a_df['spa'] == 50]
        Figure_3a_df = Figure_3a_df[Figure_3a_df['label'] == 'EXP']
        Figure_3a_df = Figure_3a_df[Figure_3a_df['buff'] == 'C']
        Figure_3a_df['xlabel'] = [buff_tag[b] for b in Figure_3a_df['buff'].values]
        print(Figure_3a_df)

        # Save the figure 1a raw plot df
        Figure_3a_outfile_csv = os.path.join(self.out_path, f'Figure_3a.csv')
        Figure_3a_df.to_csv(Figure_3a_outfile_csv)
        print(f'SAVED: {Figure_3a_outfile_csv}')
        print(Figure_3a_df.keys())

        # Extract the required values from the DataFrame
        values = [
            [Figure_3a_df['entY_nonrefoldedY'].values[0], Figure_3a_df['entY_nonrefoldedN'].values[0]],
            [Figure_3a_df['entN_nonrefoldedY'].values[0], Figure_3a_df['entN_nonrefoldedN'].values[0]],
        ]
        print(values)
        # make odds ratio and pvalue strings
        OR = Figure_3a_df['OR'].values[0]
        OR_lb, OR_ub = Figure_3a_df['OR_lb'].values[0], Figure_3a_df['OR_ub'].values[0]
        OR_str = f'Odds Ratio = {OR:.2f} ({OR_lb:.2f}, {OR_ub:.2f})'

        pvalue = Figure_3a_df['pvalue'].values[0]
        pvalue_str = f'p-value = {format_scientific(pvalue)}'

        # Define row and column labels
        row_labels = ['Yes ', 'No ']
        col_labels = ['Yes', 'No']

        # Hide the axes
        ax.axis('tight')
        ax.axis('off')

        # Create the table
        bbox = [0.3, 0.4, 0.6, 0.3]
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
        ax.text(-0.19, bbox[1] + 0.1, "Entangled\n  proteins", transform=ax.transAxes, va='center', ha='left', rotation=90, fontsize=6, fontweight='bold')  

        # Add column header text to the top of the table
        ax.text(0.6,  bbox[1] + 0.375, "Non-Refoldable", transform=ax.transAxes, va='center', ha='center', fontsize=6, fontweight='bold') 
    
        # Add OR and pvalue strings below table
        ax.text(0.5, bbox[1] - 0.15, OR_str, transform=ax.transAxes, va='center', ha='center', fontsize=6)  
        ax.text(0.5, bbox[1] - 0.25, pvalue_str, transform=ax.transAxes, va='center', ha='center', fontsize=6) 

        ax.tick_params(axis='y', labelsize=6)
        ax.tick_params(axis='x', labelsize=6)
        return ax
    #################################################################################################################

    #################################################################################################################
    def plot_Figure_3b(self, ax):
        """
        Figure 3b (row 1 column 2):
            Association between change in proteolysis suseptibility and protein region (entangled versus non-entangled)
            C, CD, CG buffers at spa50 and LiPMScov 50 and has to have an entanglement
            C buffers at spa50 and LiPMScov 50 in EXP dataset
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
        """
        buff_tag = {'C': 'cyto-serum', 'CD': '+DnaK', 'CG': '+GroEL'}
        #######################################
        ## Load Figure 1a data
        inp = f'{self.slug_path}/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/whole_proteome/ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv'
        print(f'inp: {inp}')
        self.Figure_3b_df = pd.read_csv(inp)
        print(f'Figure_3b_df:\n{self.Figure_3b_df}')
        #######################################

        Figure_3b_df = self.Figure_3b_df.copy()
        Figure_3b_df = Figure_3b_df[Figure_3b_df['spa'] == 50]
        Figure_3b_df['xlabel'] = [buff_tag[b] for b in Figure_3b_df['buff'].values]
        #print(Figure_3b_df)

        # Define the custom order for the 'buff' column
        custom_order = ['C', 'CD', 'CG']
        Figure_3b_df['buff'] = pd.Categorical(Figure_3b_df['buff'], categories=custom_order, ordered=True)
        Figure_3b_df = Figure_3b_df.sort_values('buff')
        print(Figure_3b_df)
      
        # Save the figure 1a raw plot df
        Figure_3b_outfile_csv = os.path.join(self.out_path, f'Figure_3b.csv')
        Figure_3b_df.to_csv(Figure_3b_outfile_csv)
        print(f'SAVED: {Figure_3b_outfile_csv}')
                
        # Adjust the linewidth of the axis spines
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)  # Set the linewidth for all spines

        # Adjust the linewidth of the ticks
        ax.tick_params(width=0.5)  # Both major and minor ticks

        maxy = 0
        # Calculate error bars
        yerr_lower = Figure_3b_df['OR'] - Figure_3b_df['OR_lb']
        yerr_upper = Figure_3b_df['OR_ub'] - Figure_3b_df['OR']
        yerr = [yerr_lower, yerr_upper]
        
        # Plot the trace with error bars
        ax.errorbar(
            Figure_3b_df['xlabel'], 
            Figure_3b_df['OR'], 
            yerr=yerr, 
            fmt='o',  # Line and point markers
            capsize=3,  # Add caps to the error bars
            color='black',
            markersize=3,
            markerfacecolor='black',  
            elinewidth=0.5,  
            markeredgewidth=0.5
        )

        ## mark pvalue annotations
        # Annotate significant p-values with '*'
        for x, y, err, p in zip(Figure_3b_df['xlabel'], Figure_3b_df['OR'], yerr_upper, Figure_3b_df['pvalues']):
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
        ax.axhline(y=1, color='#7E8D85', linestyle='--', linewidth=0.5)
        ax.set_ylabel("Odds ratio between entangled\nprotein regions and misfolding")
        ax.set_xticklabels(Figure_3b_df['xlabel'], rotation=45, ha='right')  # Rotate x-axis labels for better readability
            
        # Customize the plot
        print(f'Max(y) ceil: {maxy} {np.ceil(maxy)}')
        ax.set_ylim(0.5,np.ceil(maxy))
        ax.set_xlim(-0.5, 2.5)

        ax.tick_params(axis='y', labelsize=6)
        ax.tick_params(axis='x', labelsize=6)
        # Remove the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        return ax   
    #################################################################################################################

    #################################################################################################################
    def plot_Figure_3c(self, pre_ax, post_ax):
        """
    Figure 3c: 
        PSM matched distributions of SASA
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Propensity_Score_Matching/EXP/matched_residue_features_res_sasa_n1.csv (Post matched data)
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gen_proteome_features_EXP/res_features_lib/ (Pre matched feature files)
        """
        buff_tag = {'C': 'cyto-serum', 'CD': '+DnaK', 'CG': '+GroEL'}
        #######################################
        ## Load Figure 3c data
        # Load post matched res_sasa
        inp = f'{self.slug_path}/Modeling_Odds_of_Misfolding/Propensity_Score_Matching/EXP/matched_residue_features_res_sasa_n1.csv'
        print(f'inp: {inp}')
        Figure_3c_post_df = pd.read_csv(inp, sep='|', low_memory=False)
        Figure_3c_post_df = Figure_3c_post_df[['region', 'res_sasa']]
        Figure_3c_post_df['PSM'] = 'Post'
        print(Figure_3c_post_df)

        post_nonentR = Figure_3c_post_df[Figure_3c_post_df['region'] == 0]['res_sasa'].values
        print(f'post_nonentR: {post_nonentR} {post_nonentR.shape}')

        # Load pre matched res_sasa
        feature_files = glob.glob(f'{self.slug_path}/Make_Protein_Feature_Files/Gen_proteome_features_EXP/res_features_lib/*_resfeatures.csv')
        Figure_3c_pre_df = [pd.read_csv(f, sep='|', low_memory=False) for f in feature_files]
        Figure_3c_pre_df = pd.concat(Figure_3c_pre_df)
        Figure_3c_pre_df = Figure_3c_pre_df[['region', 'res_sasa']]
        Figure_3c_pre_df['PSM'] = 'Pre'
        print(f'Figure_3c_pre_df:\n{Figure_3c_pre_df}')

        pre_nonentR = Figure_3c_pre_df[Figure_3c_pre_df['region'] == 0]['res_sasa'].values
        print(f'pre_nonentR: {pre_nonentR} {pre_nonentR.shape}')

        pre_entR = Figure_3c_pre_df[Figure_3c_pre_df['region'] == 1]['res_sasa'].values
        print(f'pre_entR: {pre_entR} {pre_entR.shape}')
        #######################################

        Figure_3c_df = pd.concat([Figure_3c_post_df, Figure_3c_pre_df])
        print(f'Figure_3c_df:\n{Figure_3c_df}')

        # Save the figure 1a raw plot df
        Figure_3c_outfile_csv = os.path.join(self.out_path, f'Figure_3c.csv')
        Figure_3c_df.to_csv(Figure_3c_outfile_csv)
        print(f'SAVED: {Figure_3c_outfile_csv}')

        # Adjust the linewidth of the axis spines
        for spine in pre_ax.spines.values():
            spine.set_linewidth(0.5)  # Set the linewidth for all spines

        # Adjust the linewidth of the ticks
        pre_ax.tick_params(width=0.5)  # Both major and minor ticks

        # Compute CDFs for each dataset
        x_pre_entR, cdf_pre_entR = compute_cdf(pre_entR)
        pre_entR_df = pd.DataFrame({'SASA':x_pre_entR, 'CDF':cdf_pre_entR})
        Figure_3c_pre_entR_csv = os.path.join(self.out_path, f'Figure_3c_pre_entR.csv')
        pre_entR_df.to_csv(Figure_3c_pre_entR_csv)
        print(f'SAVED: {Figure_3c_pre_entR_csv}')

        x_pre_nonentR, cdf_pre_nonentR = compute_cdf(pre_nonentR)
        pre_nonentR_df = pd.DataFrame({'SASA':x_pre_nonentR, 'CDF':cdf_pre_nonentR})
        Figure_3c_pre_nonentR_csv = os.path.join(self.out_path, f'Figure_3c_pre_nonentR.csv')
        pre_nonentR_df.to_csv(Figure_3c_pre_nonentR_csv)
        print(f'SAVED: {Figure_3c_pre_nonentR_csv}')

        x_post_nonentR, cdf_post_nonentR = compute_cdf(post_nonentR)
        post_nonentR_df = pd.DataFrame({'SASA':x_post_nonentR, 'CDF':cdf_post_nonentR})
        Figure_3c_post_nonentR_csv = os.path.join(self.out_path, f'Figure_3c_post_nonentR.csv')
        post_nonentR_df.to_csv(Figure_3c_post_nonentR_csv)
        print(f'SAVED: {Figure_3c_post_nonentR_csv}')
       
        # Plot the CDFs
        #ax.plot(x_pre_entR, cdf_pre_entR, color='black', label='Pre entR', linewidth=0.75)
        #ax.plot(x_pre_nonentR, cdf_pre_nonentR, color='blue', label='Pre nonentR', linewidth=0.75)
        #ax.plot(x_post_nonentR, cdf_post_nonentR, color='red', linestyle='--', label='Post nonentR', linewidth=0.75)

        # Plot the PDFs
        # Compute the histogram for all three datasets using the same bins
        all_data = np.concatenate([pre_entR, pre_nonentR, post_nonentR])
        bin_edges = np.histogram_bin_edges(all_data, bins=30)

        # Plot the first PDF (pre_entR)
        pre_ax.hist(pre_entR, bins=bin_edges, density=True, alpha=0.5, label='Entangled', edgecolor='none')

        # Plot the second PDF (pre_nonentR)
        pre_ax.hist(pre_nonentR, bins=bin_edges, density=True, alpha=0.5, label='Non-Entangled', edgecolor='black', histtype='step')

        # Add labels, legend, and title
        pre_ax.set_xlabel(r'SASA ($\text{nm}^2$)')
        pre_ax.set_ylabel('Probability Density Function')
        pre_ax.legend(fontsize='small', handlelength=1, handletextpad=0.5, borderpad=0.3, frameon=False, bbox_to_anchor=[0.05, 1])

        pre_ax.tick_params(axis='y', labelsize=6)
        pre_ax.tick_params(axis='x', labelsize=6)
        # Remove the right and top spines
        pre_ax.spines['right'].set_visible(False)
        pre_ax.spines['top'].set_visible(False)
        pre_ax.set_xlim(0,2)


        # Compute the histogram for post_nonentR and plot as a cityscape trace
        hist, bin_edges = np.histogram(post_nonentR, bins=bin_edges, density=True)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        post_ax.hist(pre_entR, bins=bin_edges, density=True, alpha=0.5, label='Entangled', edgecolor='none')
        post_ax.hist(post_nonentR, bins=bin_edges, density=True, alpha=0.5, label='Non-Entangled', edgecolor='black', histtype='step')

        # Add labels, legend, and title
        post_ax.set_xlabel(r'SASA ($\text{nm}^2$)')
        post_ax.set_ylabel('Probability Density Function')
        post_ax.legend(fontsize='small', handlelength=1, handletextpad=0.5, borderpad=0.3, frameon=False, bbox_to_anchor=[0.05, 0.825])
        post_ax.text(0.1, 1, "propensity\nscore\nmatched", transform=post_ax.transAxes, va='top', ha='left', fontsize=6) 

        post_ax.tick_params(axis='y', labelsize=6)
        post_ax.tick_params(axis='x', labelsize=6)
        # Remove the right and top spines
        post_ax.spines['right'].set_visible(False)
        post_ax.spines['top'].set_visible(False)
        post_ax.set_xlim(0,2)
        #return ax
    #################################################################################################################

    #################################################################################################################
    def plot_Figure_3d(self, ax):
        """
        Figure 3d:
            Association between change in proteolysis suseptibility and protein region (entangled versus non-entangled)
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/PSM/ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv
        """
        buff_tag = {'C': 'cyto-serum'}
        #######################################
        ## Load Figure 3c data
        inp = f'{self.slug_path}/Modeling_Odds_of_Misfolding/Regressions/Plots/EXP/PSM/ent_genes_Rall_binomial_regression_results_var-region_LiPMScov50.csv'
        print(f'inp: {inp}')
        self.Figure_3d_df = pd.read_csv(inp)
        print(f'Figure_3d_df:\n{self.Figure_3d_df}')
        #######################################
    
        Figure_3d_df = self.Figure_3d_df.copy()
        Figure_3d_df = Figure_3d_df[Figure_3d_df['spa'] == 50]
        Figure_3d_df = Figure_3d_df[Figure_3d_df['buff'] == 'C']
        Figure_3d_df['xlabel'] = [buff_tag[b] for b in Figure_3d_df['buff'].values]
        print(Figure_3d_df)
      
        # Save the figure 1a raw plot df
        Figure_3d_outfile_csv = os.path.join(self.out_path, f'Figure_3d.csv')
        Figure_3d_df.to_csv(Figure_3d_outfile_csv)
        print(f'SAVED: {Figure_3d_outfile_csv}')
        
        # Adjust the linewidth of the axis spines
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)  # Set the linewidth for all spines

        # Adjust the linewidth of the ticks
        ax.tick_params(width=0.5)  # Both major and minor ticks

        maxy = 0
        # Calculate error bars
        yerr_lower = Figure_3d_df['OR'] - Figure_3d_df['OR_lb']
        yerr_upper = Figure_3d_df['OR_ub'] - Figure_3d_df['OR']
        yerr = [yerr_lower, yerr_upper]
        
        # Plot the trace with error bars
        ax.errorbar(
            Figure_3d_df['xlabel'], 
            Figure_3d_df['OR'], 
            yerr=yerr, 
            fmt='o',  # Line and point markers
            capsize=3,  # Add caps to the error bars
            color='black',
            markersize=3,
            markerfacecolor='black',  
            elinewidth=0.5,  
            markeredgewidth=0.5
        )

        ## mark pvalue annotations
        # Annotate significant p-values with '*'
        for x, y, err, p in zip(Figure_3d_df['xlabel'], Figure_3d_df['OR'], yerr_upper, Figure_3d_df['pvalues']):
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
                y + err + 0.05,  # Place the annotation slightly above the upper error bar
                p_annot, 
                fontsize=6, 
                ha='center', 
                color='black')                  
            #print(p_annot, y + err + 0.1)

            if y + err + 0.1 > maxy:
                maxy = y + err + 0.1
            
        # Add a dashed black line at y=1    
        ax.axhline(y=1, color='#7E8D85', linestyle='--', linewidth=0.5)
        ax.set_ylabel("Odds ratio between entangled\nprotein regions and misfolding")
        #ax.set_xticklabels(Figure_3d_df['xlabel'], rotation=0, ha='right')  # Rotate x-axis labels for better readability
            
        # add PSM label
        ax.text(0.05, 1, "propensity\nscore\nmatched", transform=ax.transAxes, va='top', ha='left', fontsize=6) 

        # Customize the plot
        print(f'Max(y) ceil: {maxy} {np.ceil(maxy)}')
        ax.set_ylim(0.5,np.ceil(maxy))
        ax.set_xlim(-0.5, 0.5)
        ax.tick_params(axis='y', labelsize=6)
        ax.tick_params(axis='x', labelsize=6)
        # Remove the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
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
        ax.text(0.81, 0.65, 'Positive\nAssociation', rotation=90, va='center', ha='center', color='#7E8D85')

        # Arrow pointing down (Negative Association)
        ax.annotate(
            '',  # No text for the arrow itself
            xy=(0.5, 0),  # End point of the arrow (bottom center)
            xytext=(0.5, 0.32),  # Start point of the arrow (middle)
            color='#7E8D85',
            arrowprops=dict(facecolor='#7E8D85', edgecolor='#7E8D85', arrowstyle='->', lw=0.5)
        )
        # Label for the downward arrow
        ax.text(0.81, 0.175, 'Negative\nAssociation', rotation=90, va='center', ha='center', color='#7E8D85')

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
    fig_width_mm = 89  # Width in mm
    fig_height_mm = 100  # Height in mm
    #custom_font_path = "/storage/group/epo2/default/ims86/miniconda3/envs/FtoF/fonts/Arial.ttf" # Path to your custom font
    #arial_font = fm.FontProperties(fname=custom_font_path) # Create a FontProperties object
    plt.rcParams['font.family'] = 'Arial'  # Change to your desired font, e.g., 'Times New Roman', 'DejaVu Sans', etc.
    plt.rcParams['font.size'] = 6  # Default font size
    plt.rcParams['pdf.fonttype'] = 42
    #fig, axs = plt.subplots(1, 5, figsize=(mm_to_inches(fig_width_mm), mm_to_inches(fig_height_mm)), dpi=600)  

    fig = plt.figure(figsize=(mm_to_inches(fig_width_mm), mm_to_inches(fig_height_mm)))
    gs = gridspec.GridSpec(3, 3, figure=fig)
    # Add subplots
    ax1 = fig.add_subplot(gs[0, 0])  # First row, first column
    ax2 = fig.add_subplot(gs[0, 1])  # First row, second column
    ax3 = fig.add_subplot(gs[0, 2])  # First row, second column
    ax4 = fig.add_subplot(gs[1, 0])  # Second row, first column
    ax5 = fig.add_subplot(gs[1, 1])  # Second row, second column
    ax6 = fig.add_subplot(gs[1, 2])  # Second row, third column
    axs = [ax1, ax2, ax3, ax4, ax5, ax6]

    ## make subplot figures 
    plotter.plot_Figure_3a(axs[0])
    plotter.plot_Figure_3b(axs[1])
    plotter.plot_arrows(axs[2])
    plotter.plot_Figure_3c(axs[3], axs[4])
    plotter.plot_Figure_3d(axs[5])

    #########################################################
    # Get the current position of the subplots and adjust their positions manually
    ## panel a
    axs0_position = axs[0].get_position()
    x00, y00 = axs0_position.extents[0], axs0_position.extents[1]
    width0, height0 = axs0_position.extents[2] - axs0_position.extents[0], axs0_position.extents[3] - axs0_position.extents[1]
    axs[0].set_position([0.05, 0.6, width0, 0.35])  # [left, bottom, width, height]
    print('axs[0]', axs0_position, width0, 0.35)

    bbox_in_fig_coords = axs[0].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bbox_in_fig_coords.x0, 1, 'a', fontsize=8, fontweight='bold', va='top', ha='left')

    ## panel b
    axs1_position = axs[1].get_position()
    x01, y01 = axs1_position.extents[0], axs1_position.extents[1]
    width1, height1 = axs1_position.extents[2] - axs1_position.extents[0], axs1_position.extents[3] - axs1_position.extents[1]
    axs[1].set_position([0.52, 0.6, 0.35, 0.35])
    print('axs[1]', axs1_position, width1, 0.35)

    bbox_in_fig_coords = axs[1].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bbox_in_fig_coords.x0, 1, 'b', fontsize=8, fontweight='bold', va='top', ha='left')

    ## arrow panel
    axs2_position = axs[2].get_position()
    x02, y02 = axs2_position.extents[0], axs2_position.extents[1]
    width2, height2 = axs2_position.extents[2] - axs2_position.extents[0], axs2_position.extents[3] - axs2_position.extents[1]
    axs[2].set_position([0.85, 0.6, 0.1, 0.35])
    print('axs[2]', axs2_position, 0.1, 0.35)

    ## panel c
    axs3_position = axs[3].get_position()
    x03, y03 = axs3_position.extents[0], axs3_position.extents[1]
    width3, height3 = axs3_position.extents[2] - axs3_position.extents[0], axs3_position.extents[3] - axs3_position.extents[1]
    axs[3].set_position([0.085, 0.09, width3, 0.35])
    print('axs[3]', axs3_position, width3, 0.35)

    bbox_in_fig_coords = axs[3].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bbox_in_fig_coords.x0, 0.48, 'c', fontsize=8, fontweight='bold', va='top', ha='left')


    axs4_position = axs[4].get_position()
    x04, y04 = axs4_position.extents[0], axs4_position.extents[1]
    width4, height4 = axs4_position.extents[2] - axs4_position.extents[0], axs4_position.extents[3] - axs4_position.extents[1]
    axs[4].set_position([0.4225, 0.09, width4, 0.35])
    print('axs[4]', axs4_position, width4, 0.35)

    ## panel d
    axs5_position = axs[5].get_position()
    x05, y05 = axs5_position.extents[0], axs5_position.extents[1]
    width5, height5 = axs5_position.extents[2] - axs5_position.extents[0], axs5_position.extents[3] - axs5_position.extents[1]
    axs[5].set_position([0.8175, 0.09, 0.15, 0.35])
    print('axs[5]', axs5_position, width5, 0.35)

    bbox_in_fig_coords = axs[5].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bbox_in_fig_coords.x0, 0.48, 'd', fontsize=8, fontweight='bold', va='top', ha='left')
    # Adjust the colorbar position

    #########################################################

    # final formating and output
    # Automatically adjust layout
    #fig.tight_layout(pad=2.0)  # 'pad' is the overall padding
    figure_outpath = os.path.join(args.out_path, 'Figure3.pdf')
    plt.savefig(figure_outpath)
    print(f'SAVED: {figure_outpath}')

    figure_outpath = os.path.join(args.out_path, 'Figure3.svg')
    plt.savefig(figure_outpath)
    print(f'SAVED: {figure_outpath}')

    figure_outpath = os.path.join(args.out_path, 'Figure3.png')
    plt.savefig(figure_outpath)
    print(f'SAVED: {figure_outpath}')

    print('NORMAL TERMINATION')

if __name__ == "__main__":
    main()

