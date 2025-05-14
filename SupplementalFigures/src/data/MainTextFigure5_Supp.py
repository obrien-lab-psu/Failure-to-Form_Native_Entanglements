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
    Figure 4a (row 1 column 1): 
        Association between essentiality and being a client of DnaK or GroEL
        C buffers at spa50 and LiPMScov 50 in EXP dataset
        make a contingency table with the values and put the OR and pvalue below
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Chaperone_Client_Associations/Assoc_Client_n_Essential/Plots/Assoc_Client_n_Essential_plot_data_Knockout_clientType-DnaKonly.csv
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Chaperone_Client_Associations/Assoc_Client_n_Essential/Plots/Assoc_Client_n_Essential_plot_data_Knockout_clientType-GroELonly.csv

    Figure 4b (row 1 column 2):
        Entanglement complexity as a discrimination tool
        C, CD, CG buffers at spa50 and LiPMScov 50 and has to have an entanglement
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Entanglement_Topological_Complexity_and_Discrimination/EXP/Lasso_plot_results_C_50.csv
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Entanglement_Topological_Complexity_and_Discrimination/EXP/Lasso_plot_results_CD_50.csv
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Entanglement_Topological_Complexity_and_Discrimination/EXP/Lasso_plot_results_CG_50.csv
        
    Figure 4c: 
        DnaK binding motif analysis
        ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Chaperone_Client_Associations/Dnak_simple_motif_scan/Plots/consolidated_Dnak_binding_data_EXP.csv       

    Figure 4d (row 2 column 2):
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
    def plot_MainTextFigure4a_Supp(self, chap='DnaK'):
        """
        Figure 4a (row 1 column 1): 
            Association between essentiality and being a client of DnaK or GroEL
            C buffers at spa50 and LiPMScov 50 in EXP dataset
            make a contingency table with the values and put the OR and pvalue below
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Chaperone_Client_Associations/Assoc_Client_n_Essential/Plots/Assoc_Client_n_Essential_plot_data_Knockout_clientType-DnaKonly.csv
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Chaperone_Client_Associations/Assoc_Client_n_Essential/Plots/Assoc_Client_n_Essential_plot_data_Knockout_clientType-GroELonly.csv
        """
        buff_tag = {'C': 'cyto-serum', 'CD': '+DnaK', 'CG': '+GroEL'}
        #######################################
        ## Load Figure 1a data
        inp = f'{self.slug_path}/Chaperone_Client_Associations/Assoc_Client_n_Essential/Plots/Assoc_Client_n_Essential_plot_data_Knockout_clientType-{chap}only.csv'
        print(f'inp: {inp}')
        Figure_4a_df = pd.read_csv(inp)
        print(Figure_4a_df)
        Figure_4a_EXP_df = Figure_4a_df[(Figure_4a_df['LiPMScov'] == 50) & (Figure_4a_df['label'] == 'EXP')]
        Figure_4a_EXP_df['Chap'] = f'+{chap}'
        #print(f'Figure_4a_df:\n{Figure_4a_df}')

        Figure_4a_AF_df = Figure_4a_df[(Figure_4a_df['LiPMScov'] == 50) & (Figure_4a_df['label'] == 'AF')]
        Figure_4a_AF_df['Chap'] = f'+{chap}'
        #print(f'Figure_4a_GroEL_df:\n{Figure_4a_GroEL_df}')

        # Save the figure 1a raw plot df
        Figure_4a_df = pd.concat([Figure_4a_EXP_df, Figure_4a_AF_df])
        print(Figure_4a_df)
        Figure_4a_outfile_csv = os.path.join(self.out_path, f'MainTextFigure4a_chap-{chap}_Supp.csv')
        Figure_4a_df.to_csv(Figure_4a_outfile_csv)
        print(f'SAVED: {Figure_4a_outfile_csv}')
        #######################################

        #######################################
        ## set up figure(s)
        fig_width_mm = 80  # Width in mm
        fig_height_mm = 150  # Height in mm
        fig, axes = plt.subplots(3, 2, figsize=(mm_to_inches(fig_width_mm), mm_to_inches(fig_height_mm)), gridspec_kw={'width_ratios': [1, 0.175]})
        #######################################

        #######################################
        ## plot data
        lines = []
        for buffi, buff in enumerate(['C']):
            buff_df = Figure_4a_df[Figure_4a_df['buff'] == buff]
            #print(buff_df)

            for label in ['EXP', 'AF']:
                label_df = buff_df[buff_df['label'] == label]
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
                    #axes[0, 0].text(0.05, 0.98, "propensity score\nmatched", transform=axes[0, 0].transAxes, va='top', ha='left', fontsize=6)

                ## Plot the pvalues
                y = label_df['pvalue'].values
                axes[1, buffi].plot(x, y, marker='o', ls='none', markersize=3,  markeredgewidth=0.5, label=label)
                axes[1, buffi].set_yscale('log')
                axes[1, buffi].axhline(y=0.05, color='black', linestyle='--', linewidth=0.5)

                ## Plot the pvalues
                y = label_df['n'].values
                axes[2, buffi].plot(x, y, marker='o', ls='none', markersize=3,  markeredgewidth=0.5, label=label)
        #######################################

        ######################################
        ## adjust axes locations and labels
        panel_labels = np.asarray([['a', 'b', 'c'], ['d', 'e', 'f'], ['g', 'h', 'i']])
        panel_label_ys = {0:1, 1:0.66, 2:0.325}
        x0s = {0:0.175, 1:0.365, 2:2/3}
        y0s = {0:0.74, 1:0.4, 2:0.07}
        for coli, colname in enumerate(['cyto-serum']):
            for rowi, rowname in enumerate([f'Odds ratio for association of\nessentiality and being a known {chap} client', 'p-value', 'Number of proteins in dataset']):

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
                #if rowi == 0:
                #    axes[rowi, coli].set_ylim(0.5, 2)

                axs_position = axes[rowi, coli].get_position()
                x0 = x0s[coli]
                y0 = y0s[rowi]
                width0, height0 = axs_position.extents[2] - axs_position.extents[0], axs_position.extents[3] - axs_position.extents[1]
                #print(axs_position, width0, height0)
                width0 = 0.5
                axes[rowi, coli].set_position([x0, y0, width0, height0])  # [left, bottom, width, height]

                bbox_in_fig_coords = axes[rowi, coli].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
                fig.text(bbox_in_fig_coords.x0, panel_label_ys[rowi], panel_labels[rowi, coli], fontsize=8, fontweight='bold', va='top', ha='left')
        #######################################

        axes[0, 1].set_axis_off()
        axes[1, 1].set_axis_off()
        axes[2, 1].set_axis_off()

        # Adjust layout for clarity
        labels = ['Crystal Strutures', 'AlphaFold']
        fig.legend(lines, labels, loc='upper right', bbox_to_anchor=(1.009, 1), title=f"Dataset\n+{chap}", frameon=False)

        #######################################
        figure_outpath = os.path.join(self.out_path, f'MainTextFigure4a_chap-{chap}_Supp.pdf')
        plt.savefig(figure_outpath)
        print(f'SAVED: {figure_outpath}')

        figure_outpath = os.path.join(self.out_path, f'MainTextFigure4a_chap-{chap}_Supp.png')
        plt.savefig(figure_outpath)
        print(f'SAVED: {figure_outpath}')

        figure_outpath = os.path.join(self.out_path, f'MainTextFigure4a_chap-{chap}_Supp.svg')
        plt.savefig(figure_outpath)
        print(f'SAVED: {figure_outpath}')
        ####################################### 
    #################################################################################################################

    #################################################################################################################
    def plot_MainTextFigure4b_Supp(self,):
        """
        Figure 4b (row 1 column 2):
                Entanglement complexity as a discrimination tool
                C, CD, CG buffers at spa50 and LiPMScov 50 and has to have an entanglement
                ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Entanglement_Topological_Complexity_and_Discrimination/EXP/Lasso_plot_results_C_50.csv
                ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Entanglement_Topological_Complexity_and_Discrimination/EXP/Lasso_plot_results_CD_50.csv
                ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Entanglement_Topological_Complexity_and_Discrimination/EXP/Lasso_plot_results_CG_50.csv
        """
        buff_tag = {'C': 'cyto-serum', 'CD': '+DnaK', 'CG': '+GroEL'}
        #######################################
        ## Load Figure 1a data
        Figure_4b_df = []
        for dataset in ['EXP', 'AF']:
            for spa in [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]:
                inp = f'{self.slug_path}/Entanglement_Topological_Complexity_and_Discrimination/{dataset}/Lasso_plot_results_C_{spa}.csv'
                print(f'inp: {inp}')
                Figure_4b_spa_C_df = pd.read_csv(inp, sep='|')
                Figure_4b_spa_C_df = Figure_4b_spa_C_df[np.isclose(Figure_4b_spa_C_df['C'], 1.5)]
                Figure_4b_spa_C_df['buff'] = 'C'
                #print(f'Figure_4b_spa_C_df:\n{Figure_4b_spa_C_df}')


                inp = f'{self.slug_path}/Entanglement_Topological_Complexity_and_Discrimination/{dataset}/Lasso_plot_results_CD_{spa}.csv'
                print(f'inp: {inp}')
                Figure_4b_spa_CD_df = pd.read_csv(inp, sep='|')
                Figure_4b_spa_CD_df = Figure_4b_spa_CD_df[np.isclose(Figure_4b_spa_CD_df['C'], 1.5)]
                Figure_4b_spa_CD_df['buff'] = 'CD'
                #print(f'Figure_4b_spa_CD_df:\n{Figure_4b_spa_CD_df}')

                inp = f'{self.slug_path}/Entanglement_Topological_Complexity_and_Discrimination/{dataset}/Lasso_plot_results_CG_{spa}.csv'
                print(f'inp: {inp}')
                Figure_4b_spa_CG_df = pd.read_csv(inp, sep='|')
                Figure_4b_spa_CG_df = Figure_4b_spa_CG_df[np.isclose(Figure_4b_spa_CG_df['C'], 1.5)]
                Figure_4b_spa_CG_df['buff'] = 'CG'
                #print(f'Figure_4b_spa_CG_df:\n{Figure_4b_spa_CG_df}')
            
                Figure_4b_spa_df = pd.concat([Figure_4b_spa_C_df, Figure_4b_spa_CD_df, Figure_4b_spa_CG_df])
                Figure_4b_spa_df['xlabel'] = [buff_tag[b] for b in Figure_4b_spa_df['buff'].values]
                Figure_4b_spa_df['spa'] = spa
                Figure_4b_spa_df['label'] = dataset
                #print(Figure_4b_spa_df)
                Figure_4b_df += [Figure_4b_spa_df]
           
        # Save the figure 1a raw plot df
        Figure_4b_df = pd.concat(Figure_4b_df)
        Figure_4b_df['pvalue'] = np.where(Figure_4b_df['pvalue'].values == 0, 0.0001, Figure_4b_df['pvalue'].values)
        print(f'Figure_4b_df:\n{Figure_4b_df}')
        Figure_4b_outfile_csv = os.path.join(self.out_path, f'MainTextFigure4b_Supp.csv')
        Figure_4b_df.to_csv(Figure_4b_outfile_csv)
        print(f'SAVED: {Figure_4b_outfile_csv}')
        #######################################
 
        #######################################
        ## set up figure(s)
        fig_width_mm = 180  # Width in mm
        fig_height_mm = 120  # Height in mm
        fig, axes = plt.subplots(2, 4, figsize=(mm_to_inches(fig_width_mm), mm_to_inches(fig_height_mm)), gridspec_kw={'width_ratios': [1, 1, 1, 0.175]})
        #######################################

        #######################################
        ## plot data
        lines = []
        for buffi, buff in enumerate(['C', 'CD', 'CG']):
            buff_df = Figure_4b_df[Figure_4b_df['buff'] == buff]
            #print(buff_df)

            for label in ['EXP', 'AF']:
                label_df = buff_df[buff_df['label'] == label]
                #print(label)
                x = label_df['spa'].values

                ## Plot the OR
                y = label_df['<BA>'].values
                ylb = y - label_df['BA_lb'].values
                yub = label_df['BA_ub'].values - y
                yerr = [ylb, yub]
                line1 = axes[0, buffi].errorbar(x, y, yerr=yerr, fmt='o', capsize=3, markersize=3,  elinewidth=0.5,  markeredgewidth=0.5, label=label)
                #axes[0, buffi].set_yscale('log')
                axes[0, buffi].axhline(y=0.5, color='black', linestyle='--', linewidth=0.5)
                axes[0, buffi].set_ylim(0,1)
                if buffi == 0:
                    lines += [line1]

                ## Plot the pvalues
                y = label_df['pvalue'].values
                axes[1, buffi].plot(x, y, marker='o', ls='none', markersize=3,  markeredgewidth=0.5, label=label)
                axes[1, buffi].set_yscale('log')
                axes[1, buffi].axhline(y=0.05, color='black', linestyle='--', linewidth=0.5)
        #######################################

        ######################################
        ## adjust axes locations and labels
        panel_labels = np.asarray([['a', 'b', 'c'], ['d', 'e', 'f']])
        panel_label_ys = {0:1, 1:0.49}
        x0s = {0:0.075, 1:0.365, 2:2/3}
        y0s = {0:0.6, 1:0.07}
        for coli, colname in enumerate(['cyto-serum', '+DnaK', '+GroEL']):
            for rowi, rowname in enumerate(['Average balanaced accuracy', 'p-value']):

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
                #if rowi == 0:
                #    axes[rowi, coli].set_ylim(0.5, 2)

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

        # Adjust layout for clarity
        labels = ['Crystal Strutures', 'AlphaFold']
        fig.legend(lines, labels, loc='upper right', bbox_to_anchor=(1.009, 1), title=f"Dataset", frameon=False)

        #######################################
        figure_outpath = os.path.join(self.out_path, f'MainTextFigure4b_Supp.pdf')
        plt.savefig(figure_outpath)
        print(f'SAVED: {figure_outpath}')

        figure_outpath = os.path.join(self.out_path, f'MainTextFigure4b_Supp.png')
        plt.savefig(figure_outpath)
        print(f'SAVED: {figure_outpath}')

        figure_outpath = os.path.join(self.out_path, f'MainTextFigure4b_Supp.svg')
        plt.savefig(figure_outpath)
        print(f'SAVED: {figure_outpath}')
        #######################################  
    #################################################################################################################

    #################################################################################################################
    def plot_MainTextFigure4c_Supp(self,):
        """
        Figure 4c: 
            DnaK binding motif analysis
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Chaperone_Client_Associations/Dnak_simple_motif_scan/Plots/consolidated_Dnak_binding_data_EXP.csv       
        """
        buff_tag = {'C': 'cyto-serum', 'CD': '+DnaK', 'CG': '+GroEL'}
        #######################################
        ## Load Figure 2c data
        # Load post matched res_sasa
        for dataset in ['EXP', 'AF']:        
            inp = f'{self.slug_path}/Chaperone_Client_Associations/Dnak_simple_motif_scan/Plots/consolidated_Dnak_binding_data_{dataset}.csv'
            print(f'inp: {inp}')
            Figure_4c_df = pd.read_csv(inp)
            Figure_4c_df = Figure_4c_df[(Figure_4c_df['OnlyEnt'] == True) & (Figure_4c_df['motif'] == 'Schymkowitz')]

            Figure_4c_df['D_type'] = [ l.replace('D_', '') for l in Figure_4c_df['D_type'].values]
            Figure_4c_df['label'] = dataset
            print(f'Figure_4c_df:\n{Figure_4c_df}')

            Figure_4c_outfile_csv = os.path.join(self.out_path, f'MainTextFigure4c_{dataset}_Supp.csv')
            Figure_4c_df.to_csv(Figure_4c_outfile_csv)
            print(f'SAVED: {Figure_4c_outfile_csv}')
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
            for D_typei, D_type in enumerate(['All', 'Thread', 'Loop']):
                D_type_df = Figure_4c_df[Figure_4c_df['D_type'] == D_type]
                #print(D_type_df)

                x = D_type_df['spa'].values

                ## Plot 
                for esstag, label in [('Ess', 'Essential'), ('NonEss', 'Non-essential')]:
                    y = D_type_df[f'{esstag}_mean'].values
                    ylb = y - D_type_df[f'{esstag}_lower_ci'].values
                    yub = D_type_df[f'{esstag}_upper_ci'].values - y
                    y, ylb, yub = y*100, ylb*100, yub*100
                    yerr = [ylb, yub]
                    line1 = axes[0, D_typei].errorbar(x, y, yerr=yerr, fmt='o', capsize=3, markersize=3,  elinewidth=0.5,  markeredgewidth=0.5, label=label)
                    #axes[0, buffi].set_yscale('log')
                    #axes[0, D_typei].axhline(y=1.0, color='black', linestyle='--', linewidth=0.5)
                    if D_typei == 0:
                        lines += [line1]

                ## Plot the pvalues
                y = D_type_df['p_value'].values
                axes[1, D_typei].plot(x, y, marker='o', ls='none', markersize=3,  markeredgewidth=0.5)
                axes[1, D_typei].set_yscale('log')
                axes[1, D_typei].axhline(y=0.05, color='black', linestyle='--', linewidth=0.5)

                ## Plot the n values
                y1 = D_type_df['Ess_n'].values
                y2 = D_type_df['NonEss_n'].values
                axes[2, D_typei].plot(x, y1, marker='o', ls='none', markersize=3,  markeredgewidth=0.5, label=f'Essential')
                axes[2, D_typei].plot(x, y2, marker='o', ls='none', markersize=3,  markeredgewidth=0.5, label=f'Non-essential')
        
            #######################################

            ######################################
            ## adjust axes locations and labels
            panel_labels = np.asarray([['a', 'b', 'c'], ['d', 'e', 'f'], ['g', 'h', 'i']])
            panel_label_ys = {0:1, 1:0.66, 2:0.325}
            x0s = {0:0.075, 1:0.38, 2:0.685}
            y0s = {0:0.74, 1:0.4, 2:0.07}
            for coli, colname in enumerate(['All', 'Thread', 'Loop']):
                for rowi, rowname in enumerate(['Number of motif hits per 100 residues', 'p-value', 'Number of proteins in dataset']):

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
                    #if rowi == 0:
                    #    axes[rowi, coli].set_ylim(0.5, 2)

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

            # Adjust layout for clarity
            labels = ['Essential', 'Non-essential']
            fig.legend(lines, labels, loc='upper right', bbox_to_anchor=(1.009, 1), title="Dataset", frameon=False)

            #######################################
            figure_outpath = os.path.join(self.out_path, f'MainTextFigure4c_{dataset}_Supp.pdf')
            plt.savefig(figure_outpath)
            print(f'SAVED: {figure_outpath}')

            figure_outpath = os.path.join(self.out_path, f'MainTextFigure4c_{dataset}_Supp.png')
            plt.savefig(figure_outpath)
            print(f'SAVED: {figure_outpath}')

            figure_outpath = os.path.join(self.out_path, f'MainTextFigure4c_{dataset}_Supp.svg')
            plt.savefig(figure_outpath)
            print(f'SAVED: {figure_outpath}')
            #######################################  
    #################################################################################################################

    #################################################################################################################
    def plot_MainTextFigure4d_Supp(self,):
        """
        Figure 4d (row 2 column 2):
            Loop closing contacts enrichement matrix
            C buffers at spa50 and LiPMScov 50 in EXP dataset
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Sequence_Complexity_and_Discrimination/Compare_OR_with_permutation/EXP/C_50_p100000/FrequencyGeneratorOutput/OR_FDR_pvalues.csv
            ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Sequence_Complexity_and_Discrimination/Compare_OR_with_permutation/EXP/C_50_p100000/FrequencyGeneratorOutput/OR_GT.csv
        """
        #######################################
        ## Load Figure 4d data
        #######################################
        ## set up figure(s)
        fig_width_mm = 180  # Width in mm
        fig_height_mm = 170  # Height in mm
        fig, axes = plt.subplots(3, 2, figsize=(mm_to_inches(fig_width_mm), mm_to_inches(fig_height_mm)), gridspec_kw={'width_ratios': [1, 1]})
        for dataseti, dataset in enumerate(['EXP', 'AF']):  

            #######################################
            for buffi, buff in enumerate(['C', 'CD', 'CG']):
                inp = f'{self.slug_path}/Sequence_Complexity_and_Discrimination/Compare_OR_with_permutation/{dataset}/{buff}_50_p100000/FrequencyGeneratorOutput/OR_FDR_pvalues.csv'
                print(f'inp: {inp}')
                Figure_4d_FDR_df = pd.read_csv(inp)
                print(f'Figure_4d_FDR_df:\n{Figure_4d_FDR_df}')
                Figure_4d_FDR_df.set_index('AA', inplace=True)
                print(f'Figure_4d_FDR_df:\n{Figure_4d_FDR_df}')

                inp = f'{self.slug_path}/Sequence_Complexity_and_Discrimination/Compare_OR_with_permutation/{dataset}/{buff}_50_p100000/FrequencyGeneratorOutput/OR_GT.csv'
                print(f'inp: {inp}')
                Figure_4d_OR_df = pd.read_csv(inp)
                #Figure_4d_OR_df = Figure_4d_OR_df.replace([np.inf, -np.inf], np.nan).fillna(0)
                Figure_4d_OR_df.set_index('AA', inplace=True)
                Figure_4d_OR_df.at['C', 'C'] = 1
                print(f'Figure_4d_OR_df:\n{Figure_4d_OR_df}')
                #######################################

                Figure_4d_OR_df = Figure_4d_OR_df.replace([np.inf, -np.inf], np.nan).fillna(np.max(Figure_4d_OR_df.values[np.where(Figure_4d_OR_df.values != np.inf)]))
                #print(Figure_4d_OR_df)
      
                min_OR = min(Figure_4d_OR_df.values.flatten())
                max_OR = max(Figure_4d_OR_df.values.flatten())
                max_OR = 5
                
                # Create a TwoSlopeNorm for custom normalization
                norm = TwoSlopeNorm(vmin=min_OR, vcenter=1, vmax=max_OR)  # vcenter sets the whitest point

                # Replace `np.inf` with `max_OR` for visualization
                # Create a mask for the upper triangle
                visualization_df = Figure_4d_OR_df.replace(np.inf, max_OR)
                mask = np.triu(np.ones(visualization_df.shape), k=1)  # k=1 excludes the diagonal
                visualization_df = visualization_df.mask(mask == 1)
                print(visualization_df, visualization_df.shape)

                # Create a DataFrame for annotations, replacing `np.inf` with 'inf'
                annotations = Figure_4d_OR_df.round(decimals=1).replace(np.inf, "inf").astype(str)
                annotations.at['C', 'C'] = 'NA'
                print(annotations)

                #sns.heatmap(Figure_4d_OR_df, annot=Figure_4d_OR_df.round(decimals=1).astype(str), fmt="", cmap="coolwarm", norm=norm, vmin=min_OR, vmax=max_OR, linewidths=.5, cbar_kws={"shrink": .8, "aspect": 30, "label": "Odds Ratio"}, annot_kws={"fontsize": 5}, ax=ax)
                #sns.heatmap(visualization_df, annot=annotations, fmt="", cmap="coolwarm", norm=norm, vmin=min_OR, vmax=max_OR, linewidths=.5, cbar_kws={"shrink": .8, "aspect": 30, "label": "Odds Ratio"}, annot_kws={"fontsize": 5}, ax=ax)
                heatmap = sns.heatmap(visualization_df, annot=annotations, fmt="", cmap="coolwarm", norm=norm, vmin=min_OR, vmax=max_OR, linewidths=.5, cbar_kws={"shrink": .8, "aspect": 30, "label": "Odds Ratio"}, annot_kws={"fontsize": 5}, ax=axes[buffi, dataseti])
            
                # if there is a valid pvalue dataframe color those with values below 0.05 in a yellow highlight
                if isinstance(Figure_4d_FDR_df, pd.DataFrame):
                    #print('Pvalues found')
                    # Highlight cells based on df2 values
                    for i in range(Figure_4d_FDR_df.shape[0]):
                        for j in range(Figure_4d_FDR_df.shape[1]):
                            if Figure_4d_FDR_df.iloc[i, j] < 0.05:
                                if i == 9 and j == 2:
                                    axes[buffi, dataseti].add_patch(plt.Rectangle((j, i), 1, 1, fill=False, edgecolor='black', lw=0.75, ls='--')) # make T-Y dashed
                                else:
                                    axes[buffi, dataseti].add_patch(plt.Rectangle((j, i), 1, 1, fill=False, edgecolor='black', lw=0.75))

                # Set labels
                print(Figure_4d_OR_df.index, len(Figure_4d_OR_df.index))
                axes[buffi, dataseti].set_yticks(np.arange(len(Figure_4d_OR_df.columns)) + 0.5)
                axes[buffi, dataseti].set_yticklabels(Figure_4d_OR_df.index, rotation=0)

                print(Figure_4d_OR_df.columns, len(Figure_4d_OR_df.columns))
                axes[buffi, dataseti].set_xticks(np.arange(len(Figure_4d_OR_df.columns)) + 0.5)
                axes[buffi, dataseti].set_xticklabels(Figure_4d_OR_df.columns, rotation=0)

                axes[buffi, dataseti].set_xlabel('Amino Acids')
                axes[buffi, dataseti].set_ylabel('Amino Acids')

                axes[buffi, dataseti].tick_params(axis='y', labelsize=6, width=0.5)
                axes[buffi, dataseti].tick_params(axis='x', labelsize=6, width=0.5)


        ######################################
        ## adjust axes locations and labels
        panel_labels = np.asarray([['a', 'b'], ['c', 'd'], ['e', 'f']])
        panel_label_ys = {0:1, 1:0.675, 2:0.34}
        x0s = {0:0.05, 1:0.6}
        y0s = {0:0.725, 1:0.385, 2:0.05}
        for coli, colname in enumerate(['Crystal Structres', 'AlphaFold']):
            for rowi, rowname in enumerate(['Cyto-serum', '+DnaK', '+GroEL']):
                axs_position = axes[rowi, coli].get_position()
                x0 = x0s[coli]
                y0 = y0s[rowi]
                width0, height0 = axs_position.extents[2] - axs_position.extents[0], axs_position.extents[3] - axs_position.extents[1]
                print(axs_position, width0, height0, '->', x0, y0, 0.35, 0.25)
                axes[rowi, coli].set_position([x0, y0, 0.35, 0.25])  # [left, bottom, width, height]

                bbox_in_fig_coords = axes[rowi, coli].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
                fig.text(bbox_in_fig_coords.x0, panel_label_ys[rowi], panel_labels[rowi, coli], fontsize=8, fontweight='bold', va='top', ha='left')

                # Adjust the colorbar position
                colorbar = axes[rowi, coli].collections[0].colorbar
                colorbar.ax.set_position([x0 + 0.3, y0 + 0.075, 0.02, 0.15])  # [left, bottom, width, height]
                colorbar.ax.tick_params(width=0.5)
                # Get the existing tick positions
                tick_positions = colorbar.get_ticks()
                #tick_positions = np.linspace(0, 5, 6)  # Equal spacing for ticks
                custom_labels = [' 0', ' 1', ' 2', ' 3', ' 4', '≥5']
                colorbar.set_ticks(tick_positions)
                print(tick_positions, custom_labels)
                colorbar.set_ticklabels(custom_labels)

                axes[rowi, coli].text(0.5, 0.8, f'{rowname}\n{colname}', fontsize=7, va='top', ha='left', transform=axes[rowi, coli].transAxes)

        #plt.tight_layout()
        #######################################
        figure_outpath = os.path.join(self.out_path, f'MainTextFigure4d_Supp.pdf')
        plt.savefig(figure_outpath)
        print(f'SAVED: {figure_outpath}')

        figure_outpath = os.path.join(self.out_path, f'MainTextFigure4d_Supp.png')
        plt.savefig(figure_outpath)
        print(f'SAVED: {figure_outpath}')

        figure_outpath = os.path.join(self.out_path, f'MainTextFigure4d_Supp.svg')
        plt.savefig(figure_outpath)
        print(f'SAVED: {figure_outpath}')
        ####################################### 
 
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
    fig_width_mm = 180  # Width in mm
    fig_height_mm = 150  # Height in mm
    #custom_font_path = "/storage/group/epo2/default/ims86/miniconda3/envs/FtoF/fonts/Arial.ttf" # Path to your custom font
    #arial_font = fm.FontProperties(fname=custom_font_path) # Create a FontProperties object
    plt.rcParams['font.family'] = 'Arial'  # Change to your desired font, e.g., 'Times New Roman', 'DejaVu Sans', etc.
    plt.rcParams['font.size'] = 6  # Default font size
    plt.rcParams['pdf.fonttype'] = 42

    plotter.plot_MainTextFigure4a_Supp(chap='DnaK')
    plotter.plot_MainTextFigure4a_Supp(chap='GroEL')
    plotter.plot_MainTextFigure4b_Supp()
    plotter.plot_MainTextFigure4c_Supp()
    plotter.plot_MainTextFigure4d_Supp()


    # Adjust the colorbar position
    #colorbar = axs[3].collections[0].colorbar
    #colorbar.ax.set_position([0.935, 0.5, 0.02, 0.4])  # [left, bottom, width, height]
    #colorbar.ax.tick_params(width=0.5)
    # Get the existing tick positions
    #tick_positions = colorbar.get_ticks()
    #tick_positions = np.linspace(0, 5, 6)  # Equal spacing for ticks
    #custom_labels = [' 0', ' 1', ' 2', ' 3', ' 4', '≥5']
    #colorbar.set_ticks(tick_positions)
    #print(tick_positions, custom_labels)
    #colorbar.set_ticklabels(custom_labels)

    ########################################################

    print('NORMAL TERMINATION')

if __name__ == "__main__":
    main()

