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

def mm_to_inches(mm):
    return mm / 25.4

def main():
    """
    Creates Figure 1
    """

    # Create a master figure and axes
    fig_width_mm = 183  # Width in mm
    fig_height_mm = 100  # Height in mm
    #custom_font_path = "/storage/group/epo2/default/ims86/miniconda3/envs/FtoF/fonts/Arial.ttf" # Path to your custom font
    #arial_font = fm.FontProperties(fname=custom_font_path) # Create a FontProperties object
    #plt.rcParams['font.family'] = 'Arial'  # Change to your desired font, e.g., 'Times New Roman', 'DejaVu Sans', etc.
    plt.rcParams['font.size'] = 6  # Default font size
    plt.rcParams['pdf.fonttype'] = 42
    fig = plt.figure(figsize=(mm_to_inches(fig_width_mm), mm_to_inches(fig_height_mm)))
    fig, axs = plt.subplots(1, 1, 
                                figsize=(mm_to_inches(fig_width_mm), mm_to_inches(fig_height_mm)), 
                                dpi=600, constrained_layout=True) 


    ## load gene list
    gene_list = np.loadtxt('../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/EXP/EXP_0.6g_C_Rall_spa50_LiPMScov50_ent_genes.txt', dtype=str)
    print(f'gene_list: {gene_list} {gene_list.shape}')

    ## load feature files 
    #gene|pdb|chain|uniprot_length|essential|ent_present|pdb_resid|resname|AA|nearest_neighbors|num_nearest_neighbors|region|ent_idx|res_sasa|median_sasa|NC|crossing|mapped_resid|secondary_struct|SCOP_class|IDR|cut_str|cut_C_Rall|cut_CD_Rall|cut_CG_Rall|pdb_coverage|unresolved_IDR|buried
    ent_vectors = {}
    ranking = {'gene':[], 'ratio':[]}
    for gene in gene_list:
        feature_file = glob.glob(f'../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gen_proteome_features_EXP/res_features_lib/{gene}_*_resfeatures.csv')
        feature_df = pd.read_csv(feature_file[0], sep='|')
        feature_df = feature_df[['gene', 'mapped_resid', 'region', 'uniprot_length', 'cut_str']]
        feature_df = feature_df[~feature_df['mapped_resid'].isnull()]
        #print(feature_df)

        prot_size = feature_df['uniprot_length'].values[0]
        if prot_size >= 250 and prot_size <= 290:
            vec = np.zeros(prot_size)
            norm_length_vec = []
            for i, (rowi, row) in enumerate(feature_df.iterrows()):
                region = row['region']
                cut = isinstance(row['cut_str'], str)
                #print(gene, i, region, cut)

                if cut == True and region == 1:
                    tag = 2
                elif cut == True and region == 0:
                    tag = -2
                elif cut == False and region == 1:
                    tag = 1
                elif cut == False and region == 0:
                    tag = -1
                
                vec[i] = tag
                norm_length_vec += [row['mapped_resid']/prot_size]
            
            norm_length_vec = np.asarray(norm_length_vec)
            #print(f'vec: {vec} {vec.shape}')
            #print(f'norm_length_vec: {norm_length_vec} {norm_length_vec.shape}')
            fract_NE_cut = len(vec[np.where(vec == -2)]) / len(vec[np.where(vec == -1)])
            fract_E_cut = len(vec[np.where(vec == 2)]) / len(vec[np.where(vec == 1)])
            if fract_NE_cut != 0:
                ratio = fract_E_cut/fract_NE_cut
            else: 
                ratio = np.inf
            #print(f'fract_E_cut: {fract_E_cut} | fract_NE_cut: {fract_NE_cut} | ratio: {ratio}')
            ent_vectors[gene] = (vec, norm_length_vec)
            ranking['gene'] += [gene]
            ranking['ratio'] += [ratio]

    ###
    print(len(ent_vectors))
    ranking = pd.DataFrame(ranking)
    ranking = ranking.sort_values(by=['ratio'], ascending=True, ignore_index=True)
    print(ranking)

    # Create the plot
    plt.figure(figsize=(10, 6))

    color_dict = {-1:'lightblue', -2:'royalblue', 1:'pink', 2:'red'}

    for row_idx, row in ranking.iterrows():
        color_row, x_row = ent_vectors[row['gene']]
        #print(color_row, x_row)
        for col_idx, (color_code, x_val) in enumerate(zip(color_row, x_row)):
            color = color_dict.get(color_code, "#000000")  # Default to black if not found
            plt.fill_between([x_val - 0.0025, x_val + 0.0025], row_idx - 0.5, row_idx + 0.5, color=color, alpha=0.8)

    plt.xlabel("Normalized Length")
    plt.ylabel("Gene idx")
    plt.title("Shaded Rows with Assigned Colors")
    #plt.grid(True, linestyle="--", alpha=0.5)

    # Adjust y-ticks to match row indices
    plt.yticks(np.arange(len(ent_vectors)))

    # Save and close the plot
    plt.savefig('data/Plot_ProteinLines.svg', bbox_inches="tight")
    plt.close()

    print('NORMAL TERMINATION')

if __name__ == "__main__":
    main()

