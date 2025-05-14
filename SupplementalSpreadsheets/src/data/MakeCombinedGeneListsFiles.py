import pandas as pd
import numpy as np
import glob 
import os

"""
Collect the unique entanglement feature files for both experimental and alphafold data into a single parseable dataset

paths to feature files to collect
../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/EXP/EXP_0.6g_*_spa*_LiPMScov*_*_genes.txt
../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/AF/AF_0.6g_*_spa*_LiPMScov*_*_genes.txt
"""

### Check that the output directory is already made and if not, make it
if not os.path.exists("Combined_GeneList_Files"):
    os.makedirs("Combined_GeneList_Files")

##########################################################################################
### Get the experimental unique entanglement feature files and concatenate them into a single dataframe
files = glob.glob("../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/EXP/EXP_0.6g_*_spa50_LiPMScov50_*_genes.txt")

# Dictionary to hold data for each description
data = {}

for file in files:
    filename = file.split("/")[-1]
    buff = filename.split("_")[2]
    timepoint = filename.split("_")[3]
    description = '_'.join(filename.split("_")[6:-1])
    print(filename, buff, description)

    # Read the file contents
    with open(file, 'r') as f:
        lines = f.read().splitlines()
    
    # Initialize the dictionary for the description if not already present
    if description not in data:
        data[description] = {}
    
    # Add the lines to the corresponding buff in the description
    data[description]['-'.join([buff, timepoint])] = lines

# Create a Pandas Excel writer using XlsxWriter as the engine
with pd.ExcelWriter("Combined_GeneList_Files/EXP_combined_gene_lists_spa50_LiPMScov50.xlsx", engine='xlsxwriter') as writer:
    for description, buff_data in data.items():
        # Create a DataFrame from the dictionary
        df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in buff_data.items()]))
        # Write the DataFrame to an Excel sheet
        df.to_excel(writer, sheet_name=description, index=False)

print(f'SAVED: Combined_GeneList_Files/EXP_combined_gene_lists_spa50_LiPMScov50.xlsx')

##########################################################################################
### Get the experimental unique entanglement feature files and concatenate them into a single dataframe
files = glob.glob("../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/AF/AF_0.6g_*_spa50_LiPMScov50_*_genes.txt")

# Dictionary to hold data for each description
data = {}

for file in files:
    filename = file.split("/")[-1]
    buff = filename.split("_")[2]
    timepoint = filename.split("_")[3]
    description = '_'.join(filename.split("_")[6:-1])
    print(filename, buff, description)

    # Read the file contents
    with open(file, 'r') as f:
        lines = f.read().splitlines()
    
    # Initialize the dictionary for the description if not already present
    if description not in data:
        data[description] = {}
    
    # Add the lines to the corresponding buff in the description
    data[description]['-'.join([buff, timepoint])] = lines

# Create a Pandas Excel writer using XlsxWriter as the engine
with pd.ExcelWriter("Combined_GeneList_Files/AF_combined_gene_lists_spa50_LiPMScov50.xlsx", engine='xlsxwriter') as writer:
    for description, buff_data in data.items():
        # Create a DataFrame from the dictionary
        df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in buff_data.items()]))
        # Write the DataFrame to an Excel sheet
        df.to_excel(writer, sheet_name=description, index=False)

print(f'SAVED: Combined_GeneList_Files/AF_combined_gene_lists_spa50_LiPMScov50.xlsx')

print('NORMAL TERMINATION')
