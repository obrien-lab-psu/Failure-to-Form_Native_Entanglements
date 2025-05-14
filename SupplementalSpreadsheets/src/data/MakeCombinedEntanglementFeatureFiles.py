import pandas as pd
import numpy as np
import glob 
import os

"""
Collect the unique entanglement feature files for both experimental and alphafold data into a single parseable dataset

paths to feature files to collect
../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gen_proteome_features_EXP/uent_features_lib/
../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gen_proteome_features_AF/uent_features_lib/
"""

### Check that the output directory is already made and if not, make it
if not os.path.exists("Combined_Entanglement_FeatureFiles"):
    os.makedirs("Combined_Entanglement_FeatureFiles")

### Get the experimental unique entanglement feature files and concatenate them into a single dataframe
files = glob.glob("../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gen_proteome_features_EXP/uent_features_lib/*.csv")
df = []
for file in files:
    gene = file.split('/')[-1].split('_')[0]
    gene_df = pd.read_csv(file, sep='|')
    gene_df.insert(0, 'Gene', gene)  # Insert the 'Gene' column as the first column
    df.append(gene_df)
df = pd.concat(df)
print(df)
df.to_csv("Combined_Entanglement_FeatureFiles/EXP_combined_uent_features.csv", index=False)
df.to_excel("Combined_Entanglement_FeatureFiles/EXP_combined_uent_features.xlsx", index=False)

### Get the alphafold unique entanglement feature files and concatenate them into a single dataframe
files = glob.glob("../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gen_proteome_features_AF/uent_features_lib/*.csv")
df = []
for file in files:
    gene = file.split('/')[-1].split('_')[0]
    gene_df = pd.read_csv(file, sep='|')
    gene_df.insert(0, 'Gene', gene)  # Insert the 'Gene' column as the first column
    df.append(gene_df)
df = pd.concat(df)
print(df)
df.to_csv("Combined_Entanglement_FeatureFiles/AF_combined_uent_features.csv", index=False)
df.to_excel("Combined_Entanglement_FeatureFiles/AF_combined_uent_features.xlsx", index=False)

print("SAVED: Combined entanglement feature files for both experimental and alphafold data")
print("NORMAL TERMINATION")

