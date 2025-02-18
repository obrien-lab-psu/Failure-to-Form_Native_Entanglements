import time, sys
import multiprocessing as mp
from scipy.stats import bootstrap
import logging
import argparse
import glob

"""
python src/data/Regression_Holdout.py 
-f ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gen_proteome_features_EXP/res_features_lib/\*.csv 
-o ./test_output/ 
-g ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/EXP/EXP_0.6g_C_Rall_spa50_LiPMScov50_essential_ent_genes.txt 
-t essential_ent_genes_270000 -b C -s 50 -c 50 -r "cut_C_Rall ~ AA + region" -l True -v region -sub_g test_output/holdouts_essential_ent_genes_270000.txt
"""

files = glob.glob('test_output/holdouts_essential_ent_genes_*.txt')

cmds = []
for f in files: 

    tag = f.split('holdouts_')[1].split('.txt')[0]
    #print(f, tag)
    cmd = f'python src/data/Regression_Holdout.py -f ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gen_proteome_features_EXP/res_features_lib/\*.csv -o ./test_output/ -g ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/EXP/EXP_0.6g_C_Rall_spa50_LiPMScov50_essential_ent_genes.txt -t {tag} -b C -s 50 -c 50 -r "cut_C_Rall ~ AA + region" -l True -v region -sub_g test_output/holdouts_{tag}.txt'
    print(cmd)

