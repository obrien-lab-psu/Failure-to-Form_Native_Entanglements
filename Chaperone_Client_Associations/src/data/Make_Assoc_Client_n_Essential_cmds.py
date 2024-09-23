import sys,os
import numpy as np
#python codes/ClientContingency_PROD.py 
#-Ag ../../make_gene_lists/Gene_lists_FLiPPR_SPAwCOV_prod/EXP/EXP_0.6g_C_Rall_spa0_LiPMScov50_all_genes.txt 
#-Eg inpfiles/Knockout_Essential.txt 
#-Entg ../../make_gene_lists/Gene_lists_FLiPPR_SPAwCOV_prod/EXP/EXP_0.6g_C_Rall_spa0_LiPMScov50_ent_genes.txt 
#-c inpfiles/Combined_ChapClient_DB.csv -l logs/EXP_0.6g_C_Rall_spa0_LiPMScov50_Client_Knockout.log 
#-o ClientContingency_PROD2/EXP/Knockout/ -t EXP_Knockout -b C -s 0 --LiPMScov 50

path_to_slug = '../../../git_slugs/Failure-to-Form_Native_Entanglements_slug'
for set_type in ['EXP', 'AF']:
    for spa in np.arange(0,100,10):

        script = f'python src/data/Assoc_Client_n_Essential.py'
        outpath = f'-o {path_to_slug}/Chaperone_Client_Associations/Assoc_Client_n_Essential/{set_type}/'
        log = f'-l {path_to_slug}/Chaperone_Client_Associations/Assoc_Client_n_Essential/logs/{set_type}_C_Rall_spa{spa}_LiPMScov50.log'
        Ag = f'-Ag {path_to_slug}/Make_Protein_Feature_Files/Gene_lists/{set_type}/{set_type}_0.6g_C_Rall_spa{spa}_LiPMScov50_all_genes.txt'
        Eg = f'-Eg data/Knockout_Essential.txt'
        Entg = f'-Entg {path_to_slug}/Make_Protein_Feature_Files/Gene_lists/{set_type}/{set_type}_0.6g_C_Rall_spa{spa}_LiPMScov50_ent_genes.txt '
        misc = f'-t {set_type}_Knockout -b C -s {spa} --LiPMScov 50 -c data/Combined_ChapClient_DB.csv'


        cmd = ' '.join([script, outpath, misc, log, Ag, Eg, Entg])
        print(cmd)
  
        



