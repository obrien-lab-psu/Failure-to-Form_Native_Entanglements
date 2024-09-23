import sys,os
import numpy as np
#python codes/Compare_ent_complexity_metrics.py -o Compare_ent_complexity_metrics_FLiPPR_SPAwLiPcov/EXP/ 
#-Eg ../make_gene_lists/Gene_lists_FLiPPR_SPAwCOV_prod/EXP/EXP_0.6g_C_Rall_spa0_LiPMScov50_essential_ent_genes.txt 
#-NEg ../make_gene_lists/Gene_lists_FLiPPR_SPAwCOV_prod/EXP/EXP_0.6g_C_Rall_spa0_LiPMScov50_nonessential_ent_genes.txt 
#-l logs/Compare_ent_complexity_metrics_FLiPPR_EXP_SPA0_LiPMScov50.log 
#-e ../proteome_features/Gen_proteome_features_EXP/FLiPPR_logRN_1_FDR-True_0.05/uent_features_lib/\* -p 100000 -b C -s 0

path_to_slug = '../../../git_slugs/Failure-to-Form_Native_Entanglements_slug'
for set_type in ['EXP', 'AF']:
    for buff in ['C', 'CD', 'CG']:
        for spa in np.arange(0,100,10):

            script = f'python src/data/Compare_ent_complexity_metrics.py'
            outpath = f'-o {path_to_slug}/Entanglement_Topological_Complexity_and_Discrimination/{set_type}/'
            log = f'-l {path_to_slug}/Entanglement_Topological_Complexity_and_Discrimination/logs/{set_type}_{buff}_Rall_spa{spa}_LiPMScov50.log'
            Eg = f'-Eg {path_to_slug}/Make_Protein_Feature_Files/Gene_lists/{set_type}/{set_type}_0.6g_{buff}_Rall_spa{spa}_LiPMScov50_essential_ent_genes.txt'
            NEg = f'-NEg {path_to_slug}/Make_Protein_Feature_Files/Gene_lists/{set_type}/{set_type}_0.6g_{buff}_Rall_spa{spa}_LiPMScov50_nonessential_ent_genes.txt'
            misc = f'-e {path_to_slug}/Make_Protein_Feature_Files/Gen_proteome_features_{set_type}/uent_features_lib/ -p 100000 -b {buff} -s {spa}'

            cmd = ' '.join([script, outpath, misc, Eg, NEg, log])
            print(cmd)
  
        



