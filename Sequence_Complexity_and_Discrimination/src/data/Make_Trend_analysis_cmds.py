import sys,os
import numpy as np
#python codes/Trend_analysis.py 
#-g ../make_gene_lists/Gene_lists_FLiPPR_SPAwCOV_prod/EXP/EXP_0.6g_C_Rall_spa0_LiPMScov50_ent_genes.txt 
#-Eg ../make_gene_lists/Gene_lists_FLiPPR_SPAwCOV_prod/EXP/EXP_0.6g_C_Rall_spa0_LiPMScov50_essential_ent_genes.txt 
#-NEg ../make_gene_lists/Gene_lists_FLiPPR_SPAwCOV_prod/EXP/EXP_0.6g_C_Rall_spa0_LiPMScov50_nonessential_ent_genes.txt 
#-r ../proteome_features/Gen_proteome_features_EXP/FLiPPR_logRN_1_FDR-True_0.05/res_features_lib/ 
#-c inpfiles/LoopformingContactsClassification.csv 
#-o Trend_analysis/EXP/C_spa0_LiPMScov50/ 
#-t EXP -b C -s 0 --LiPMScov 50 -l logs/Trend_analysis_EXP_0.6g_C_Rall_spa0_LiPMScov50.log

path_to_slug = '../../../git_slugs/Failure-to-Form_Native_Entanglements_slug'
for set_type in ['EXP', 'AF']:
    for buff in ['C', 'CG', 'CD']:
        for spa in np.arange(0, 100, 10):

            script = f'python src/data/Trend_analysis.py'
            outpath = f'-o {path_to_slug}/Sequence_Complexity_and_Discrimination/Trend_analysis/{set_type}/'

            res_feats = f'-r {path_to_slug}/Make_Protein_Feature_Files/Gen_proteome_features_{set_type}/res_features_lib/'       
            g = f'-g {path_to_slug}/Make_Protein_Feature_Files/Gene_lists/{set_type}/{set_type}_0.6g_{buff}_Rall_spa{spa}_LiPMScov50_ent_genes.txt' 
            Eg = f'-Eg {path_to_slug}/Make_Protein_Feature_Files/Gene_lists/{set_type}/{set_type}_0.6g_{buff}_Rall_spa{spa}_LiPMScov50_essential_ent_genes.txt'
            NEg = f'-NEg {path_to_slug}/Make_Protein_Feature_Files/Gene_lists/{set_type}/{set_type}_0.6g_{buff}_Rall_spa{spa}_LiPMScov50_nonessential_ent_genes.txt'
            logs = f'-l {path_to_slug}/Sequence_Complexity_and_Discrimination/Trend_analysis/logs/{set_type}_{buff}_{spa}.log'
            tag  = f'-t {set_type}_{buff}_spa{spa} -b {buff} -s {spa} --LiPMScov 50'
            contact_df = f'-c data/SigContact_df/{set_type}/LoopformingContactsClassification.csv'
            cmd = ' '.join([script, outpath, res_feats, g, Eg, NEg, tag, contact_df, logs])
            print(cmd)
        
        



