import sys,os
import numpy as np
#python codes/compare_hydropathy.py 
#-o compare_hydropathy_FLiPPR_SPAwCOV_PROD/EXP/C_0/ 
#-l logs/compare_hydropathy_EXP_C_0.log 
#-c ../proteome_features/Gen_proteome_features_EXP/FLiPPR_logRN_1_FDR-True_0.05/contact_type2_lib/\*.csv 
#-eg ../make_gene_lists/Gene_lists_FLiPPR_SPAwCOV_prod/EXP/EXP_0.6g_C_Rall_spa0_LiPMScov50_essential_ent_genes.txt 
#-neg ../make_gene_lists/Gene_lists_FLiPPR_SPAwCOV_prod/EXP/EXP_0.6g_C_Rall_spa0_LiPMScov50_nonessential_ent_genes.txt 
#-e ../get_entanglements/entanglements_0.6g/clustered_mapped_GE/\*



path_to_slug = '../../../git_slugs/Failure-to-Form_Native_Entanglements_slug'
for set_type in ['EXP', 'AF']:
    for buff in ['C', 'CG', 'CD']:
        for spa in np.arange(0, 100, 10):

            script = f'python src/data/compare_hydropathy.py'
            outpath = f'-o {path_to_slug}/Sequence_Complexity_and_Discrimination/hydropathy/{set_type}/'
            contacts = f'-c {path_to_slug}/Make_Protein_Feature_Files/Gen_proteome_features_{set_type}/contact_type2_lib/'        
            eg = f'-eg {path_to_slug}/Make_Protein_Feature_Files/Gene_lists/{set_type}/{set_type}_0.6g_{buff}_Rall_spa{spa}_LiPMScov50_essential_ent_genes.txt'
            neg = f'-neg {path_to_slug}/Make_Protein_Feature_Files/Gene_lists/{set_type}/{set_type}_0.6g_{buff}_Rall_spa{spa}_LiPMScov50_nonessential_ent_genes.txt'
            logs = f'-l {path_to_slug}/Sequence_Complexity_and_Discrimination/hydropathy/logs/{set_type}_{buff}_{spa}.log'
            tag  = f'-t {set_type}_{buff}_spa{spa}'
            if set_type == 'EXP':
                ents = f'-e {path_to_slug}/Native_Entanglements_in_PDBs/Entanglements/Ecoli/EXP/mapped_NoSlipKNots_clustered_GE/'
            elif set_type == 'AF':
                ents = f'-e {path_to_slug}/Native_Entanglements_in_PDBs/Entanglements/Ecoli/AF/unmapped_HQ-AF_HQ-GE_clustered/'
            cmd = ' '.join([script, outpath, contacts, eg, neg, tag, ents, logs])
            print(cmd)
        
        



