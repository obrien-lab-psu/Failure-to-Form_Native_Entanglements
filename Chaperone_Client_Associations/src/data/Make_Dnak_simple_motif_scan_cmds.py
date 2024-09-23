import sys,os
import numpy as np
#python codes/Dnak_simple_motif_scan_FASTAprod.py 
#-o Dnak_simple_motif_scan_FASTA_FLiPPR_SPAwLiPMScov_prod/EXP/ 
#-f FASTA/ -m None -r inpfiles/rep_genes_pdb_chain.txt 
#-Eg ../../make_gene_lists/Gene_lists_FLiPPR_SPAwCOV_prod/EXP/EXP_0.6g_C_Rall_spa0_LiPMScov50_essential_genes.txt 
#-NEg ../../make_gene_lists/Gene_lists_FLiPPR_SPAwCOV_prod/EXP/EXP_0.6g_C_Rall_spa0_LiPMScov50_nonessential_genes.txt 
#-c ../../get_entanglements/entanglements_0.6g/clustered_unmapped_GE/ -b C -s 0 > logs/Dnak_simple_motif_scan_FASTAprod_EXP_C_spa0_LiPMScov50.log

path_to_slug = '../../../git_slugs/Failure-to-Form_Native_Entanglements_slug'
for set_type in ['EXP', 'AF']:
    for spa in np.arange(0,100,10):

        script = f'python src/data/Dnak_simple_motif_scan.py'
        outpath = f'-o {path_to_slug}/Chaperone_Client_Associations/Dnak_simple_motif_scan/{set_type}/'
        log = f'> {path_to_slug}/Chaperone_Client_Associations/Dnak_simple_motif_scan/logs/{set_type}_C_Rall_spa{spa}_LiPMScov50.log'
        Eg = f'-Eg {path_to_slug}/Make_Protein_Feature_Files/Gene_lists/{set_type}/{set_type}_0.6g_C_Rall_spa{spa}_LiPMScov50_essential_genes.txt'
        NEg = f'-NEg {path_to_slug}/Make_Protein_Feature_Files/Gene_lists/{set_type}/{set_type}_0.6g_C_Rall_spa{spa}_LiPMScov50_nonessential_genes.txt'
        misc = f'-f {path_to_slug}/Chaperone_Client_Associations/Dnak_simple_motif_scan/FASTA/ -r data/{set_type}_genes_pdb_chain.txt -b C -s {spa}'
        if set_type == 'EXP':
            ents = f'-c {path_to_slug}/Native_Entanglements_in_PDBs/Entanglements/Ecoli/EXP/mapped_NoSlipKNots_clustered_GE/'
        elif set_type == 'AF':
            ents = f'-c {path_to_slug}/Native_Entanglements_in_PDBs/Entanglements/Ecoli/AF/unmapped_HQ-AF_HQ-GE_clustered/'

        cmd = ' '.join([script, outpath, misc, Eg, NEg, ents])
        print(cmd)
  
        



