import sys,os
import numpy as np

path_to_slug = '../../../git_slugs/Failure-to-Form_Native_Entanglements_slug'
gene_tags = ['all_genes',
'ent_genes',
'essential_ent_genes',
'essential_genes',
'nonessential_ent_genes',
'nonessential_genes']
#for cutT in [1,2,3]:
for set_type in ['EXP', 'AF']:
    for timepoint in ['Rall', 'R1min', 'R5min', 'R2hr']:
        for buff in ['C', 'CD', 'CG']:
            for spa in np.arange(0,100,10):
                for cov in np.arange(0,100,10):
                    for gene_tag in gene_tags:

                        script = f'python src/data/Regression.py'
                        res_feats = f'-f {path_to_slug}/Make_Protein_Feature_Files/Gen_proteome_features_{set_type}/res_features_lib/\*.csv'
                        outpath = f'-o {path_to_slug}/Modeling_Odds_of_Misfolding/Regressions/{set_type}/whole_proteome/'
                        gene_list = f'-g {path_to_slug}/Make_Protein_Feature_Files/Gene_lists/{set_type}/{set_type}_0.6g_{buff}_{timepoint}_spa{spa}_LiPMScov{cov}_{gene_tag}.txt'
                        misc = f'-t {gene_tag} -b {buff} -s {spa} -c {cov} -r "cut_{buff}_{timepoint} ~ AA + region" -l True -v region'
                        log = f'> {path_to_slug}/Modeling_Odds_of_Misfolding/Regressions/logs/{set_type}_whole_{buff}_{timepoint}_spa{spa}_LiPMScov{cov}_{gene_tag}_region.log'

                        whole_cmd = ' '.join([script, res_feats, outpath, gene_list, misc, log])
                        print(whole_cmd)
                        

                        #PSM_cmd = f'python codes/Regression.py -f ../propensity_score_matching/PSM_PROD/{set_type}/FLiPPR_logRN_1_FDR-True_0.05/res_sasa/matched_residue_features_res_sasa_n1.csv -o Regression_PROD/FLiPPR/{set_type}/Propensity_score_matched/logRN_1_FDR-True_0.05/ -g ../make_gene_lists/Gene_lists_FLiPPR_SPAwCOV_prod/{set_type}/{set_type}_0.6g_{buff}_Rall_spa{spa}_LiPMScov{cov}_{gene_tag}.txt -t {gene_tag}_M1 -b {buff} -s {spa} -c {cov} -r "cut_{buff}_Rall ~ AA + region" -l False -v region > logs/{set_type}_PSM_FLiPPR_{buff}_spa{spa}_LiPMScov{cov}_{gene_tag}_region_M1.log'
                        #print(PSM_cmd)


