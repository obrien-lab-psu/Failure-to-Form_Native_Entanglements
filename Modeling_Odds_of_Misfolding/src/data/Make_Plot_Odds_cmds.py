import sys,os
import numpy as np
#python codes/Plot_Regression_results.py -f Regression_Rall_PROD/FLiPPR/EXP/whole_proteome/logRN_1_FDR-True_0.05/regression_results_region_all_genes_C\*.csv -o Plot_Regression_Rall_PROD/FLiPPR/EXP/whole_proteome/ -t all_genes -r region

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
        for gene_tag in gene_tags:

            script = f'python src/data/Plot_Odds_results.py'
            reg_res = f'-f {path_to_slug}/Modeling_Odds_of_Misfolding/Regressions/{set_type}/whole_proteome/regression_odds_holding-region_{gene_tag}_C\*_{timepoint}_\*.csv'
            outpath = f'-o {path_to_slug}/Modeling_Odds_of_Misfolding/Regressions/Plots/{set_type}/whole_proteome/'
            misc = f'-t {gene_tag}_{timepoint} -r region'
            #log = f'> {path_to_slug}/Modeling_Odds_of_Misfolding/Regressions/logs/{set_type}_whole_{buff}_{timepoint}_spa{spa}_LiPMScov{cov}_{gene_tag}_region.log'

            whole_cmd = ' '.join([script, reg_res, outpath, misc])
            print(whole_cmd)
            

            #PSM_cmd = f'python codes/Regression.py -f ../propensity_score_matching/PSM_PROD/{set_type}/FLiPPR_logRN_1_FDR-True_0.05/res_sasa/matched_residue_features_res_sasa_n1.csv -o Regression_PROD/FLiPPR/{set_type}/Propensity_score_matched/logRN_1_FDR-True_0.05/ -g ../make_gene_lists/Gene_lists_FLiPPR_SPAwCOV_prod/{set_type}/{set_type}_0.6g_{buff}_Rall_spa{spa}_LiPMScov{cov}_{gene_tag}.txt -t {gene_tag}_M1 -b {buff} -s {spa} -c {cov} -r "cut_{buff}_Rall ~ AA + region" -l False -v region > logs/{set_type}_PSM_FLiPPR_{buff}_spa{spa}_LiPMScov{cov}_{gene_tag}_region_M1.log'
            #print(PSM_cmd)


