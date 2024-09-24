import sys,os
import numpy as np
#python codes/Plot_energetics_results_with_permutation.py 
#-En compare_energetics_with_permutation_entOnly_floor0_PROD/EXP/C_50_p10000/FrequencyGeneratorOutput/Enorm_ab_Ess_GT.csv 
#-EsdE compare_energetics_with_permutation_entOnly_floor0_PROD/EXP/C_50_p10000/FrequencyGeneratorOutput/deltaE_Ess_GT.csv 
#-EsdEp compare_energetics_with_permutation_entOnly_floor0_PROD/EXP/C_50_p10000/FrequencyGeneratorOutput/deltaE_Ess_FDR_pvalues.csv 
#-NEsdE compare_energetics_with_permutation_entOnly_floor0_PROD/EXP/C_50_p10000/FrequencyGeneratorOutput/deltaE_NonEss_GT.csv 
#-NEsdEp compare_energetics_with_permutation_entOnly_floor0_PROD/EXP/C_50_p10000/FrequencyGeneratorOutput/deltaE_NonEss_FDR_pvalues.csv 
#-dDE compare_energetics_with_permutation_entOnly_floor0_PROD/EXP/C_50_p10000/FrequencyGeneratorOutput/deltaDeltaE_GT.csv 
#-dDEp compare_energetics_with_permutation_entOnly_floor0_PROD/EXP/C_50_p10000/FrequencyGeneratorOutput/deltaDeltaE_FDR_pvalues.csv 
#-o Plot_compare_energetics_with_permutation_entOnly_SPAwLiPMScov_floor0_PROD/EXP/C_50_p10000/ -t C_50_p10000 


path_to_slug = '../../../git_slugs/Failure-to-Form_Native_Entanglements_slug'
for set_type in ['EXP', 'AF']:
    for tag in ['C_50_p100000', 'CD_50_p100000', 'CG_50_p100000', 'Total_p100000']:
    

        script = f'python src/data/Plot_energetics_results_with_permutation.py'
        outpath = f'-o {path_to_slug}/Sequence_Complexity_and_Discrimination/compare_energetics_with_permutation_entOnly_floor0/Plots/{set_type}/'
        En = f'-En {path_to_slug}/Sequence_Complexity_and_Discrimination/compare_energetics_with_permutation_entOnly_floor0/{set_type}/{tag}/FrequencyGeneratorOutput/Enorm_ab_Ess_GT.csv'
        EsdE = f'-EsdE {path_to_slug}/Sequence_Complexity_and_Discrimination/compare_energetics_with_permutation_entOnly_floor0/{set_type}/{tag}/FrequencyGeneratorOutput/deltaE_Ess_GT.csv'
        EsdEp = f'-EsdEp {path_to_slug}/Sequence_Complexity_and_Discrimination/compare_energetics_with_permutation_entOnly_floor0/{set_type}/{tag}/FrequencyGeneratorOutput/deltaE_Ess_FDR_pvalues.csv'
        NEsdE = f'-NEsdE {path_to_slug}/Sequence_Complexity_and_Discrimination/compare_energetics_with_permutation_entOnly_floor0/{set_type}/{tag}/FrequencyGeneratorOutput/deltaE_NonEss_GT.csv'
        NEsdEp = f'-NEsdEp {path_to_slug}/Sequence_Complexity_and_Discrimination/compare_energetics_with_permutation_entOnly_floor0/{set_type}/{tag}/FrequencyGeneratorOutput/deltaE_NonEss_FDR_pvalues.csv'
        dDE = f'-dDE {path_to_slug}/Sequence_Complexity_and_Discrimination/compare_energetics_with_permutation_entOnly_floor0/{set_type}/{tag}/FrequencyGeneratorOutput/deltaDeltaE_GT.csv'
        dDEp = f'-dDEp {path_to_slug}/Sequence_Complexity_and_Discrimination/compare_energetics_with_permutation_entOnly_floor0/{set_type}/{tag}/FrequencyGeneratorOutput/deltaDeltaE_FDR_pvalues.csv'
        plot_tag = f'-t {tag}_{set_type}'
        logs = f'> {path_to_slug}/Sequence_Complexity_and_Discrimination/compare_energetics_with_permutation_entOnly_floor0/Plots/logs/{tag}_{set_type}.log'
        cmd = ' '.join([script, outpath, En, EsdE, EsdEp, NEsdE, NEsdEp, dDE, dDEp, plot_tag, logs])
        print(cmd)
        
        



