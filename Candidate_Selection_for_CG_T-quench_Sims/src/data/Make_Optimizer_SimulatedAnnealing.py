#python codes/Optimizer_SimulatedAnnealing_v4.1.py 
#-f ../proteome_features/Gen_proteome_features_EXP/FLiPPR_logRN_1_FDR-True_0.05/res_features_lib/\*.csv 
#-o Optimizer_SimulatedAnnealing_v4.1n2_NoiseTest_AllBuff/EXP/0/ 
#--ent_gene_list ../make_gene_lists/Gene_lists_FLiPPR_SPAwCOV_prod/EXP/EXP_0.6g_C_Rall_spa50_LiPMScov50_ent_genes.txt 

#--nonRefold_gene_list ../make_gene_lists/Gene_lists_FLiPPR_SPAwCOV_prod/EXP/EXP_0.6g_C_Rall_spa50_LiPMScov50_nonrefolded_genes.txt 
#-l logs/Optimizer_SimulatedAnnealing_v4.1n2_C_spa50_LiPMScov50_ent_genes_EXP0.log 
#--restart_path Optimizer_SimulatedAnnealing_v3.6n2_NoiseTest_AllBuff/EXP/0/ 
#-t ent_genes_M1 -b C -s 50 -c 50 -r "cut_C_Rall ~ AA + region" --random False -n 2 -C1 1 -C2 1 --steps 300 -beta 1000

import sys,os
import numpy as np

path_to_slug = '../../../git_slugs/Failure-to-Form_Native_Entanglements_slug'
for buff in ['C']:
    for random in ['False', 'True']:
        for set_type in ['EXP', 'AF']:
            for i in range(0, 20):

                #script = f'python src/data/Optimizer_SimulatedAnnealing.py'
                script = f'python src/data/Optimizer_SimulatedAnnealing.py'
                outpath = f'-o {path_to_slug}/Candidate_Selection_for_CG_T-quench_Sims/prod/Rand-{random}/{set_type}/{i}/'
                f = f'-f {path_to_slug}/Make_Protein_Feature_Files/Gen_proteome_features_{set_type}/res_features_lib/'       
                ent_gene_list = f'--ent_gene_list {path_to_slug}/Make_Protein_Feature_Files/Gene_lists/{set_type}/{set_type}_0.6g_{buff}_Rall_spa50_LiPMScov50_ent_genes.txt' 
                nonRefold_gene_list = f'--nonRefold_gene_list {path_to_slug}/Make_Protein_Feature_Files/Gene_lists/{set_type}/{set_type}_0.6g_{buff}_Rall_spa50_LiPMScov50_nonrefolded_genes.txt'

                logs = f'-l {path_to_slug}/Candidate_Selection_for_CG_T-quench_Sims/prod/logs/{set_type}_{buff}_{i}_Rand-{random}_linearT-False.log'
                #tag  = f'-t {set_type}_{buff}_traj{i} -b {buff} -s 50 -c 50 -r "cut_{buff}_Rall ~ AA + region" --random {random} -n 5 -beta 0.05 -C1 1 -C2 2.5 --steps 100000'
                tag  = f'-t {set_type}_{buff}_traj{i} -b {buff} -s 50 -c 50 -r "cut_{buff}_Rall ~ AA + region" --random {random} -n 4 -beta 0.05 -C1 1 -C2 2.5 -C3 0 --steps 200000 -linearT False'

                cmd = ' '.join([script, outpath, f, ent_gene_list, nonRefold_gene_list, tag, logs])
                print(cmd)



quit()   
for set_type in ['EXP', 'AF']:

        for i in range(0,10):

            script = f'python src/data/Optimizer_SimulatedAnnealing.py'
            outpath = f'-o {path_to_slug}/Candidate_Selection_for_CG_T-quench_Sims/Rand/{set_type}/{i}/'
            f = f'-f {path_to_slug}/Make_Protein_Feature_Files/Gen_proteome_features_{set_type}/res_features_lib/'       
            ent_gene_list = f'--ent_gene_list {path_to_slug}/Make_Protein_Feature_Files/Gene_lists/{set_type}/{set_type}_0.6g_{buff}_Rall_spa50_LiPMScov50_ent_genes.txt' 
            nonRefold_gene_list = f'--nonRefold_gene_list {path_to_slug}/Make_Protein_Feature_Files/Gene_lists/{set_type}/{set_type}_0.6g_{buff}_Rall_spa50_LiPMScov50_nonrefolded_genes.txt'

            logs = f'-l {path_to_slug}/Candidate_Selection_for_CG_T-quench_Sims/logs/{set_type}_{buff}_{i}_Rand.log'
            tag  = f'-t {set_type}_{buff}_traj{i}_Rand -b {buff} -s 50 -c 50 -r "cut_{buff}_Rall ~ AA + region" --random True -n 5 -beta 0.05 -C1 1 -C2 1 --steps 100000'

            cmd = ' '.join([script, outpath, f, ent_gene_list, nonRefold_gene_list, tag, logs])
            print(cmd)


