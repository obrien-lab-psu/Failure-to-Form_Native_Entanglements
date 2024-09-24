#!/bin/bash
#SBATCH -J c15Trend
#SBATCH --partition=mgc-open
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --account=epo2_nih
#SBATCH --mem=10G

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Sequence_Complexity_and_Discrimination
python src/data/Trend_analysis.py -o ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Sequence_Complexity_and_Discrimination/Trend_analysis/EXP/ -r ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gen_proteome_features_EXP/res_features_lib/ -g ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/EXP/EXP_0.6g_CG_Rall_spa50_LiPMScov50_ent_genes.txt -Eg ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/EXP/EXP_0.6g_CG_Rall_spa50_LiPMScov50_essential_ent_genes.txt -NEg ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/EXP/EXP_0.6g_CG_Rall_spa50_LiPMScov50_nonessential_ent_genes.txt -t EXP_CG_spa50 -b CG -s 50 --LiPMScov 50 -c data/SigContact_df/EXP/LoopformingContactsClassification.csv -l ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Sequence_Complexity_and_Discrimination/Trend_analysis/logs/EXP_CG_50.log

