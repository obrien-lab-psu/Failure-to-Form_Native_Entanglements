#!/bin/bash
#SBATCH -J c56Trend
#SBATCH --partition=mgc-open
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --account=epo2_nih
#SBATCH --mem=10G

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Sequence_Complexity_and_Discrimination
python src/data/Trend_analysis.py -o ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Sequence_Complexity_and_Discrimination/Trend_analysis/AF/ -r ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gen_proteome_features_AF/res_features_lib/ -g ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/AF/AF_0.6g_CD_Rall_spa60_LiPMScov50_ent_genes.txt -Eg ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/AF/AF_0.6g_CD_Rall_spa60_LiPMScov50_essential_ent_genes.txt -NEg ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/AF/AF_0.6g_CD_Rall_spa60_LiPMScov50_nonessential_ent_genes.txt -t AF_CD_spa60 -b CD -s 60 --LiPMScov 50 -c data/SigContact_df/AF/LoopformingContactsClassification.csv -l ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Sequence_Complexity_and_Discrimination/Trend_analysis/logs/AF_CD_60.log

