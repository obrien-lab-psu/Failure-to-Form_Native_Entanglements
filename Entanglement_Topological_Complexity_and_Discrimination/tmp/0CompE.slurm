#!/bin/bash
#SBATCH -J c0CompE
#SBATCH --partition=open
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 40:00:00
#SBATCH --account=open
#SBATCH --mem=10G

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Entanglement_Topological_Complexity_and_Discrimination
python src/data/Compare_ent_complexity_metrics.py -o ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Entanglement_Topological_Complexity_and_Discrimination/EXP_022025/ -e ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gen_proteome_features_EXP/uent_features_lib/ -p 10000 -b C -s 50 -Eg ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/EXP/EXP_0.6g_C_Rall_spa50_LiPMScov50_essential_ent_genes.txt -NEg ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/EXP/EXP_0.6g_C_Rall_spa50_LiPMScov50_nonessential_ent_genes.txt -l ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Entanglement_Topological_Complexity_and_Discrimination/logs/EXP_022025_C_Rall_spa50_LiPMScov50.log

