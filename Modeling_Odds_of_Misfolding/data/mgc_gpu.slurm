#!/bin/bash
#SBATCH -J PSMAF 
#SBATCH --partition=mgc-nih
#SBATCH -N 1
#SBATCH -n 30
#SBATCH -t 128:00:00
#SBATCH --account=epo2_nih
#SBATCH --mem=100G
#SBATCH --gres=gpu:1

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Modeling_Odds_of_Misfolding
python src/data/PSM.py -f ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gen_proteome_features_AF/res_features_lib/\*.csv -o ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Propensity_Score_Matching/AF/ -m res_sasa -l ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Propensity_Score_Matching/logs/PSM_AF.log -n 1
