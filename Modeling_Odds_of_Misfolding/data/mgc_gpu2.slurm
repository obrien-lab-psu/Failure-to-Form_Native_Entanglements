#!/bin/bash
#SBATCH -J PSMEXPT 
#SBATCH --partition=mgc-nih
#SBATCH -N 1
#SBATCH -n 30
#SBATCH -t 128:00:00
#SBATCH --account=epo2_nih
#SBATCH --mem=100G
#SBATCH --gres=gpu:1

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Modeling_Odds_of_Misfolding
conda activate FtoF
python src/data/PSM.py -f ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gen_proteome_features_EXP/res_features_lib/\*.csv -o ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Propensity_Score_Matching/EXP_test/ -m res_sasa -l ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Modeling_Odds_of_Misfolding/Propensity_Score_Matching/logs/PSM_EXP_test.log -n 1
