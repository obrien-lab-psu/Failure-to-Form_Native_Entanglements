#!/bin/bash
#SBATCH -J 
#SBATCH --partition=mgc-open
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH --account=epo2_nih
#SBATCH --mem=25G

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Make_Protein_Feature_Files
cmd
