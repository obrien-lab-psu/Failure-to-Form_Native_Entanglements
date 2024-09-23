#!/bin/bash
#SBATCH -J 
#SBATCH --partition=open
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 40:00:00
#SBATCH --account=open
#SBATCH --mem=6G

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Modeling_Odds_of_Misfolding
cmd
