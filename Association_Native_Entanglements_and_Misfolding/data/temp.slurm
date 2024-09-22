#!/bin/bash
#SBATCH -J 
#SBATCH --partition=open
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH --account=open
#SBATCH --mem=4G

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Association_Native_Entanglements_and_Misfolding
cmd
