#!/bin/bash
#SBATCH -J  
#SBATCH --partition=open
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --account=open
#SBATCH --mem=30G

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Candidate_Selection_for_CG_T-quench_Sims
cmd
