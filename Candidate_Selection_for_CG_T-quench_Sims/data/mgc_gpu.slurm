#!/bin/bash
#SBATCH -J 
#SBATCH --partition=mgc-nih
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 48:00:00
#SBATCH --account=epo2_nih
#SBATCH --mem=20G
#SBATCH --gres=gpu:1

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Candidate_Selection_for_CG_T-quench_Sims
cmd
