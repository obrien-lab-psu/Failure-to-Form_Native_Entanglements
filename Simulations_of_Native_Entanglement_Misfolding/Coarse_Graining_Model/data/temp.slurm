#!/bin/bash
#SBATCH -J 
#SBATCH --partition=mgc-nih
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH -t 100:00:00
#SBATCH --account=epo2_nih

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Simulations_of_Native_Entanglement_Misfolding/Coarse_Graining_Model
conda activate FtoF
cmd
