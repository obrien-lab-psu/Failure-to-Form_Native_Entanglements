#!/bin/bash
#SBATCH -J 
#SBATCH --partition=mgc-nih
#SBATCH -N 1
#SBATCH -n 6
#SBATCH -t 72:00:00
#SBATCH --account=epo2_nih
#SBATCH --mem=30G
#SBATCH --gres=gpu:1

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis
conda activate FtoF
cmd
