#!/bin/bash
#SBATCH -J 
#SBATCH --partition=mgc-open
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 24:00:00
#SBATCH --account=epo2_nih
#SBATCH --mem=12G

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis
conda activate FtoF
cmd
