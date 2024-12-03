#!/bin/bash
#SBATCH -J  
#SBATCH --partition=open
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 48:00:00
#SBATCH --account=open
#SBATCH --mem=10G

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis
conda activate MSM
cmd
