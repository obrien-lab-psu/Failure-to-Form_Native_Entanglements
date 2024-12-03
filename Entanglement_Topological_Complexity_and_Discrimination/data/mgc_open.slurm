#!/bin/bash
#SBATCH -J  
#SBATCH --partition=mgc-open
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 72:00:00
#SBATCH --account=epo2_nih
#SBATCH --mem=10G

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Entanglement_Topological_Complexity_and_Discrimination
cmd
