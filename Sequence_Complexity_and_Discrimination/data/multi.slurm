#!/bin/bash
#SBATCH -J  
#SBATCH --partition=mgc-open
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 48:00:00
#SBATCH --account=epo2_nih
#SBATCH --mem=10G

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Sequence_Complexity_and_Discrimination
cmd
