#!/bin/bash
#SBATCH -J 
#SBATCH --partition=open
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 40:00:00
#SBATCH --account=open
#SBATCH --mem=10G

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Entanglement_Topological_Complexity_and_Discrimination
cmd
