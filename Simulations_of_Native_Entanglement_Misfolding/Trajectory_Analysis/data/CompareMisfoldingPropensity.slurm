#!/bin/bash
#SBATCH -J CMPx 
#SBATCH --partition=mgc-open
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --account=epo2_nih
#SBATCH --mem=10G

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis
conda activate FtoF
python src/data/CompareMisfoldingPropensity.py --outpath ../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CompareMisfoldingPropensity/ --candidates data/simulation_candidates_ids.csv --GQfiles data/GQfiles.txt --outname CompareMisfoldingPropensity
