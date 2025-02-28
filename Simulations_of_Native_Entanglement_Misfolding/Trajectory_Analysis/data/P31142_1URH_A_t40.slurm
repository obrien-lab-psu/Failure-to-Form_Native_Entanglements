#!/bin/bash
#SBATCH -J c28Reb
#SBATCH --partition=mgc-open
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 48:00:00
#SBATCH --account=epo2_nih
#SBATCH --mem=20G

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis
conda activate FtoF
python src/data/backmap.py --aa_pdb ../Rebuild_AllAtom_structures/data/post_rebuilt/P31142_1URH_A_rebuilt.pdb --cg_pdb ../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/P31142_1URH_A/Quenching/P31142_1URH_A_t40_quench_finalframe26667.pdb -p 1

