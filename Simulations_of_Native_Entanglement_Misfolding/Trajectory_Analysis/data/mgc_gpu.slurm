#!/bin/bash
#SBATCH -J rebuild 
#SBATCH --partition=mgc-nih
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --account=epo2_nih
#SBATCH --mem=10G
#SBATCH --gres=gpu:1

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis
conda activate FtoF
python src/data/backmap.py --aa_pdb ../Rebuild_AllAtom_structures/data/post_rebuilt/Q46856_1OJ7_D_rebuilt.pdb --cg_pdb ../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/Q46856_1OJ7_D/Quenching/Q46856_1OJ7_D_t49_quench_finalframe26667.pdb -p 1
