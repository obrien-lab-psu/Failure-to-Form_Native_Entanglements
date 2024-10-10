#!/bin/bash
#SBATCH -J test 
#SBATCH --partition=mgc-nih
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH -t 24:00:00
#SBATCH --account=epo2_nih

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Simulations_of_Native_Entanglement_Misfolding/Coarse_Graining_Model
conda activate FtoF
python src/data/opt_nscal.py -i /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Simulations_of_Native_Entanglement_Misfolding/Rebuild_AllAtom_structures/data/post_rebuilt/P0A6B4_4WR3_A_rebuilt.pdb -d data/domain_files/P0A6B4_4WR3_A.dat -o ../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Coarse_Graining_Model/P0A6B4_4WR3_A/ -p -1 
