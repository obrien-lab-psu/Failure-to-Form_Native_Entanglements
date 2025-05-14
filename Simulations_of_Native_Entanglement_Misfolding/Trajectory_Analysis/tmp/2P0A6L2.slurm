#!/bin/bash
#SBATCH -J c2P0A6L2
#SBATCH --partition=mgc-open
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 72:00:00
#SBATCH --account=epo2_nih
#SBATCH --mem=20G

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis
conda activate FtoF
python src/data/GQ.py --psf /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/P0A6L2_1DHP_A/setup/P0A6L2_1DHP_A_rebuilt_clean_ca.psf --cor /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/P0A6L2_1DHP_A/setup/P0A6L2_1DHP_A_rebuilt_clean_ca.cor --dcd /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/P0A6L2_1DHP_A/Quenching/P0A6L2_1DHP_A_t42_quench.dcd --sec_elements /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/P0A6L2_1DHP_A/setup/secondary_struc_defs.txt  --outpath /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis_setID3Full/P0A6L2_1DHP_A/ --outname P0A6L2_1DHP_A_t42 --topoly False

