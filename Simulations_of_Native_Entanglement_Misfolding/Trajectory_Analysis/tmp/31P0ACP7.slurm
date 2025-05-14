#!/bin/bash
#SBATCH -J c31P0ACP7
#SBATCH --partition=mgc-open
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --account=epo2_nih
#SBATCH --mem=12G

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis
conda activate FtoF
python src/data/GQ.py --psf /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/P0ACP7_1QPZ_A/setup/P0ACP7_1QPZ_A_rebuilt_clean_ca.psf --cor /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/P0ACP7_1QPZ_A/setup/P0ACP7_1QPZ_A_rebuilt_clean_ca.cor --dcd /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/P0ACP7_1QPZ_A/Quenching/P0ACP7_1QPZ_A_t31_quench.dcd --sec_elements /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/P0ACP7_1QPZ_A/setup/secondary_struc_defs.txt  --outpath /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis_setID3Full/P0ACP7_1QPZ_A/ --outname P0ACP7_1QPZ_A_t31 --topoly False

