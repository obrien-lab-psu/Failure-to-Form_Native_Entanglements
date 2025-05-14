#!/bin/bash
#SBATCH -J c16Q46856
#SBATCH --partition=mgc-open
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --account=epo2_nih
#SBATCH --mem=12G

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis
conda activate FtoF
python src/data/GQ.py --psf /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/Q46856_1OJ7_D/setup/Q46856_1OJ7_D_rebuilt_clean_ca.psf --cor /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/Q46856_1OJ7_D/setup/Q46856_1OJ7_D_rebuilt_clean_ca.cor --dcd /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/Q46856_1OJ7_D/Quenching/Q46856_1OJ7_D_t16_quench.dcd --sec_elements /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/Q46856_1OJ7_D/setup/secondary_struc_defs.txt  --outpath /storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis_setID3Full/Q46856_1OJ7_D/ --outname Q46856_1OJ7_D_t16 --topoly False

