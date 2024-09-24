#!/bin/bash
#SBATCH -J c5CompE
#SBATCH --partition=mgc-open
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --account=epo2_nih
#SBATCH --mem=15G

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Sequence_Complexity_and_Discrimination
python src/data/compare_energetics_with_permutation_with_mp_full_prod.py -f ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Sequence_Complexity_and_Discrimination/AF_FASTA/\*.fasta -o ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Sequence_Complexity_and_Discrimination/compare_energetics_with_permutation_entOnly_floor0/AF/CD_50_p100000/ -l ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Sequence_Complexity_and_Discrimination/logs/compare_energetics_CD_spa50_LiPMS50_p100000.log -c ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gen_proteome_features_AF/contact_type2_lib/\*.csv -r ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gen_proteome_features_AF/res_features_lib/\*.csv -g data/AF_ecoli_uniprotIDs.txt -e ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Native_Entanglements_in_PDBs/Entanglements/Ecoli/AF/unmapped_HQ-AF_HQ-GE_clustered/\* -Eg ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/AF/AF_0.6g_CD_Rall_spa50_LiPMScov50_essential_ent_genes.txt -NEg ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/AF/AF_0.6g_CD_Rall_spa50_LiPMScov50_nonessential_ent_genes.txt -p 100000 --floor 0  

