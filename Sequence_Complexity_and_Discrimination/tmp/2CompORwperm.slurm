#!/bin/bash
#SBATCH -J c2CompORwperm
#SBATCH --partition=mgc-open
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 48:00:00
#SBATCH --account=epo2_nih
#SBATCH --mem=10G

cd /storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Sequence_Complexity_and_Discrimination
python src/data/compare_energetics_results_with_permutation.py -f ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Sequence_Complexity_and_Discrimination/EXP_FASTA/\*.fasta -o ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Sequence_Complexity_and_Discrimination/Compare_OR_with_permutation/EXP/CG_50_p100000/ -l ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Sequence_Complexity_and_Discrimination/logs/compare_energetics_CG_spa50_LiPMS50_p100000.log -c ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gen_proteome_features_EXP/contact_type2_lib/\*.csv -r ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gen_proteome_features_EXP/res_features_lib/\*.csv -g data/EXP_ecoli_uniprotIDs.txt -e ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Native_Entanglements_in_PDBs/Entanglements/Ecoli/EXP/mapped_NoSlipKNots_clustered_GE/\* -Eg ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/EXP/EXP_0.6g_CG_Rall_spa50_LiPMScov50_essential_ent_genes.txt -NEg ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/EXP/EXP_0.6g_CG_Rall_spa50_LiPMScov50_nonessential_ent_genes.txt -p 100000 --floor 0  

