
for set_type in ['EXP', 'AF']:
    for buff in ['C', 'Total']:
        script = "python src/data/CalcHydrophobicFraction.py"
        feature_files =  f"-f ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gen_proteome_features_{set_type}/res_features_lib/"
        outpath = f"-o ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Sequence_Complexity_and_Discrimination/CalcHydrophobicFraction/{set_type}/"
        log = f'-l ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Sequence_Complexity_and_Discrimination/CalcHydrophobicFraction/logs/{set_type}_{buff}.log'
        if buff == 'Total':
            Es = f'-Es ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/{set_type}/{set_type}_0.6g_{buff}_Rall_spa{buff}_LiPMScov{buff}_essential_ent_genes.txt'
            NEs = f'-NEs ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/{set_type}/{set_type}_0.6g_{buff}_Rall_spa{buff}_LiPMScov{buff}_nonessential_ent_genes.txt'
            tag = f'-t {set_type}_{buff}_Rall_spa{buff}_LiPMScov{buff}'
        else:
            Es = f'-Es ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/{set_type}/{set_type}_0.6g_{buff}_Rall_spa50_LiPMScov50_essential_ent_genes.txt'
            NEs = f'-NEs ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gene_lists/{set_type}/{set_type}_0.6g_{buff}_Rall_spa50_LiPMScov50_nonessential_ent_genes.txt'
            tag = f'-t {set_type}_{buff}_Rall_spa50_LiPMScov50'

        cmd = ' '.join([script, feature_files, outpath, log, Es, NEs, tag])
        print(cmd)
