python ../../src/data/gaussian_entanglement.py --PDB PDBs/P0AD61-4YNG_C.pdb > logs/P0AD61-4YNG_C_gaussian_entanglement.log 
python ../../src/data/gaussian_entanglement.py --PDB PDBs/P00934-1VB3_A.pdb > logs/P00934-1VB3_A_gaussian_entanglement.log 
python ../../src/data/get_mapped_ent.py -e unmapped_GE/ -o mapped_GE/ -m PDBs/ > logs/get_mapped_ent.log
python ../../src/data/remove_EXP_slipknots.py -e mapped_GE/ -o mapped_NoSlipKNots_GE/ > logs/remove_EXP_slipknots.log
