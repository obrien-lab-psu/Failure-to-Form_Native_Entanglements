python ../../src/data/gaussian_entanglement.py --PDB PDBs/P0AD61-4YNG_C.pdb 
python ../../src/data/gaussian_entanglement.py --PDB PDBs/P00934-1VB3_A.pdb 
python ../../src/data/get_mapped_ent.py -e unmapped_GE/ -o mapped_GE/ -m PDBs/ 
python ../../src/data/remove_EXP_slipknots.py -e mapped_GE/ -o mapped_NoSlipKNots_GE/ 
