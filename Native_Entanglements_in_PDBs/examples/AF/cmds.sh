python ../../src/data/gaussian_entanglement.py --PDB PDBs/P0AD61.pdb 
python ../../src/data/gaussian_entanglement.py --PDB PDBs/P00934.pdb 
python ../../src/data/get_HQ_AF_structures.py -e unmapped_GE/ -o unmapped_HQ-AF_HQ-GE/ -p PDBs/
python ../../src/data/clustering.py --prot_GE_file unmapped_HQ-AF_HQ-GE/P0AD61_GE.txt -o unmapped_HQ-AF_HQ-GE_clustered/ --organism Ecoli
python ../../src/data/clustering.py --prot_GE_file unmapped_HQ-AF_HQ-GE/P00934_GE.txt -o unmapped_HQ-AF_HQ-GE_clustered/ --organism Ecoli
