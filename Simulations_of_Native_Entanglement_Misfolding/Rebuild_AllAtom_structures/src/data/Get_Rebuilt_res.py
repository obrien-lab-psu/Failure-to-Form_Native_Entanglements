#!/usr/bin/env python3
import requests, logging, os, sys
import time
import argparse
import pandas as pd
import numpy as np
import glob

pd.set_option('display.max_rows', 500)

############## MAIN #################
def main():
    
    script_name = f'Get_Rebuilt_res'
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--inpfiles", type=str, required=True, help="Path to candidates file")
    parser.add_argument("--mapping", type=str, required=True, help="Path to mapping files to check for insertions or deletions")
    parser.add_argument("--outpath", type=str, required=True, help="Path to output directory")
    args = parser.parse_args()

    files = glob.glob(os.path.join(args.inpfiles, '*'))
    dfs = []
    for f in files:
        dfs += [pd.read_csv(f)]
    dfs = pd.concat(dfs)

    dfs = dfs[(dfs['MISSING'] == True) | (dfs['MUTATION'] == True)]
    print(f'Rebuilt residues:\n{dfs.to_string()}')

    dfs.to_csv(args.outpath, index=False)
    print(f'SAVED: {args.outpath}')


    ## find mapping files and check for insertions or deletions
    Insertions = []
    Deletions = []
    for f in files:
        gene, pdb, chain = f.split('/')[-1].replace('.missing', '').split('_')
        #print(f, gene, pdb, chain)
        map_file = os.path.join(args.mapping, f'{gene}-{pdb}_{chain}_resid_mapping.txt')
        #print(f'map_file: {map_file}')

        # grep of Insertion and Deletion lines
        Ires = list(os.popen(f'grep Insertion {map_file}'))
        #print(f'grep Insertion {map_file}')
        #print(f'Ires: {Ires}')
        if len(Ires) != 0:
            Insertions += [f'{gene}-{pdb}_{chain} {I}' for I in Ires]


        Dres = list(os.popen(f'grep Deletion {map_file}'))
        #print(f'grep Deletion {map_file}')
        #print(f'Dres: {Dres}')
        if len(Dres) != 0:
            Deletions += [f'{gene}-{pdb}_{chain} {D}' for D in Dres]

    if len(Insertions) != 0:
        print(f'Insertions:\n{"".join(Insertions)}')
    else:
        print(f'No Insertions found in dataset')

    if len(Deletions) != 0:
        print(f'Deletions:\n{"".join(Deletions)}')
    else:
        print(f'No Deletions found in dataset')   
   

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
print(f'NORMAL TERMINATION: {end_time - start_time}')
logging.info(f'NORMAL TERMINATION: {end_time - start_time}')