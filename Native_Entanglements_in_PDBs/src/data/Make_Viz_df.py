import sys,os
import pandas as pd
import argparse
import glob
import ast
import numpy as np
import shutil

#pd.set_option('display.max_rows', 4000)
if len(sys.argv) != 5:
    print('[1] path to feature files')
    print('[2] path to outfile')
    print('[3] path to PDBs')
    print('[4] path to clustered GE files')
    quit()

featfiles = glob.glob(os.path.join(sys.argv[1], '*'))
outfile = sys.argv[2]
pdbfiles = glob.glob(os.path.join(sys.argv[3], '*'))
centfiles = glob.glob(os.path.join(sys.argv[4], '*'))
print(f'Number of feature files: {len(featfiles)}')
print(f'Number of clustered: {len(centfiles)}')

outdf = {'UniprotID':[], 
        'PDB':[],
        'Chain':[],
        'ENT-ID':[],
        'i':[],
        'i-resname':[],
        'j':[],
        'j-resname':[],
        'Gn':[],
        'Gc':[],
        'crossing_res':[],
        'crossing_res_w_buff':[],
        'surrounding_res':[],
        'C_cut':[],
        'CD_cut':[],
        'CG_cut':[],
        'C_Ecut':[],
        'CD_Ecut':[],
        'CG_Ecut':[],
        'entR_size':[],
        'nonentR_size':[],
        'entR_cov':[]}


keys = ['gene', 'pdb', 'chain', 'essential', 'ent_present', 'pdb_resid', 'resname', 'AA', 'region', 'ent_idx', 'NC', 'crossing', 'mapped_resid', 'cut_C_Rall', 'cut_CD_Rall', 'cut_CG_Rall']
for fi, featf in enumerate(featfiles):
    gene = featf.split('/')[-1].split('_')[0]
    #if gene != 'P0A799':
    #    continue
    print('\n',fi, featf, gene)
    df = pd.read_csv(featf, sep='|', dtype={'pdb':str})
    df = df[keys]
    print(df)
    df_all = df.copy()
    if True not in df['ent_present'].values:
        continue
    df = df[~df['ent_idx'].isnull()]
    df = df[df['region'] == 1]
    entR_size = len(df)
    nonentR_size = len(df_all) - entR_size
    entR_cov = entR_size/len(df_all)

    if len(df) == 0:
        print(f'ERROR: There should be an entanglement here for {gene}!')
        continue
    print(np.unique(df['pdb'].values))
    pdb = df['pdb'].unique()[0]
    chain = df['chain'].unique()[0]
    print(gene, pdb, chain)

    # get unique entIDs
    ENT_IDs = []
    for x in df['ent_idx'].values:
        ENT_IDs += ast.literal_eval(x)
    ENT_IDs = np.unique(ENT_IDs)
    #print(f'ENT_IDs: {ENT_IDs}')

    # get clustered ent file if it exists
    cfile = [f for f in centfiles if gene in f]
    if len(cfile) == 0:
        print('No clustered ent')
        continue
    else:
        uent_df = pd.read_csv(cfile[0], sep='|')
        print(uent_df)
        # copy pdb file
        #pfile = [p for p in pdbfiles if f'{gene}-{pdb}_{chain}' in p] ## for experimental PDBs
        pfile = [p for p in pdbfiles if f'{gene}' in p] ## for AF
        if len(pfile) == 0:
            raise ValueError("PDB file missing")
        pfile = pfile[0]
        shutil.copy(pfile, f'/storage/group/epo2/default/ims86/git_repos/viz_AF_entanglements/PDB/{gene}.pdb')
        print(f'Copied {pfile} -> /storage/group/epo2/default/ims86/git_repos/viz_AF_entanglements/PDB/{gene}.pdb')

        for erowi, erow in uent_df.iterrows():
            ijr = ast.literal_eval(erow['ijr'])
            i = ijr[0]
            j = ijr[1]
            crossings = []
            for c in ijr[2:]:
                if c.count('-') == 2:
                    c = c.replace('--', '-')
                crossings += [abs(ast.literal_eval(c))]
            #crossings = [ast.literal_eval(c) for c in ijr[2:]]
            Gn = erow['gn']
            Gc = erow['gc']
            print('\n', erowi, ijr, i, j, crossings, Gn, Gc)

            core_res = [i,j] + crossings
            erow_df = df[df['pdb_resid'].isin(core_res)]
            print(erow_df)
            if len(erow_df) == 0:
                print(f'No entangled information found meaning this was a non mapped entanglement')
                continue
            
            if len(erow_df) != len(core_res):
                print(f'No entangled information found meaning this was a non mapped entanglement')
                continue

            # check if any of the key entanglement info is unmapped
            mapped_erow_df = erow_df[erow_df['mapped_resid'].notna()]
            print(mapped_erow_df)
            if len(mapped_erow_df) != len(erow_df):
                print(f'This entanglement contains unmapped information, skipping')
                continue

            # get the resname for native contacts
            iresname = df[df['pdb_resid'] == i]['AA'].values[0]
            jresname = df[df['pdb_resid'] == j]['AA'].values[0]

            # get surrounding res
            surrounding = []
            for rowi, row in df.iterrows():
                row_ENTID = ast.literal_eval(row['ent_idx'])
                if str(erowi + 1) in row_ENTID:
                    row_mapped_reisd = row['mapped_resid']
                    
                    if isinstance(row_mapped_reisd, int) or isinstance(row_mapped_reisd, float):
                        row_pdb_reisd = row['pdb_resid']
                        row_pdb_reisd = int(row_pdb_reisd)
                        surrounding += [row_pdb_reisd]
       
            print(f'surrounding: {surrounding}')

            # get and all cuts
            C_cut = []
            CD_cut = []
            CG_cut = []
            for rowi, row in df_all.iterrows():
                row_mapped_reisd = row['mapped_resid']
                
                if isinstance(row_mapped_reisd, int) or isinstance(row_mapped_reisd, float):
                    row_pdb_reisd = row['pdb_resid']
                    row_pdb_reisd = int(row_pdb_reisd)
                    row_cut_C_Rall = row['cut_C_Rall']
                    row_cut_CD_Rall = row['cut_CD_Rall']
                    row_cut_CG_Rall = row['cut_CG_Rall']

                    if row_cut_C_Rall:
                        C_cut += [row_pdb_reisd]
                    if row_cut_CD_Rall:
                        CD_cut += [row_pdb_reisd]
                    if row_cut_CG_Rall:
                        CG_cut += [row_pdb_reisd]

            # get and entangled region cuts
            C_Ecut = []
            CD_Ecut = []
            CG_Ecut = []
            for rowi, row in df.iterrows():
                row_mapped_reisd = row['mapped_resid']
                
                if isinstance(row_mapped_reisd, int) or isinstance(row_mapped_reisd, float):
                    row_pdb_reisd = row['pdb_resid']
                    row_pdb_reisd = int(row_pdb_reisd)
                    row_cut_C_Rall = row['cut_C_Rall']
                    row_cut_CD_Rall = row['cut_CD_Rall']
                    row_cut_CG_Rall = row['cut_CG_Rall']

                    if row_cut_C_Rall:
                        C_Ecut += [row_pdb_reisd]
                    if row_cut_CD_Rall:
                        CD_Ecut += [row_pdb_reisd]
                    if row_cut_CG_Rall:
                        CG_Ecut += [row_pdb_reisd]

            # get crossing, surrounding, and cut site strings
            crossings_str = ';'.join([str(c) for c in crossings])
            crossings_w_buff_str = np.hstack([np.arange(c-5, c+6) for c in crossings])
            crossings_w_buff_str = ';'.join([str(c) for c in crossings_w_buff_str])

            surrounding_str = ';'.join([str(s) for s in surrounding])


            # get those cuts not in entangled regins
            C_cut_str = ';'.join([str(c) for c in C_cut if c not in C_Ecut])
            CD_cut_str = ';'.join([str(c) for c in CD_cut if c not in CD_Ecut])
            CG_cut_str = ';'.join([str(c) for c in CG_cut if c not in CG_Ecut])
            print(f'C_cut_str: {C_cut_str}')
            print(f'CD_cut_str: {CD_cut_str}')
            print(f'CG_cut_str: {CG_cut_str}')

            C_Ecut_str = ';'.join([str(c) for c in C_Ecut])
            CD_Ecut_str = ';'.join([str(c) for c in CD_Ecut])
            CG_Ecut_str = ';'.join([str(c) for c in CG_Ecut])
            print(f'C_Ecut_str: {C_Ecut_str}')
            print(f'CD_Ecut_str: {CD_Ecut_str}')
            print(f'CG_Ecut_str: {CG_Ecut_str}')

            outdf['UniprotID'] += [gene]            
            outdf['PDB'] += [pdb]            
            outdf['Chain'] += [chain]            
            outdf['ENT-ID'] += [erowi + 1]            
            outdf['i'] += [i]            
            outdf['i-resname'] += [iresname]            
            outdf['j'] += [j]            
            outdf['j-resname'] += [jresname]            
            outdf['Gn'] += [Gn]            
            outdf['Gc'] += [Gc]            
            outdf['crossing_res'] += [crossings_str]            
            outdf['crossing_res_w_buff'] += [crossings_w_buff_str]            
            outdf['surrounding_res'] += [surrounding_str]            
            outdf['C_cut'] += [C_cut_str]            
            outdf['CD_cut'] += [CD_cut_str]            
            outdf['CG_cut'] += [CG_cut_str] 
            outdf['C_Ecut'] += [C_Ecut_str]            
            outdf['CD_Ecut'] += [CD_Ecut_str]            
            outdf['CG_Ecut'] += [CG_Ecut_str]     
            outdf['entR_size'] += [entR_size]  
            outdf['nonentR_size'] += [nonentR_size]     
            outdf['entR_cov'] += [entR_cov]              
            

outdf = pd.DataFrame(outdf)
print(outdf)
outdf.to_csv(outfile, sep='|', index=False)
print(f'SAVED: {outfile}')
