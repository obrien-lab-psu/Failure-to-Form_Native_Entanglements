import sys,os
import glob
import numpy as np
import pandas as pd
import argparse
pd.set_option('display.max_rows', 1000)  # Display 10 rows


###Load user defined arguments
parser = argparse.ArgumentParser(description="Process user specified arguments")
parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory")
parser.add_argument("-e", "--express_file", type=str, required=True, help='path to expression control file or Total')
parser.add_argument("-l", "--lipms_covs_file", type=str, required=True, help='path to lipms coverage file')
args = parser.parse_args()

outpath = args.outpath
express_file = args.express_file
lipms_covs_file = args.lipms_covs_file


## Load the file containg the list of knot genes
EXP_knots = pd.read_csv('data/KnotProt2.0.txt', sep=';')
EXP_knots = [f'{pdb.upper()}{chain}' for pdb, chain in zip(EXP_knots['pdbid'], EXP_knots['chain'])]
print(f'EXP_knots:\n{EXP_knots}')

AF_knots = pd.read_csv('data/AlphaKnotProt.txt', sep=';')
AF_knots = AF_knots['gene'].values
print(f'AF_knots:\n{AF_knots}')


## make outdirectory if it doesnt exists yet
if not os.path.exists(f'{outpath}'):
    os.makedirs(f'{outpath}')
    print(f'Made output directories {outpath}')


## Loading expression level SPA data
print(f'LOADING express_file ...')
print(express_file)
SPA_df = pd.read_pickle(express_file)
print(SPA_df.keys())


## Loading coverage level data
print(f'LOADING lipms_covs_file ...')
print(lipms_covs_file)
COV_df = pd.read_pickle(lipms_covs_file)
print(COV_df.keys())


## Load Deg genes
deg_data = [x.strip() for x in open(f'data/deg_annotation_p.csv', 'r').readlines() if 'MG1655 II' in x]
essential_genes = []
for entry in deg_data:
    #print(entry, entry.split(';'))
    essential_genes += [entry.split(';')[-2].replace('"','')]
print(f'essential_genes: {essential_genes} {len(essential_genes)}')


#### Create supplemental file
"""
Each sheet of this excel spread sheet will be a different buffer and timepoint combination
where:
    the buffers can include: C, CD, CG
    the timepoints can include: R1min, R5min, R2hr, Rall
for a total of 12 sheets
Each sheet will have one row for each gene observed in that buffer-timepoint combination
the columns will be as follows
1. uniprot ID
2. SPA
3. highest SPA CDF percentile where the SPA(i) >= threshold (i.e. a value of 50 indicates this protein had an SPA in the untreated sample equal to or greater than the 50th percentile of the CDF)
4. LiP-MS coverage
5. high quality crystal structure
6. non-covalent entanglement present in crystal structure
7. high quality AF structure
8. non-covalent entanglement present in high quality AF structure
9. essential or not
"""
master_dict = {}
for buff in ['C', 'CD', 'CG']:
    for timepoint in ['R1min', 'R5min', 'R2hr', 'Rall']:

        df = {'UniprotID':[], 'SPA':[], 'Highest SPA percentile':[], 'LiP-MS coverage':[], 'Essential':[],
              'HQ EXP structure':[], 'NCLE in HQ EXP structure':[], 'Knot in HQ EXP structure':[],
              'HQ AF structure':[], 'NCLE in HQ AF structure':[], 'Knot in HQ AF structure':[]}
        key = (buff, timepoint)
        print(key)

        key_SPA_df = SPA_df[key]
        #print(key_SPA_df.keys())

        for k,v in COV_df[key].items():
            if k[0] == 0:
                key_COV_df = v
        #key_COV_df = COV_df[key][(0,0)]
        #print(key_COV_df)

        for gene, N_SPA in key_SPA_df[(0,0)].values:
            cov = key_COV_df[key_COV_df['Accession'] == gene]['coverage'].values[0]


            # Find the SPA threshold where this genes SPA is >= threshold
            threshold = np.nan
            for k,v in key_SPA_df.items():
                if gene in v['Accession'].values:
                    threshold = k[0]
                if gene not in v['Accession'].values:
                    break
            

            # Determine if there is an EXP structure present and whether it has an NCLE and a knot
            EXP_feature_file = glob.glob(f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gen_proteome_features_EXP/res_features_lib/{gene}_*')
            #print(EXP_feature_file)

            if len(EXP_feature_file) != 0:
                EXP_df = pd.read_csv(EXP_feature_file[0], sep="|")
                EXP = True
                EXP_NCLE = EXP_df['ent_present'].any()
                pdb = EXP_feature_file[0].split('/')[-1].split('_')[1]
                chain = EXP_feature_file[0].split('/')[-1].split('_')[2]
                EXP_knot = f'{pdb.upper()}{chain}' in EXP_knots
            else:
                EXP = False
                EXP_NCLE = np.nan
                EXP_knot = np.nan



            # Determine if there is an EXP structure present and whether it has an NCLE
            AF_feature_file = glob.glob(f'/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Make_Protein_Feature_Files/Gen_proteome_features_AF/res_features_lib/{gene}_*')

            if len(AF_feature_file) != 0:
                AF_df = pd.read_csv(AF_feature_file[0], sep="|")
                AF = True
                AF_NCLE = AF_df['ent_present'].any()
                AF_knot = gene in AF_knots
            else:
                AF = False
                AF_NCLE = np.nan
                AF_knot = np.nan


            # Determine if it is essential
            if gene in essential_genes:
                ESS = True
            else:
                ESS = False
            #print(gene, N_SPA, threshold, cov, EXP, EXP_NCLE, AF, AF_NCLE, ESS)

            # Determine if the gene is a knot


            df['UniprotID'] += [gene]
            df['SPA'] += [N_SPA]
            df['Highest SPA percentile'] += [threshold]
            df['LiP-MS coverage'] += [cov] 
            df['Essential'] += [ESS] 
            df['HQ EXP structure'] += [EXP]
            df['NCLE in HQ EXP structure'] += [EXP_NCLE]
            df['HQ AF structure'] += [AF]
            df['NCLE in HQ AF structure'] += [AF_NCLE]
            df['Knot in HQ EXP structure'] += [EXP_knot]
            df['Knot in HQ AF structure'] += [AF_knot]

        master_dict[key] = pd.DataFrame(df)

### SAVE the excel file
with pd.ExcelWriter(f'{outpath}/Master_combined_gene_lists.xlsx', engine='xlsxwriter') as writer:
    for key, df in master_dict.items():
        df.to_excel(writer, sheet_name='_'.join(key), index=False)

print(f'SAVED: {outpath}/Master_combined_gene_lists.xlsx')
print('NORMAL TERMINATION')



