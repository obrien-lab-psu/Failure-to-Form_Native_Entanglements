import sys,os
import glob
import numpy as np
import pandas as pd
import argparse
pd.set_option('display.max_rows', 1000)  # Display 10 rows

def load_express_data(express_file, buff, timepoint, threshold):
    """
    This function takes the SPA pickle file made in the FLiPPR analysis and gets the set of genes that corresponds to 
    (1) the desired buffer: C, CD, CG
    (2) the desired timepoint: R1min, R5min, R2hr, Rall(the spa for this is the per gene average across all timepoints)
    (3) the desired CDF threshold: 0, 10, ..., 90
    """

    print(f'load_express_data...')
    print(express_file, buff, timepoint)
    df = pd.read_pickle(express_file)

    buff_timepoint_df = df[(buff, timepoint)]
    #buff_timepoint_df = df[(buff, 'Rall')]
    for percentile_key, percentile_df in buff_timepoint_df.items():
        print(percentile_key, len(percentile_df))
        if percentile_key[0] == threshold:

            express_genes = set(percentile_df['Accession'])
    print(f'Number of genes in {buff} {timepoint} {threshold} = {len(express_genes)}')
    return express_genes


def load_LiPMS_cov(filepath, timepoint):
    """
    This function takes the LiPMS coverage pickle file made in the FLiPPR analysis and gets the set of genes that corresponds to 
    (1) the desired buffer: C, CD, CG
    (2) the desired timepoint: R1min, R5min, R2hr, Rall(the spa for this is the per gene average across all timepoints)
    it then gets all proteins and makes a dataframe of their coverage as we do not want to select based on the CDF but on the raw value of the coverage
    """
    df = pd.read_pickle(filepath)
    coverage = {}
    for b in ['C', 'CD', 'CG']:
        buff_df = df[(b, timepoint)]
        buff_cov = pd.DataFrame()
        for percentile_key, percentile_df in buff_df.items():
            #print(buff, percentile_key, len(percentile_df))
            buff_cov = pd.concat((buff_cov, percentile_df))

        buff_cov = buff_cov.drop_duplicates()
        buff_cov = buff_cov.reset_index()
        coverage[b] = buff_cov

    for b,v in coverage.items():
        print(b,v)
    return coverage


###Load user defined arguments
parser = argparse.ArgumentParser(description="Process user specified arguments")
parser.add_argument("-f", "--feature_files", type=str, required=True, help="Path to residue feature files")
parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory")
parser.add_argument("-t", "--tag", type=str, required=True, help="Tag for output filenames")
parser.add_argument("-e", "--express_file", type=str, required=True, help='path to expression control file or Total')
parser.add_argument("-l", "--lipms_covs_file", type=str, required=True, help='path to lipms coverage file')
parser.add_argument("-b", "--buffer", type=str, required=True, help="Buffer system to use: C, CD, CG")
parser.add_argument("-th", "--threshold", type=int, required=True, help="threshold for either SPA: 0 - 90 in 10 int steps")
parser.add_argument("-tp", "--timepoint", type=str, required=True, help="R1min, R5min, R2hr, Rall")
parser.add_argument("-m", "--mask", type=str, required=False, help='Global mask for genes. if present no gene will pass unless in mask regardless of other filters')
parser.add_argument("-c", "--lipms_cov_threshold", type=int, required=True, help="threshold for either LiPMS-COV: 0 - 90 in 10 int steps")
parser.add_argument("-k", "--knots", type=str, required=True, help="List of genes to ignore that contain knots")
args = parser.parse_args()

feature_files = args.feature_files
outpath = args.outpath
tag = args.tag
express_file = args.express_file
buff = args.buffer
timepoint = args.timepoint
threshold = args.threshold
mask = args.mask
lipms_covs_file = args.lipms_covs_file
lipms_cov_threshold = args.lipms_cov_threshold
knots = pd.read_csv(args.knots, sep=';')
print('Knots:\n{knots}')

if timepoint == 'Rall':
    timepoints = ['R1min', 'R5min', 'R2hr']
else:
    timepoints = [timepoint]

if mask != None:
    print(f'User specified a mask that will be used: {mask}')
    mask = np.loadtxt(mask, dtype=str)
else:
    print(f'User DID NOT specified a mask that will be use. All feature file genes will be used')
    mask = [f.split('/')[-1].split('_')[0] for f in glob.glob(os.path.join(feature_files, '*'))]
print(f'Mask: {mask} {len(mask)}')

## make outdirectory if it doesnt exists yet
if not os.path.exists(f'{outpath}'):
    os.makedirs(f'{outpath}')
    print(f'Made output directories {outpath}')


## if buff != None then load expression data. else if it does == None then set expression data to a list of all 1294 genes.
###Load PSM data and get genes that only meet the PSM percentile specified by the user
if buff in ['C', 'CD', 'CG']:
    express_data = load_express_data(express_file, buff, timepoint, threshold)
    #print(f'express_data: {express_data}, {len(express_data)}')

elif buff in ['Total']:
    threshold = 'Total'
    lipms_cov_threshold = 'Total'
else:
    raise ValueError("Buffer argument is not valid: C, CD, CG, Total are the only valid options")
    quit()

## Load LiPMS coverage data
LiPMS_cov = load_LiPMS_cov(lipms_covs_file, timepoint)


### Gene list definitions (only those in the expression data at the given buffer and threshold threshold)
all_genes = [] # all genes no conditions

ent_genes = [] # gene has ent
nonent_genes = [] # gene does not have ent

refolded_genes = [] # gene has LESS than {num_cuts} significatnt cutsites
nonrefolded_genes = [] # gene has atleast {num_cuts} significatnt cutsites

essential_genes = [] # gene is essential
nonessential_genes = [] # gene is NOT essential

essential_ent_genes = [] # gene has ent AND is essential
nonessential_ent_genes = [] # gene has ent AND is NOT essential

essential_nonent_genes = [] # gene has NO ent AND is essential
nonessential_nonent_genes = [] #gene has NO ent AND is not essential

Total_ent_genes = 0
for f_i, feature_file in enumerate(glob.glob(os.path.join(feature_files, '*'))):
    print(f'{"-"*30}\n{f_i} {feature_file}\n{"-"*30}')

    # load the feature datafile for this gene
    feature_data = pd.read_csv(feature_file, sep='|')
    #print(feature_data[['gene', 'ent_present', 'essential', 'cut_str']])

    gene = feature_file.split('/')[-1].split('_')[0]
    pdb = feature_file.split('/')[-1].split('_')[1]
    chain = feature_file.split('/')[-1].split('_')[2]
    print(f'gene: {gene} pdb: {pdb} chain: {chain}')
    if pdb != 'AF':
        knots_df = knots[(knots['pdbid'] == pdb.lower()) & (knots['chain'] == chain)]
        if len(knots_df) != 0:
            print(f'Knot found for {gene} {pdb} {chain} and will be ignored\n{knots_df}')
            continue
    else:
        knots_df = knots[knots['gene'] == gene]
        if len(knots_df) != 0:
            print(f'Knot found for {gene} {pdb} {chain} and will be ignored\n{knots_df}')
            continue       

    ## Check if gene is in mask
    if gene in mask:
        print(f'{gene} is in mask')
    else:
        print(f'{gene} not in mask and will be skipped')
        continue


    ############################################################
    # Get gene essentiality
    essential = feature_data['essential'].sum()
    if essential > 0:
        essential = True
    elif essential == 0:
        essential = False
    print(f'essential: {essential}')


    ############################################################
    # Get whether gene has entanglement or not
    ent_present = feature_data['ent_present'].sum()
    if ent_present > 0:
        ent_present = True
        Total_ent_genes += 1
    elif ent_present == 0:
        ent_present = False
    print(f'ent_present: {ent_present}')


    ############################################################
    # check if the protein was refolded or nonrefolded by looking for a number of cut sites in the current buffer greater than or equal to the num_cuts
    cut_strings = [c for c in feature_data['cut_str'].values if isinstance(c, str)]
    print(f'Gene: {gene} cut string is {cut_strings}')
    cuts = []
    for res_cut_strs in cut_strings:
        cut_strs = res_cut_strs.split('/')
        #print(cut_strs)

        for cut_str in cut_strs:
            cut_str = cut_str.split(',')
            #print(cut_str)
            if cut_str[0] == buff and cut_str[1] in timepoints:
                if cut_str[2] not in cuts:
                    cuts += [cut_str[2]]
    ncuts = len(cuts)

    if ncuts >= 1:
        nonrefoldable = True
    else:
        nonrefoldable = False
    print(f'Total unique cuts in buff: {buff} for gene: {gene} is {ncuts} {cuts}\nnonrefoldable = {nonrefoldable}')

    ############################################################
    ## if gene is in expression data record it into buffer, timepoint, and SPA threshold specific lists
    ## if buff == Total disregard expression and count all genes
    if buff != 'Total':
        if gene in express_data:
            print(f'Gene {gene} is in express cut off list with percentile cutoff: {threshold}')
            express_gene = True
        else:
            express_gene = False
    else:
        print(f'Buffer is specified as Total so all genes are valid no matter SPA')
        express_gene = True


    if express_gene:

        ############################################################
        ## Get LiPMS coverage for gene
        if buff != 'Total':
            cov = LiPMS_cov[buff]
            cov = cov[cov['Accession'] == gene]['coverage'].values[0]
            print(f'LiPMS coverage of gene {gene} in buffer {buff} is {cov}')
            if cov >= lipms_cov_threshold:
                print(f'Gene {gene} has coverage in LiPMS experiment greater than or equal to the specified threshold {lipms_cov_threshold}')
                coverage = True
            else:
                coverage = False
        else:
            print(f'Buffer is specified as Total so all genes are valid no matter LiPMS coverage')
            coverage = True

        if coverage:
            all_genes += [gene]

            ############################################################
            ## check if gene is nonrefolable
            if nonrefoldable == True:
                nonrefolded_genes += [gene]
            else:
                refolded_genes += [gene]


            ############################################################
            #check if gene has entanglement or not
            if ent_present:
                print(f'Gene {gene} has an entanglement')
                ent_genes += [gene]

                #check if it is a significant gene
                if essential:
                    essential_ent_genes += [gene]
                else:
                    nonessential_ent_genes += [gene]

            else:
                print(f'Gene {gene} does NOT have an entanglement')
                nonent_genes += [gene]

                #check if it is a significant gene
                if essential:
                    print(f'Gene {gene} is considered essential')
                else:
                    print(f'Gene {gene} is NOT considered essential')


            ############################################################
            #check if it is a significant gene
            if essential == True:
                print(f'Gene {gene} is considered essential')
                essential_genes += [gene]
            else:
                print(f'Gene {gene} is NOT considered essential')
                nonessential_genes += [gene]


#save lists for buffer, timepoint, expression percentile specific lists
np.savetxt(f'{outpath}{tag}_{buff}_{timepoint}_spa{threshold}_LiPMScov{lipms_cov_threshold}_all_genes.txt', all_genes, fmt='%s')
print(f'Saved: {outpath}{tag}_{buff}_{timepoint}_spa{threshold}_LiPMScov{lipms_cov_threshold}_all_genes.txt {len(all_genes)}')

np.savetxt(f'{outpath}{tag}_{buff}_{timepoint}_spa{threshold}_LiPMScov{lipms_cov_threshold}_refolded_genes.txt', refolded_genes, fmt='%s')
print(f'Saved: {outpath}{tag}_{buff}_{timepoint}_spa{threshold}_LiPMScov{lipms_cov_threshold}_refolded_genes.txt {len(refolded_genes)}')

np.savetxt(f'{outpath}{tag}_{buff}_{timepoint}_spa{threshold}_LiPMScov{lipms_cov_threshold}_nonrefolded_genes.txt', nonrefolded_genes, fmt='%s')
print(f'Saved: {outpath}{tag}_{buff}_{timepoint}_spa{threshold}_LiPMScov{lipms_cov_threshold}_nonrefolded_genes.txt {len(nonrefolded_genes)}')

np.savetxt(f'{outpath}{tag}_{buff}_{timepoint}_spa{threshold}_LiPMScov{lipms_cov_threshold}_essential_genes.txt', essential_genes, fmt='%s')
print(f'Saved: {outpath}{tag}_{buff}_{timepoint}_spa{threshold}_LiPMScov{lipms_cov_threshold}_essential_genes.txt {len(essential_genes)}')
np.savetxt(f'{outpath}{tag}_{buff}_{timepoint}_spa{threshold}_LiPMScov{lipms_cov_threshold}_nonessential_genes.txt', nonessential_genes, fmt='%s')
print(f'Saved: {outpath}{tag}_{buff}_{timepoint}_spa{threshold}_LiPMScov{lipms_cov_threshold}_nonessential_genes.txt {len(nonessential_genes)}')

np.savetxt(f'{outpath}{tag}_{buff}_{timepoint}_spa{threshold}_LiPMScov{lipms_cov_threshold}_essential_ent_genes.txt', essential_ent_genes, fmt='%s')
print(f'Saved: {outpath}{tag}_{buff}_{timepoint}_spa{threshold}_LiPMScov{lipms_cov_threshold}_essential_ent_genes.txt {len(essential_ent_genes)}')
np.savetxt(f'{outpath}{tag}_{buff}_{timepoint}_spa{threshold}_LiPMScov{lipms_cov_threshold}_nonessential_ent_genes.txt', nonessential_ent_genes, fmt='%s')
print(f'Saved: {outpath}{tag}_{buff}_{timepoint}_spa{threshold}_LiPMScov{lipms_cov_threshold}_nonessential_ent_genes.txt {len(nonessential_ent_genes)}')

essential_nonent_genes = [g for g in essential_genes if g not in ent_genes]
nonessential_nonent_genes = [g for g in nonessential_genes if g not in ent_genes]
np.savetxt(f'{outpath}{tag}_{buff}_{timepoint}_spa{threshold}_LiPMScov{lipms_cov_threshold}_essential_nonent_genes.txt', essential_nonent_genes, fmt='%s')
print(f'Saved: {outpath}{tag}_{buff}_{timepoint}_spa{threshold}_LiPMScov{lipms_cov_threshold}_essential_nonent_genes.txt {len(essential_nonent_genes)}')
np.savetxt(f'{outpath}{tag}_{buff}_{timepoint}_spa{threshold}_LiPMScov{lipms_cov_threshold}_nonessential_nonent_genes.txt', nonessential_nonent_genes, fmt='%s')
print(f'Saved: {outpath}{tag}_{buff}_{timepoint}_spa{threshold}_LiPMScov{lipms_cov_threshold}_nonessential_nonent_genes.txt {len(nonessential_nonent_genes)}')

np.savetxt(f'{outpath}{tag}_{buff}_{timepoint}_spa{threshold}_LiPMScov{lipms_cov_threshold}_ent_genes.txt', ent_genes, fmt='%s')
print(f'Saved: {outpath}{tag}_{buff}_{timepoint}_spa{threshold}_LiPMScov{lipms_cov_threshold}_ent_genes.txt {len(ent_genes)}')
np.savetxt(f'{outpath}{tag}_{buff}_{timepoint}_spa{threshold}_LiPMScov{lipms_cov_threshold}_nonent_genes.txt', nonent_genes, fmt='%s')
print(f'Saved: {outpath}{tag}_{buff}_{timepoint}_spa{threshold}_LiPMScov{lipms_cov_threshold}_nonent_genes.txt {len(nonent_genes)}')


outdf = {'buff':[buff],
        'timepoint':[timepoint],
        'spa_threshold':[threshold],
        'LiPMScov_threshold':[lipms_cov_threshold],
        'all_genes_n':[len(all_genes)],
        'essential_genes_n':[len(essential_genes)],
        'nonessential_genes_n':[len(nonessential_genes)],
        'essential_ent_genes_n':[len(essential_ent_genes)],
        'nonessential_ent_genes_n':[len(nonessential_ent_genes)],
        'essential_nonent_genes_n':[len(essential_nonent_genes)],
        'nonessential_nonent_genes_n':[len(nonessential_nonent_genes)],
        'refolded_genes_n':[len(refolded_genes)],
        'nonrefolded_genes_n':[len(nonrefolded_genes)],
        'ent_genes_n':[len(ent_genes)],
        'nonent_genes_n':[len(nonent_genes)]}

for k,v in outdf.items():
    print(k,v)
outdf = pd.DataFrame(outdf)
outdf_file = f'{outpath}{tag}_{buff}_{timepoint}_spa{threshold}_LiPMScov{lipms_cov_threshold}_stats.csv'
outdf.to_csv(outdf_file, sep='|', index=False)
print(f'SAVED: {outdf_file}')

print(f'Quality check: Total ENT genes: {Total_ent_genes} should be 921')
