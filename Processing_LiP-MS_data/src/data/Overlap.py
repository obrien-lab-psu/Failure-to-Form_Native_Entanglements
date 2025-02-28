import pandas as pd
from matplotlib_venn import venn3, venn2
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, combine_pvalues
import pickle
import numpy as np
import sys,os
import glob
import re
import argparse
import warnings
from tqdm import tqdm


class DataAnalysis:
    """
    A class to handle the data analysis process including encoding, regression, and statistical tests.
    """

    #######################################################################################
    def __init__(self, inpfiles, outpath, tag, gene_lists):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
        - resFeat_files (str): Path to residue feature files.
        - outpath (str): Path to the output directory.
        """
        if gene_lists != 'None':
            gene_list_files = glob.glob(gene_lists)
        else:
            gene_list_files = None

        files = {}
        for f in glob.glob(inpfiles):
            buff = f.split('/')[-1].split('_')[0]
            timepoint = f.split('/')[-1].split('_')[1]

            total_df = pd.read_csv(f)
            df = total_df[['Accession', 'proteinaseKsite']]

            ## if a gene list was specified mask the data
            if gene_list_files != None:
                gene_list = [f for f in gene_list_files if f'{buff}_Rall' in f]
                print(buff, timepoint, gene_list)
                if len(gene_list) == 1:
                    genes = np.loadtxt(gene_list[0], dtype=str)
                else:
                    raise ValueError(f'Multiple gene lists found: {gene_list}')

                df = df[df['Accession'].isin(genes)]

            # Extract the numeric part and create a new column
            df['site_number'] = df['proteinaseKsite'].str.extract('(\d+)')

            # Convert the extracted numbers to integers (optional, if needed)
            df['site_number'] = df['site_number'].astype(int)

            # Check for consecutive numbers
            df['consecutive'] = df['site_number'].diff() == 1

            print(f, buff, timepoint)
            print(df)
            
            if buff not in files:
                files[buff] = {timepoint:df}
            else:
                files[buff][timepoint] = df

    
        self.inpfiles = files
        self.outpath = outpath
        self.tag = tag

        if not os.path.exists(f'{self.outpath}'):
            os.makedirs(f'{self.outpath}')
            print(f'Made output directories {self.outpath}')

        self.data_path = os.path.join(f'{self.outpath}', 'Data/')
        if not os.path.exists(f'{self.data_path}'):
            os.makedirs(f'{self.data_path}')
            print(f'Made output directories {self.data_path}')

        self.plots_path = os.path.join(f'{self.outpath}', 'Plots/')
        if not os.path.exists(f'{self.plots_path}'):
            os.makedirs(f'{self.plots_path}')
            print(f'Made output directories {self.plots_path}')

        self.reps = 10000

    #######################################################################################
    def setup_logging(self, logfile):
        """
        Sets up the logging configuration.

        Returns:
        - logger (logging.Logger): Configured logger.
        """
        logging.basicConfig(filename=logfile, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        logger = logging.getLogger(__name__)
        return logger

    #######################################################################################
    def gene_overlap_3timepoint(self,):    
        for buff in ['C', 'CD', 'CG']:
            print(f'Calcutating gene overlap between {buff} timepoints')
            R1min = set(self.inpfiles[buff]['R1min']['Accession'].values)
            R1min_size = len(R1min)
            #print(R1min, R1min_size)
            R5min = set(self.inpfiles[buff]['R5min']['Accession'].values)
            R5min_size = len(R5min)
            #print(R5min, R5min_size)
            R2hr = set(self.inpfiles[buff]['R2hr']['Accession'].values)
            R2hr_size = len(R2hr)
            #print(R2hr, R2hr_size)

            # Calculate the total number of unique elements
            total_elements = len(R1min.union(R5min).union(R2hr))
            print(f'total_elements: {total_elements}')
            total_genes = R1min.union(R5min).union(R2hr)

            ## Get the confidence intervals on the overlap percentages and prob of occuring by random chance
            gene_overlap = {'R1min':{'boot':[], 'perm':[], 'GT':len(R1min.difference(R5min).difference(R2hr))/total_elements, 'sets':[R1min]}, 
                            'R5min':{'boot':[], 'perm':[], 'GT':len(R5min.difference(R1min).difference(R2hr))/total_elements, 'sets':[R5min]},
                            'R2hr':{'boot':[], 'perm':[], 'GT':len(R2hr.difference(R5min).difference(R1min))/total_elements, 'sets':[R2hr]},
                            'R1min-R5min':{'boot':[], 'perm':[], 'GT':len(R1min.intersection(R5min).difference(R2hr))/total_elements, 'sets':[R1min, R5min]}, 
                            'R1min-R2hr':{'boot':[], 'perm':[], 'GT':len(R1min.intersection(R2hr).difference(R5min))/total_elements, 'sets':[R1min, R2hr]},
                            'R5min-R2hr':{'boot':[], 'perm':[], 'GT':len(R5min.intersection(R2hr).difference(R1min))/total_elements, 'sets':[R5min, R2hr]}, 
                            'R1min-R5min-R2hr':{'boot':[], 'perm':[], 'GT':len(R1min.intersection(R5min).intersection(R2hr))/total_elements, 'sets':[R1min, R5min, R2hr]}}


            #############################
            ## do permutation testing
            print('PERMUTATION TESTING: Gene overlap')
            total_dataset = list(R1min) + list(R5min) + list(R2hr)
            total_dataset_size = len(total_dataset)
            for p in tqdm(range(self.reps)):
                ptotal_dataset = np.random.permutation(total_dataset)
                pR1min = ptotal_dataset[:R1min_size]
                pR5min = ptotal_dataset[R1min_size:R1min_size + R5min_size]
                pR2hr = ptotal_dataset[R1min_size + R5min_size:]
                #print(len(pR1min), len(pR5min), len(pR2hr))

                #pR1min_R5min = len(set(pR1min).intersection(pR5min))/total_elements
                #pR1min_R2hr = len(set(pR1min).intersection(pR2hr))/total_elements
                #pR5min_R2hr = len(set(pR5min).intersection(pR2hr))/total_elements
                #pR1min_R5min_R2hr = len(set(pR1min).intersection(pR5min).intersection(pR2hr))/total_elements

                gene_overlap['R1min']['perm'] += [len(set(pR1min).difference(pR5min).difference(pR2hr))/total_elements]
                gene_overlap['R5min']['perm'] += [len(set(pR5min).difference(pR1min).difference(pR2hr))/total_elements]
                gene_overlap['R2hr']['perm'] += [len(set(pR2hr).difference(pR5min).difference(pR1min))/total_elements]
                gene_overlap['R1min-R5min']['perm'] += [len(set(pR1min).intersection(pR5min).difference(pR2hr))/total_elements]
                gene_overlap['R1min-R2hr']['perm'] += [len(set(pR1min).intersection(pR2hr).difference(pR5min))/total_elements]
                gene_overlap['R5min-R2hr']['perm'] += [len(set(pR5min).intersection(pR2hr).difference(pR1min))/total_elements]
                gene_overlap['R1min-R5min-R2hr']['perm'] += [len(set(pR1min).intersection(pR5min).intersection(pR2hr))/total_elements]

            for setstr, set_data in gene_overlap.items():
                pvalue = 1
                for perm in set_data['perm']:
                    if perm >= set_data['GT']:
                        pvalue += 1

                pvalue /= (self.reps + 1)
                print(buff, setstr, set_data['GT'], pvalue)
                gene_overlap[setstr]['pvalue'] = pvalue

            #############################
            ## do bootstrapping
            # Inialize a dataframe that contains unique genes and their classifications
            df = {'gene':[], 'sets_str':[]}
            for genei, gene in enumerate(total_genes):
                sets = []
                if gene in R1min:
                    sets += ['R1min']
                if gene in R5min:
                    sets += ['R5min']
                if gene in R2hr:
                    sets  += ['R2hr']
            
                sets_str = ''
                if len(sets) == 3:
                    sets_str = 'R1min-R5min-R2hr'
                
                if len(sets) == 1:
                    sets_str = sets[0]

                if len(sets) == 2:
                    sets_str = '-'.join(sets)
                #print(genei, gene, sets, sets_str)
                
                df['gene'] += [gene]
                df['sets_str'] += [sets_str]
            df = pd.DataFrame(df)

            counts = df['sets_str'].value_counts()
            counts /= total_elements
            
            # bootstrap the fraction of these classes
            print('BOOTSTRAPPIG BOUNDS: Gene overlap')
            for b in tqdm(range(self.reps)):
                bdf = df.sample(frac=1, replace=True)
                bcounts = bdf['sets_str'].value_counts()
                bcounts /= total_elements
                #print(bcounts)
                for k,v in bcounts.items():
                    gene_overlap[k]['boot'] += [v]

            outdf = {'buff':[], 'OverlapClass':[], 'Overlap':[], 'Bounds':[], 'pvalue':[], 'count':[]}
            for setstr, set_data in gene_overlap.items():
                if len(set_data['boot']) != 0:
                    lb = np.percentile(set_data['boot'], 2.5)
                    ub = np.percentile(set_data['boot'], 97.5)
                else:
                    lb = 0.0
                    ub = 0.0
                gene_overlap[setstr]['bounds'] = (lb, ub)
                GT = set_data['GT']
                pvalue = set_data['pvalue']

                outdf['buff'] += [buff]
                outdf['OverlapClass'] += [setstr]
                outdf['Overlap'] += [f'{GT:.3f}']
                outdf['Bounds'] += [f'({lb:.3f}, {ub:.3f})']
                outdf['pvalue'] += [f'{pvalue:.3e}']
                outdf['count'] += [int(GT*total_elements)]

            for setstr, set_data in gene_overlap.items():
                print(buff, setstr, set_data['GT'], set_data['bounds'], set_data['pvalue'])

            outdf = pd.DataFrame(outdf)
            print(outdf)

            # Save the dataframe to the data path
            outfile = os.path.join(self.data_path, f'GeneOverlapTableVenn3Diagram_{self.tag}_{buff}.csv')
            outdf.to_csv(outfile, index=False)
            print(f'SAVED: {outfile}')

            # Create a larger figure
            plt.figure(figsize=(12, 8))  # Increase the size by setting width=10 inches and height=8 inches

            # Calculate the Venn diagram
            venn = venn3([R1min, R5min, R2hr], ('R1min', 'R5min', 'R2hr'))
            
            # Update the labels with the percentages
            for idx, label in enumerate(venn.subset_labels):
                print(idx, label)
                if label:  # Check if the label exists
                    size = int(label.get_text())
                    percentage = size / total_elements * 100
                    label.set_text(f'{size} ({percentage:.3f}%)')
            

            # Add the table to the right of the Venn diagram
            plt.subplots_adjust(left=0.3, right=0.8)
            ax = plt.gca()
            ax_table = plt.table(cellText=outdf.values, colLabels=outdf.columns, colWidths=[0.15, 0.15, 0.15, 0.15, 0.15, 0.15], cellLoc='center', loc='right')
            
            # Increase font size
            ax_table.auto_set_font_size(False)  # Disable automatic font size scaling
            ax_table.set_fontsize(8)  # Set the desired font size

            # Display the plot
            plt.suptitle(f'{self.tag} | {buff} | GeneOverlap')
            plt.tight_layout()
            #plt.show()
            outfile = os.path.join(self.plots_path, f'GeneOverlapVenn3Diagram_{self.tag}_{buff}.pdf')
            plt.savefig(outfile)
            print(f'SAVED: {outfile}')
            plt.close()

    #######################################################################################
    def gene_overlap_2timepoint(self,):    
        for buff in ['C', 'CD', 'CG']:
            print(f'Calcutating gene overlap between {buff} timepoints')

            R5min = set(self.inpfiles[buff]['R5min']['Accession'].values)
            R5min_size = len(R5min)
            #print(R5min, R5min_size)
            R2hr = set(self.inpfiles[buff]['R2hr']['Accession'].values)
            R2hr_size = len(R2hr)
            #print(R2hr, R2hr_size)

            # Calculate the total number of unique elements
            total_elements = len(R5min.union(R2hr))
            print(f'total_elements: {total_elements}')
            total_genes = R5min.union(R2hr)

            ## Get the confidence intervals on the overlap percentages and prob of occuring by random chance
            gene_overlap = {'R5min':{'boot':[], 'perm':[], 'GT':len(R5min.difference(R2hr))/total_elements, 'sets':[R5min]},
                            'R2hr':{'boot':[], 'perm':[], 'GT':len(R2hr.difference(R5min))/total_elements, 'sets':[R2hr]},
                            'R5min-R2hr':{'boot':[], 'perm':[], 'GT':len(R5min.intersection(R2hr))/total_elements, 'sets':[R5min, R2hr]}}

            #############################
            ## do permutation testing
            print('PERMUTATION TESTING: Gene overlap')
            total_dataset = list(R5min) + list(R2hr)
            total_dataset_size = len(total_dataset)
            for p in tqdm(range(self.reps)):
                ptotal_dataset = np.random.permutation(total_dataset)
                pR5min = ptotal_dataset[:R5min_size]
                pR2hr = ptotal_dataset[R5min_size:]

                gene_overlap['R5min']['perm'] += [len(set(pR5min).difference(pR2hr))/total_elements]
                gene_overlap['R2hr']['perm'] += [len(set(pR2hr).difference(pR5min))/total_elements]
                gene_overlap['R5min-R2hr']['perm'] += [len(set(pR5min).intersection(pR2hr))/total_elements]

            for setstr, set_data in gene_overlap.items():
                pvalue = 1
                for perm in set_data['perm']:
                    if perm >= set_data['GT']:
                        pvalue += 1

                pvalue /= (self.reps + 1)
                print(buff, setstr, set_data['GT'], pvalue)
                gene_overlap[setstr]['pvalue'] = pvalue

            #############################
            ## do bootstrapping
            # Inialize a dataframe that contains unique genes and their classifications
            df = {'gene':[], 'sets_str':[]}
            for genei, gene in enumerate(total_genes):
                sets = []

                if gene in R5min:
                    sets += ['R5min']
                if gene in R2hr:
                    sets  += ['R2hr']
            
                sets_str = ''
                
                if len(sets) == 1:
                    sets_str = sets[0]

                if len(sets) == 2:
                    sets_str = '-'.join(sets)
                #print(genei, gene, sets, sets_str)
                
                df['gene'] += [gene]
                df['sets_str'] += [sets_str]
            df = pd.DataFrame(df)

            counts = df['sets_str'].value_counts()
            counts /= total_elements
            
            # bootstrap the fraction of these classes
            print('BOOTSTRAPPIG BOUNDS: Gene overlap')
            for b in tqdm(range(self.reps)):
                bdf = df.sample(frac=1, replace=True)
                bcounts = bdf['sets_str'].value_counts()
                bcounts /= total_elements
                #print(bcounts)
                for k,v in bcounts.items():
                    gene_overlap[k]['boot'] += [v]

            outdf = {'buff':[], 'OverlapClass':[], 'Overlap':[], 'Bounds':[], 'pvalue':[], 'count':[]}
            for setstr, set_data in gene_overlap.items():
                if len(set_data['boot']) != 0:
                    lb = np.percentile(set_data['boot'], 2.5)
                    ub = np.percentile(set_data['boot'], 97.5)
                else:
                    lb = 0.0
                    ub = 0.0
                gene_overlap[setstr]['bounds'] = (lb, ub)
                GT = set_data['GT']
                pvalue = set_data['pvalue']

                outdf['buff'] += [buff]
                outdf['OverlapClass'] += [setstr]
                outdf['Overlap'] += [f'{GT:.3f}']
                outdf['Bounds'] += [f'({lb:.3f}, {ub:.3f})']
                outdf['pvalue'] += [f'{pvalue:.3e}']
                outdf['count'] += [int(GT*total_elements)]

            for setstr, set_data in gene_overlap.items():
                print(buff, setstr, set_data['GT'], set_data['bounds'], set_data['pvalue'])

            outdf = pd.DataFrame(outdf)
            print(outdf)

            # Save the dataframe to the data path
            outfile = os.path.join(self.data_path, f'GeneOverlapTableVenn2Diagram_{self.tag}_{buff}.csv')
            outdf.to_csv(outfile, index=False)
            print(f'SAVED: {outfile}')

            # Create a larger figure
            plt.figure(figsize=(12, 8))  # Increase the size by setting width=10 inches and height=8 inches

            # Calculate the Venn diagram
            venn = venn2([R5min, R2hr], ('R5min', 'R2hr'))
            
            # Update the labels with the percentages
            for idx, label in enumerate(venn.subset_labels):
                print(idx, label)
                if label:  # Check if the label exists
                    size = int(label.get_text())
                    percentage = size / total_elements * 100
                    label.set_text(f'{size} ({percentage:.3f}%)')
            

            # Add the table to the right of the Venn diagram
            plt.subplots_adjust(left=0.3, right=0.8)
            ax = plt.gca()
            ax_table = plt.table(cellText=outdf.values, colLabels=outdf.columns, colWidths=[0.15, 0.15, 0.15, 0.15, 0.15, 0.15], cellLoc='center', loc='right')
            
            # Increase font size
            ax_table.auto_set_font_size(False)  # Disable automatic font size scaling
            ax_table.set_fontsize(8)  # Set the desired font size

            # Display the plot
            plt.suptitle(f'{self.tag} | {buff} | GeneOverlap')
            plt.tight_layout()
            #plt.show()
            outfile = os.path.join(self.plots_path, f'GeneOverlapVenn2Diagram_{self.tag}_{buff}.pdf')
            plt.savefig(outfile)
            print(f'SAVED: {outfile}')
            plt.close()

    #######################################################################################
    def PKsite_overlap_3timepoint(self,): 
        offset = 5   
        for buff in ['C', 'CD', 'CG']:
            print(f'Calcutating PKsite overlap between {buff} timepoints')

            # Get R1min cut sites with +/- 3 residue buffer
            def get_extended_cutsites(df):
                out = {}
                for rowi, row in df.iterrows():
                    gene = row['Accession']
                    site = row['site_number']
                    if gene not in out:
                        out[gene] = {}
                    out[gene][site] = []

                    for s in range(row['site_number'] - offset, row['site_number'] + offset + 1):
                        out[gene][site] += [s]
                return out

            R1min = get_extended_cutsites(self.inpfiles[buff]['R1min'])
            R1min_size = len(R1min)
            print(f'R1min_size: {R1min_size}')           

            R5min = get_extended_cutsites(self.inpfiles[buff]['R5min'])
            R5min_size = len(R5min)
            print(f'R5min_size: {R5min_size}')
            
            R2hr = get_extended_cutsites(self.inpfiles[buff]['R2hr'])
            R2hr_size = len(R2hr)
            print(f'R2hr_size: {R2hr_size}')

            # Calculate the total number of unique elements
            total_elements = R1min_size + R5min_size + R2hr_size
            print(f'total_elements: {total_elements}')
            
            ###########################################################################
            ## Get the confidence intervals on the overlap percentages and prob of occuring by random chance
            # Define a general function to do the overlap calcs so i dont have to repeate logic
            def get_site_overlap(R1min, R5min, R2hr):
                
                # Define a general function to check overlap of two timepoints
                def check_overlap(t1, t1str, t2, t2str):
                    overlap_strs = []
                    for gene, t1_gene_sites in t1.items():
                        #print('\n', gene, t1_gene_sites)
                        gene_count = 0
                        gene_overlap_strs = []

                        if gene in t2:
                            #print(f'{gene} in second timepoint')
                            t2_gene_sites = t2[gene]
                            #print(t2_gene_sites)

                            # check for atleast 1 overlap 
                            for t1_PKsite, t1_PKsite_ext in t1_gene_sites.items():
                                for t2_PKsite, t2_PKsite_ext in t2_gene_sites.items():
                                    if len(set(t1_PKsite_ext).intersection(t2_PKsite_ext)) != 0:
                                        #print(f'Found overlap between {t1_PKsite} and {t2_PKsite}')
                                        overlap_str = f'{gene}_{t1str}-{t1_PKsite}_{t2str}-{t2_PKsite}'
                                        gene_count += 1
                                        gene_overlap_strs += [overlap_str]
                                        
                            #print(f'{gene} Overlap {t1str}-{t2str} count: {gene_count}\n{gene_overlap_strs}')
                        overlap_strs += gene_overlap_strs
                    return overlap_strs
                    
                # get those sites in R1min and R5min
                R1min_R5min_overlap_strs = check_overlap(R1min, 'R1min', R5min, 'R5min')
                #print(R1min_R5min_overlap_strs)

                # get those sites in R1min and R2hr
                R1min_R2hr_overlap_strs = check_overlap(R1min, 'R1min', R2hr, 'R2hr')
                #print(R1min_R2hr_overlap_strs) 

                # get those sites in R5min and R2hr
                R5min_R2hr_overlap_strs = check_overlap(R5min, 'R5min', R2hr, 'R2hr')
                #print(R5min_R2hr_overlap_strs) 

                # make a master dictionary of all genes and their cutsites across all three timepoints
                master = {}
                for t in [R1min, R5min, R2hr]:
                    for gene, gene_sites in t.items():
                        pk_sites = gene_sites.keys()
                        if gene not in master:
                            master[gene] = {}
                        for site in pk_sites:
                            master[gene][site] = []
                
                # For each gene and its sites in the master determine their classifications
                for gene, sites in master.items():
                    for site in sites:
                        #print('\n', gene, site)

                        # check if site had overlap between R1min R5min
                        check15 = [s for s in R1min_R5min_overlap_strs if gene in s and f'-{site}' in s]
                        #print(f'check15: {check15}')
                        check12 = [s for s in R1min_R2hr_overlap_strs if gene in s and f'-{site}' in s]
                        #print(f'check12: {check12}')
                        check52 = [s for s in R5min_R2hr_overlap_strs if gene in s and f'-{site}' in s]
                        #print(f'check52: {check52}')

                        # if site has overlap at (R1min-R5min and R1min-R2hr) or (R1min-R5min and R5min-R2hr) then it has overlap at all three timepoints
                        if (len(check15) != 0 and len(check12) != 0) or (len(check15) != 0 and len(check52) != 0) or (len(check12) != 0 and len(check52) != 0):
                            master[gene][site] = ['R1min', 'R5min', 'R2hr']
                        
                        if len(check15) != 0 and len(check12) == 0 and len(check52) == 0:
                            master[gene][site] = ['R1min', 'R5min']

                        if len(check15) == 0 and len(check12) != 0 and len(check52) == 0:
                            master[gene][site] = ['R1min', 'R2hr']
                        
                        if len(check15) == 0 and len(check12) == 0 and len(check52) != 0:
                            master[gene][site] = ['R5min', 'R2hr']
                        
                        # No overlaps detected. Find the timepoint it is add
                        if len(check15) == 0 and len(check12) == 0 and len(check52) == 0:
                            if gene in R1min:
                                check1 = [s for s in R1min[gene] if s == site]
                            else:
                                check1 = []
                            
                            if gene in R5min:
                                check5 = [s for s in R5min[gene] if s == site]
                            else:
                                check5 = []
                            
                            if gene in R2hr:
                                check2 = [s for s in R2hr[gene] if s == site]
                            else:
                                check2 = []

                            # QC there should only be one timepoint the site is present at
                            if len(check1) + len(check5) + len(check2) != 1:
                                raise ValueError(f'There are multiple timepoints this PK cut site is present at but no overlap was detected. Must be an error in the code.')

                            #print(check1, check5, check2)
                            if len(check1) != 0:
                                master[gene][site] = ['R1min']
                            if len(check5) != 0:
                                master[gene][site] = ['R5min']
                            if len(check2) != 0:
                                master[gene][site] = ['R2hr']

                        #print(f'master[gene][site] = {master[gene][site]}')
                        if len(master[gene][site]) == 0:
                            quit()

                master_df = {'gene':[], 'PKsite':[], 'OverlapClass':[]}
                for gene in master:
                    for site in master[gene]:
                        master_df['gene'] += [gene]
                        master_df['PKsite'] += [site]
                        master_df['OverlapClass'] += ['-'.join(master[gene][site])]
                master_df = pd.DataFrame(master_df)
                #print(f'master_df:\n{master_df}')
                return master_df
            ###########################################################################

            # Get initial overlap percentages
            Overlap_df = get_site_overlap(R1min, R5min, R2hr)
            PKsite_overlap = {'R1min':{'boot':[], 'perm':[], 'GT':0, 'count':0}, 
                            'R5min':{'boot':[], 'perm':[], 'GT':0, 'count':0},
                            'R2hr':{'boot':[], 'perm':[], 'GT':0, 'count':0},
                            'R1min-R5min':{'boot':[], 'perm':[], 'GT':0, 'count':0}, 
                            'R1min-R2hr':{'boot':[], 'perm':[], 'GT':0, 'count':0},
                            'R5min-R2hr':{'boot':[], 'perm':[], 'GT':0, 'count':0}, 
                            'R1min-R5min-R2hr':{'boot':[], 'perm':[], 'GT':0, 'count':0}}

            total_count = len(Overlap_df)
            for OC, count in Overlap_df['OverlapClass'].value_counts().items():
                print(OC, count, count/total_count)
                PKsite_overlap[OC]['GT'] = count/total_count
                PKsite_overlap[OC]['count'] = count

            #############################
            ## do permutation testing
            print('PERMUTATION TESTING: PK cutsite overlap')
            total_df = pd.concat([self.inpfiles[buff]['R1min'], self.inpfiles[buff]['R5min'], self.inpfiles[buff]['R2hr']], ignore_index=True)
            #print(f'total_df:\n{total_df}')
            R1min_numcuts = len(self.inpfiles[buff]['R1min'])
            R5min_numcuts = len(self.inpfiles[buff]['R5min'])
            R2hr_numcuts = len(self.inpfiles[buff]['R2hr'])

            for p in tqdm(range(self.reps)):
                
                ptotal_df = total_df.sample(frac=1, ignore_index=True)
                #print(ptotal_df)

                pR1min = get_extended_cutsites(ptotal_df.iloc[ : R1min_numcuts])        
                pR5min = get_extended_cutsites(ptotal_df.iloc[R1min_numcuts : R1min_numcuts + R5min_numcuts])
                pR2hr = get_extended_cutsites(ptotal_df.iloc[R1min_numcuts + R5min_numcuts:])

                pOverlap_df = get_site_overlap(pR1min, pR5min, pR2hr)
                total_count = len(pOverlap_df)
                for OC, count in pOverlap_df['OverlapClass'].value_counts().items():
                    #print(OC, count, count/total_count)
                    PKsite_overlap[OC]['perm'] += [count/total_count]
                
            # calculate the pvalue
            for setstr, set_data in PKsite_overlap.items():
                pvalue = 1
                for perm in set_data['perm']:
                    if perm >= set_data['GT']:
                        pvalue += 1

                pvalue /= (self.reps + 1)
                print(buff, setstr, set_data['GT'], pvalue)
                PKsite_overlap[setstr]['pvalue'] = pvalue

            #############################
            ## do bootstrapping
            # bootstrap the fraction of these classes
            print('BOOTSTRAPPIG BOUNDS: PK cutsite overlap')
            for b in tqdm(range(self.reps)):
                bdf = Overlap_df.sample(frac=1, replace=True)
                bcounts = bdf['OverlapClass'].value_counts()
                bcounts /= total_count
                #print(bcounts)
                for k,v in bcounts.items():
                    PKsite_overlap[k]['boot'] += [v]

            outdf = {'buff':[], 'OverlapClass':[], 'Count':[], 'Overlap':[], 'Bounds':[], 'pvalue':[]}
            for setstr, set_data in PKsite_overlap.items():
                if len(set_data['boot']) != 0:
                    lb = np.percentile(set_data['boot'], 2.5)
                    ub = np.percentile(set_data['boot'], 97.5)
                else:
                    lb = 0.0
                    ub = 0.0
                PKsite_overlap[setstr]['bounds'] = (lb, ub)
                GT = set_data['GT']
                count = set_data['count']
                pvalue = set_data['pvalue']
                
                outdf['buff'] += [buff]
                outdf['OverlapClass'] += [setstr]
                outdf['Count'] += [count]
                outdf['Overlap'] += [f'{GT:.3f}']
                outdf['Bounds'] += [f'({lb:.3f}, {ub:.3f})']
                outdf['pvalue'] += [f'{pvalue:.3e}']

            outdf = pd.DataFrame(outdf)
            print(outdf)

            for setstr, set_data in PKsite_overlap.items():
                print(buff, setstr, set_data['GT'], set_data['bounds'], set_data['pvalue'])

            # Save the dataframe to the data path
            outfile = os.path.join(self.data_path, f'PKsiteOverlap3Table_{self.tag}_{buff}.csv')
            outdf.to_csv(outfile, index=False)
            print(f'SAVED: {outfile}')

            # Create a figure with GridSpec to place the Venn diagram and table in different regions
            fig = plt.figure(figsize=(12, 8))  # Adjust figure size as needed

            # Define a GridSpec with two rows: one for the Venn diagram, one for the table
            gs = fig.add_gridspec(2, 1, height_ratios=[3, 1])  # Venn diagram larger, table smaller

            # Create the Venn diagram in the first row
            ax1 = fig.add_subplot(gs[0])
            #(Abc, aBc, ABc, abC, AbC, aBC, ABC)
            class_order = ['R1min', 'R5min', 'R1min-R5min', 'R2hr', 'R1min-R2hr', 'R5min-R2hr', 'R1min-R5min-R2hr']
            print(class_order)
            counts = [outdf[outdf['OverlapClass'] == c]['Count'].values[0] for c in class_order]
            print(counts)
            overlaps = [outdf[outdf['OverlapClass'] == c]['Overlap'].values[0] for c in class_order]
         
            venn = venn3(counts, ('R1min', 'R5min', 'R2hr'))

            # Update the labels with the percentages
            for idx, label in enumerate(venn.subset_labels):
                if label:  # Check if the label exists
                    size = int(label.get_text())
                    percentage = float(overlaps[idx]) * 100
                    label.set_text(f'{size} ({percentage:.1f}%)')

            # Create the table in the second row
            ax2 = fig.add_subplot(gs[1])
            ax2.axis('off')  # Turn off the axis for the table

            # Add the table
            ax_table = ax2.table(cellText=outdf.values, colLabels=outdf.columns, colWidths=[0.075, 0.2, 0.1, 0.1, 0.2, 0.1], cellLoc='center', loc='center')

            # Increase font size for the table
            ax_table.auto_set_font_size(False)
            ax_table.set_fontsize(8)

            # Adjust layout to avoid overlap
            plt.tight_layout()

            # Set the title for the figure
            plt.suptitle(f'{self.tag} | {buff} | GeneOverlap')

            # Save the figure
            outfile = os.path.join(self.plots_path, f'PKsiteOverlapVenn3Diagram_{self.tag}_{buff}.pdf')
            plt.savefig(outfile)
            print(f'SAVED: {outfile}')

            # Close the figure
            plt.close()

    #######################################################################################
    def PKsite_overlap_2timepoint(self,): 
        offset = 5   
        for buff in ['C', 'CD', 'CG']:
            print(f'Calcutating PKsite overlap between {buff} timepoints')

            # Get R1min cut sites with +/- 3 residue buffer
            def get_extended_cutsites(df):
                out = {}
                for rowi, row in df.iterrows():
                    gene = row['Accession']
                    site = row['site_number']
                    if gene not in out:
                        out[gene] = {}
                    out[gene][site] = []

                    for s in range(row['site_number'] - offset, row['site_number'] + offset + 1):
                        out[gene][site] += [s]
                return out
      

            R5min = get_extended_cutsites(self.inpfiles[buff]['R5min'])
            R5min_size = len(R5min)
            print(f'R5min_size: {R5min_size}')
            
            R2hr = get_extended_cutsites(self.inpfiles[buff]['R2hr'])
            R2hr_size = len(R2hr)
            print(f'R2hr_size: {R2hr_size}')

            # Calculate the total number of unique elements
            total_elements = R5min_size + R2hr_size
            print(f'total_elements: {total_elements}')
            
            ###########################################################################
            ## Get the confidence intervals on the overlap percentages and prob of occuring by random chance
            # Define a general function to do the overlap calcs so i dont have to repeate logic
            def get_site_overlap(R5min, R2hr):
                
                # Define a general function to check overlap of two timepoints
                def check_overlap(t1, t1str, t2, t2str):
                    overlap_strs = []
                    for gene, t1_gene_sites in t1.items():
                        #print('\n', gene, t1_gene_sites)
                        gene_count = 0
                        gene_overlap_strs = []

                        if gene in t2:
                            #print(f'{gene} in second timepoint')
                            t2_gene_sites = t2[gene]
                            #print(t2_gene_sites)

                            # check for atleast 1 overlap 
                            for t1_PKsite, t1_PKsite_ext in t1_gene_sites.items():
                                for t2_PKsite, t2_PKsite_ext in t2_gene_sites.items():
                                    if len(set(t1_PKsite_ext).intersection(t2_PKsite_ext)) != 0:
                                        #print(f'Found overlap between {t1_PKsite} and {t2_PKsite}')
                                        overlap_str = f'{gene}_{t1str}-{t1_PKsite}_{t2str}-{t2_PKsite}'
                                        gene_count += 1
                                        gene_overlap_strs += [overlap_str]
                                        
                            #print(f'{gene} Overlap {t1str}-{t2str} count: {gene_count}\n{gene_overlap_strs}')
                        overlap_strs += gene_overlap_strs
                    return overlap_strs
                    

                # get those sites in R5min and R2hr
                R5min_R2hr_overlap_strs = check_overlap(R5min, 'R5min', R2hr, 'R2hr')
                #print(R5min_R2hr_overlap_strs) 

                # make a master dictionary of all genes and their cutsites across all three timepoints
                master = {}
                for t in [R5min, R2hr]:
                    for gene, gene_sites in t.items():
                        pk_sites = gene_sites.keys()
                        if gene not in master:
                            master[gene] = {}
                        for site in pk_sites:
                            master[gene][site] = []
                
                # For each gene and its sites in the master determine their classifications
                for gene, sites in master.items():
                    for site in sites:
                        #print('\n', gene, site)

                        # check if site had overlap between R1min R5min
                        check52 = [s for s in R5min_R2hr_overlap_strs if gene in s and f'-{site}' in s]
                        #print(f'check52: {check52}')

                        # check for overlaps
                        if len(check52) != 0:
                            master[gene][site] = ['R5min', 'R2hr']
                        
                        # No overlaps detected. Find the timepoint it is add
                        if len(check52) == 0:
                            
                            if gene in R5min:
                                check5 = [s for s in R5min[gene] if s == site]
                            else:
                                check5 = []
                            
                            if gene in R2hr:
                                check2 = [s for s in R2hr[gene] if s == site]
                            else:
                                check2 = []

                            # QC there should only be one timepoint the site is present at
                            if len(check5) + len(check2) != 1:
                                raise ValueError(f'There are multiple timepoints this PK cut site is present at but no overlap was detected. Must be an error in the code.')

                            #print(check5, check2)
                            if len(check5) != 0:
                                master[gene][site] = ['R5min']
                            if len(check2) != 0:
                                master[gene][site] = ['R2hr']

                        #print(f'master[gene][site] = {master[gene][site]}')
                        if len(master[gene][site]) == 0:
                            quit()

                master_df = {'gene':[], 'PKsite':[], 'OverlapClass':[]}
                for gene in master:
                    for site in master[gene]:
                        master_df['gene'] += [gene]
                        master_df['PKsite'] += [site]
                        master_df['OverlapClass'] += ['-'.join(master[gene][site])]
                master_df = pd.DataFrame(master_df)
                #print(f'master_df:\n{master_df}')
                return master_df
            ###########################################################################

            # Get initial overlap percentages
            Overlap_df = get_site_overlap(R5min, R2hr)
            PKsite_overlap = {'R5min':{'boot':[], 'perm':[], 'GT':0, 'count':0},
                            'R2hr':{'boot':[], 'perm':[], 'GT':0, 'count':0},
                            'R5min-R2hr':{'boot':[], 'perm':[], 'GT':0, 'count':0}}

            total_count = len(Overlap_df)
            for OC, count in Overlap_df['OverlapClass'].value_counts().items():
                print(OC, count, count/total_count)
                PKsite_overlap[OC]['GT'] = count/total_count
                PKsite_overlap[OC]['count'] = count

            #############################
            ## do permutation testing
            print('PERMUTATION TESTING: PK cutsite overlap')
            total_df = pd.concat([self.inpfiles[buff]['R5min'], self.inpfiles[buff]['R2hr']], ignore_index=True)
            #print(f'total_df:\n{total_df}')
            R5min_numcuts = len(self.inpfiles[buff]['R5min'])
            R2hr_numcuts = len(self.inpfiles[buff]['R2hr'])

            for p in tqdm(range(self.reps)):
                
                ptotal_df = total_df.sample(frac=1, ignore_index=True)
                #print(ptotal_df)

                pR5min = get_extended_cutsites(ptotal_df.iloc[ : R5min_numcuts])
                pR2hr = get_extended_cutsites(ptotal_df.iloc[R5min_numcuts : ])

                pOverlap_df = get_site_overlap(pR5min, pR2hr)
                total_count = len(pOverlap_df)
                for OC, count in pOverlap_df['OverlapClass'].value_counts().items():
                    #print(OC, count, count/total_count)
                    PKsite_overlap[OC]['perm'] += [count/total_count]
                
            # calculate the pvalue
            for setstr, set_data in PKsite_overlap.items():
                pvalue = 1
                for perm in set_data['perm']:
                    if perm >= set_data['GT']:
                        pvalue += 1

                pvalue /= (self.reps + 1)
                print(buff, setstr, set_data['GT'], pvalue)
                PKsite_overlap[setstr]['pvalue'] = pvalue

            #############################
            ## do bootstrapping
            # bootstrap the fraction of these classes
            print('BOOTSTRAPPIG BOUNDS: PK cutsite overlap')
            for b in tqdm(range(self.reps)):
                bdf = Overlap_df.sample(frac=1, replace=True)
                bcounts = bdf['OverlapClass'].value_counts()
                bcounts /= total_count
                #print(bcounts)
                for k,v in bcounts.items():
                    PKsite_overlap[k]['boot'] += [v]

            outdf = {'buff':[], 'OverlapClass':[], 'Count':[], 'Overlap':[], 'Bounds':[], 'pvalue':[]}
            for setstr, set_data in PKsite_overlap.items():
                if len(set_data['boot']) != 0:
                    lb = np.percentile(set_data['boot'], 2.5)
                    ub = np.percentile(set_data['boot'], 97.5)
                else:
                    lb = 0.0
                    ub = 0.0
                PKsite_overlap[setstr]['bounds'] = (lb, ub)
                GT = set_data['GT']
                count = set_data['count']
                pvalue = set_data['pvalue']
                
                outdf['buff'] += [buff]
                outdf['OverlapClass'] += [setstr]
                outdf['Count'] += [count]
                outdf['Overlap'] += [f'{GT:.3f}']
                outdf['Bounds'] += [f'({lb:.3f}, {ub:.3f})']
                outdf['pvalue'] += [f'{pvalue:.3e}']

            outdf = pd.DataFrame(outdf)
            print(outdf)

            for setstr, set_data in PKsite_overlap.items():
                print(buff, setstr, set_data['GT'], set_data['bounds'], set_data['pvalue'])

            # Save the dataframe to the data path
            outfile = os.path.join(self.data_path, f'PKsiteOverlap2Table_{self.tag}_{buff}.csv')
            outdf.to_csv(outfile, index=False)
            print(f'SAVED: {outfile}')

            
            # Create a figure with GridSpec to place the Venn diagram and table in different regions
            fig = plt.figure(figsize=(12, 8))  # Adjust figure size as needed

            # Define a GridSpec with two rows: one for the Venn diagram, one for the table
            gs = fig.add_gridspec(2, 1, height_ratios=[3, 1])  # Venn diagram larger, table smaller

            # Create the Venn diagram in the first row
            ax1 = fig.add_subplot(gs[0])
            venn = venn2(outdf['Count'].values, ('R5min', 'R2hr'))

            # Update the labels with the percentages
            for idx, label in enumerate(venn.subset_labels):
                if label:  # Check if the label exists
                    size = int(label.get_text())
                    percentage = float(outdf['Overlap'].values[idx]) * 100
                    label.set_text(f'{size} ({percentage:.1f}%)')

            # Create the table in the second row
            ax2 = fig.add_subplot(gs[1])
            ax2.axis('off')  # Turn off the axis for the table

            # Add the table
            ax_table = ax2.table(cellText=outdf.values, colLabels=outdf.columns, colWidths=[0.075, 0.2, 0.1, 0.1, 0.2, 0.1], cellLoc='center', loc='center')

            # Increase font size for the table
            ax_table.auto_set_font_size(False)
            ax_table.set_fontsize(8)

            # Adjust layout to avoid overlap
            plt.tight_layout()

            # Set the title for the figure
            plt.suptitle(f'{self.tag} | {buff} | GeneOverlap')

            # Save the figure
            outfile = os.path.join(self.plots_path, f'PKsiteOverlapVenn2Diagram_{self.tag}_{buff}.pdf')
            plt.savefig(outfile)
            print(f'SAVED: {outfile}')

            # Close the figure
            plt.close()



########################################################################################################
def main():
    """
    Main function to parse arguments and run the DataAnalysis class.
    """

    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-i", "--inpfiles", type=str, required=True, help="Path to FLiPPR processed files")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("-t", "--tag", type=str, required=True, help="Thresholding tag: spa or cov")
    parser.add_argument("-g", "--gene_list", type=str, required=True, help="gene list to mask with")
    args = parser.parse_args()

    inpfiles = args.inpfiles
    outpath = args.outpath
    tag = args.tag
    gene_list = args.gene_list

    Analyzer = DataAnalysis(inpfiles, outpath, tag, gene_list)
    print(Analyzer)

    # threshold SPA data
    #Analyzer.gene_overlap_3timepoint()
    #Analyzer.PKsite_overlap_3timepoint()
    Analyzer.gene_overlap_2timepoint()
    Analyzer.PKsite_overlap_2timepoint()
    print('NORMAL TERMINATION')


if __name__ == "__main__":
    main()

