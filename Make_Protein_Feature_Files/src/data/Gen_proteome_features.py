import os
import math
import requests
from operator import itemgetter
import re
import logging
import argparse
import numpy as np
import pandas as pd
from glob import glob
import mdtraj as md
from Bio.PDB import PDBParser, is_aa
from Bio import PDB
from scipy.spatial.distance import pdist, squareform
np.set_printoptions(linewidth=np.inf, precision=4)
pd.set_option('display.max_rows', None)

class BioDataProcessor:
    """
    Processes biological data including PDB files, sequence data, and interaction potentials.
    """
    #############################################################################################################
    def __init__(self, control_file, uniprot_id, pdb_id, chain, outpath):
        self.uniprot_id = uniprot_id
        self.pdb_id = pdb_id
        self.chain = chain
        self.parse_control_file(control_file)
        self.outpath = outpath
        self.setup_directories(self.outpath)

    #############################################################################################################
    def parse_control_file(self, control_file):
        ### Get control file data
        cntrl_data = {x.split(' = ')[0]:x.split(' = ')[1].strip(' \n') for x in open(control_file,'r').readlines()}
        cntrl_data['gene'] = self.uniprot_id
        cntrl_data['pdb'] = self.pdb_id
        cntrl_data['chain'] = self.chain
        for k,v in cntrl_data.items():
            logging.info(f'{k} = {v}')
            print(f'{k} = {v}')
        self.cntrl_data = cntrl_data
    
    #############################################################################################################
    def setup_directories(self, path_to_outdir):
        """
        Setup output directories for storing results.
        """
        #self.outpath = path_to_outdir
        required_dirs = ['contact_type1_lib', 'contact_type2_lib', 'contact_type3_lib', 'res_features_lib', 'uent_features_lib']
        for directory in required_dirs:
            full_path = os.path.join(path_to_outdir, directory)
            print(full_path)
            if not os.path.exists(full_path):
                os.makedirs(full_path)
                logging.info(f'Created directory: {full_path}')
                print(f'Created directory: {full_path}')

    #############################################################################################################
    def load_contact_potentials(self, path_to_contact_pt):
        """
        Load contact potential data from specified directory.
        """
        self.contact_potentials = {}
        for file_path in glob(path_to_contact_pt):
            name = os.path.basename(file_path).split('_')[0]

            # load potential as a pandas df
            cont_pot = pd.read_csv(file_path)
            cont_pot = cont_pot.rename(columns={"Unnamed: 0":"AA"})
            upper_tri = np.triu(cont_pot.loc[:, cont_pot.columns != 'AA'])

            # Mirror the upper triangle to the lower triangle
            symmetric_df = upper_tri + upper_tri.T - np.diag(np.diag(upper_tri))
            symmetric_df = pd.DataFrame(symmetric_df)
            symmetric_df.columns = cont_pot.loc[:, cont_pot.columns != 'AA'].columns
            symmetric_df['AA'] = cont_pot['AA']

            self.contact_potentials[name] = symmetric_df

            logging.info(f'Loaded contact potential from {file_path}')
            print(f'Loaded contact potential from {file_path}')


    #############################################################################################################
    def get_stride(self, cntrl_data):
        """
        Obtain the STRIDE info on the fly for the PDB file with a call to the cmd line and parse the data into a callable dictionary where the key is a residue and the value is the secondary structure class
        """
        ## Get pdb file name
        if cntrl_data['pdb'] != 'AF':
            self.pdb_file = f'{cntrl_data["pdb_dir"]}{cntrl_data["gene"]}-{cntrl_data["pdb"]}_{cntrl_data["chain"]}.pdb'
        else:
            self.pdb_file = f'{cntrl_data["pdb_dir"]}{cntrl_data["gene"]}.pdb'
        print(f'pdb_file: {self.pdb_file}')
        
        ## Call stride and collect output
        get_elements = itemgetter(*[1,3,6])
        stride_output = [get_elements(x.strip('\n').split()) for x in os.popen(f'stride -o {self.pdb_file}') if x.startswith('LOC') and f' {cntrl_data["chain"]} ' in x]
        self.stride_data = {}
        ## make output dictionary where key is an unmapped resid and the value is the secondary structure element
        for line in stride_output:
            print(line)
            sec = line[0]
            start_PDBres = int(line[1])
            end_PDBres = int(line[2])
            logging.info(f'sec: {sec} | start_PDBres: {start_PDBres} | end_PDBres: {end_PDBres}')
            #print(f'sec: {sec} | start_PDBres: {start_PDBres} | end_PDBres: {end_PDBres}')

            for res in range(start_PDBres, end_PDBres + 1):
                self.stride_data[res] = sec
                #print(f'res: {res} | sec: {sec}')


    #############################################################################################################
    def load_LipMS_data(self, files, uniprotid):
        """
        Load the Limited proteolysis data and select those peptides that have a pvalue below the threshold provided. 
        Use the FDR corrected pvalues if requested.
        """
        lip_data = {}
        lip_resid_df = {}
        for i,f in enumerate(files):
            buff = f.split('/')[-1].split('_')[0]
            timepoint = f.split('/')[-1].split('_')[1]
            print(i, f, buff, timepoint)

            excel_df = pd.read_csv(f)
            df = excel_df[excel_df['Accession'] == uniprotid]
            print(f'Sig LiP-MS cut sites:\n{df}')

            mapped_resids = df['proteinaseKsite'].values
            ratios = df['PeptideRatio1'].values
            
            for mapped_res, ratio in zip(mapped_resids, ratios):
                mapped_resid = re.sub(re.compile(r'\D+'), '', mapped_res)
                #print(mapped_resid, mapped_res, ratio)

                cstr = ",".join([buff, timepoint, mapped_res, f'{ratio:.2f}'])
                if mapped_resid not in lip_resid_df:
                    lip_resid_df[mapped_resid] = [cstr]
                else:
                    lip_resid_df[mapped_resid] += [cstr]

        ## concatenate all cuts if there are multiple per residue
        if len(lip_resid_df) != 0:
            for res, cut_strs in lip_resid_df.items():
                lip_data[int(res)] = "/".join(cut_strs)
        self.lip_data = lip_data
        print(self.lip_data)


    #############################################################################################################
    def get_canonical_sequence_by_accession(self, accession_id):
        """
        Query the uniprot website for the canonical seuqnece for the provdided accession_id.
        If no longer a valid uniprot ID it will kill the script and warn the user.
        """
        # UniProt API URL
        api_url = f"https://www.uniprot.org/uniprot/{accession_id}.json"

        # Make the API request
        response = requests.get(api_url)

        # Check if the request was successful (status code 200)
        if response.status_code == 200:
            # Parse the JSON response
            data = response.json()

            # Check if the entry exists
            if data['primaryAccession'] == accession_id and 'sequence' in data:

                # Extract relevant information

                sequence = data['sequence']['value']
                length = data['sequence']['length']
                self.seq = sequence
                self.length = length

            else:
                print(f"ERROR: No entry found for the accession ID: {accession_id}")
                logging.info(f"ERROR: No entry found for the accession ID: {accession_id}")
                quit()

        else:
            print(f"ERROR: in API request. Status code: {response.status_code}")
            logging.info(f"ERROR: in API request. Status code: {response.status_code}")
            quit()

    #############################################################################################################
    def get_mapping(self, cntrl_data, prot_size):
        """
        Get the mapping from PDB residue ID to Uniprot sequence ID and vis versa. Mappings done by Viraj for JMB paper.
        """
        if cntrl_data['pdb'] != 'AF':
            mapping_file = f'{cntrl_data["mapping_dir"]}{cntrl_data["gene"]}-{cntrl_data["pdb"]}_{cntrl_data["chain"]}_resid_mapping.txt'
            print(f'mapping_file: {mapping_file}')
            logging.info(f'mapping_file: {mapping_file}')

            mapping = np.loadtxt(mapping_file, dtype='O')
            mapping = np.vstack([x[1:] for x in mapping if ('Mapped' in x[0] or 'Modifed_Residue' in x[0] or 'Missense' in x[0])]).astype(int)

            self.mapping_pdb2uniprot = {pdb:uni for pdb, uni in mapping}
            self.mapping_uniprot2pdb = {uni:pdb for pdb, uni in mapping}

        else:
            print(f'Alphafold structures have no mapping. A 1to1 is assumed.')
            logging.info(f'Alphafold structures have no mapping. A 1to1 is assumed.')

            self.mapping_pdb2uniprot = {r:r for r in range(1, prot_size + 1)}
            self.mapping_uniprot2pdb = {r:r for r in range(1, prot_size + 1)}

    #############################################################################################################
    def get_Disprot(self, cntrl_data):
        """
        Get Disprot database information for this gene
        """
        if cntrl_data['path_to_disprot'] != 'None':
            disprot_df = pd.read_csv(cntrl_data['path_to_disprot'], sep='\t')
            disprot_df = disprot_df[['acc', 'start', 'end', 'term_name']]
            disprot_df = disprot_df[disprot_df['term_name'] == 'disorder']
            IDR_dict = {}
            for gene, gene_disprot_df in disprot_df.groupby('acc'):
                if gene == self.uniprot_id:
                    for row_i, row in gene_disprot_df.iterrows():
                        IDR_res = np.arange(row['start'], row['end'] + 1)
                        if gene not in IDR_dict:
                            IDR_dict[gene] = IDR_res
                        else:
                            IDR_dict[gene] = np.unique(np.hstack([IDR_dict[gene], IDR_res]))
        else:
            IDR_dict = {}
            print(f"No IDR file found: {cntrl_data['path_to_disprot']}")
        self.IDR = IDR_dict


    #############################################################################################################
    def get_scop_data(self, cntrl_data):
        """
        Load the SCOP database info for this protein
        #FA-DOMID - SCOP representative family domain identifier
        #FA-PDBID - representative PDB identifier
        #FA-PDBREG - family domain region in PDB resudue numbering
        #FA-UNIID - representative UniProt accession number
        #FA-UNIREG - family domain region in the Uniprot sequence

        #SF-DOMID - SCOP representative superfamily domain identifier
        #SF-PDBID - representative PDB identifier
        #SF-PDBREG - superfamily domain region in PDB resudue numbering
        #SF-UNIID - representative UniProt accession number
        #SF-UNIREG - superfamily domain region in the Uniprot sequence
        #SCOPCLA - SCOP domain classification. The abbreviations denote: TP=protein type, CL=protein class, CF=fold, SF=superfamily, FA=family

        #8024820 1DWK A:1-86   P00816 1-86   8037199 1DWK A:1-86    P00816 1-86    TP=1,CL=1000000,CF=2000087,SF=3000119,FA=4000898
        #8023412 1DWK A:87-156 P00816 87-156 8035792 1DWK A:87-156  P00816 87-156  TP=1,CL=1000003,CF=2000985,SF=3001610,FA=4001450
        """
        if cntrl_data['scope_defs_file'] != 'None':
            scop_defs = {int(x.strip('\n').split(' ')[0]):" ".join(x.strip('\n').split(' ')[1:]) for x in open(cntrl_data['scope_defs_file'], 'r').readlines()}
            scop_data = {'gene':[], 'TP':[], 'CL':[], 'CF':[], 'SF':[], 'FA':[], 'FA_UNIID':[], 'SF_UNIID':[]}
            scop = open(cntrl_data['scope_file'], 'r').readlines()
            for s in scop:
                if not s.startswith('#'):
                    if cntrl_data['gene'] in s:
                        #print(f'{uniprotid} {s}')
                        FA_UNIID = s.split(' ')[4]
                        SF_UNIID = s.split(' ')[9]
                        #print(gene, FA_UNIID, SF_UNIID)
                        
                        #TP=1,CL=1000000,CF=2000087,SF=3000119,FA=4000898
                        TP_id = int(s.split(',')[-5].split(' ')[-1].strip('TP='))
                        CL_id = int(s.split(',')[-4].strip('CL='))
                        CF_id = int(s.split(',')[-3].strip('CF='))
                        SF_id = int(s.split(',')[-2].strip('SF='))
                        FA_id = int(s.split(',')[-1].strip('FA='))

                        TP = scop_defs[TP_id]
                        CL = scop_defs[CL_id]
                        CF = scop_defs[CF_id]
                        SF = scop_defs[SF_id]
                        FA = scop_defs[FA_id]

                        scop_data['gene'] += [cntrl_data['gene']]
                        scop_data['TP'] += [TP]
                        scop_data['CL'] += [CL]
                        scop_data['CF'] += [CF]
                        scop_data['SF'] += [SF]
                        scop_data['SF_UNIID'] += [SF_UNIID]
                        scop_data['FA'] += [FA]
                        scop_data['FA_UNIID'] += [FA_UNIID]

            scop_data = pd.DataFrame(scop_data)
            self.scop_data = scop_data

            scop_resid_data = {}
            for row_i, row in scop_data[['CL', 'FA_UNIID']].iterrows():
                res_class = row['CL']
                res_range = row['FA_UNIID']
                res_range = res_range.split(',')
                res_range = np.hstack([np.arange(int(rr.split('-')[0]), int(rr.split('-')[1]) + 1) for rr in res_range])
                for res in res_range:
                    if res not in scop_resid_data:
                        scop_resid_data[res] = res_class
                    else:
                        scop_resid_data[res] = scop_resid_data[res] + f' ~ {res_class}'
        else:
            scop_resid_data = {}
            print(f"No scop file found: {cntrl_data['scope_defs_file']}")
        self.scop_resid_data = scop_resid_data


    #############################################################################################################
    def get_gene_essentiality(self, cntrl_data): 
        """
        Determine if a gene is essential based on the DEG database MG1655 II entry for E. coli
        """
        if cntrl_data['essential_genes_file'] != 'None':
            deg_data = [x.strip() for x in open(cntrl_data['essential_genes_file'], 'r').readlines() if 'MG1655 II' in x]
            essential_gene = False
            for entry in deg_data:
                if cntrl_data['gene'] in entry:
                    essential_gene = True
        else:
            print(f"No Deg file found file found: {cntrl_data['essential_genes_file']}")
            essential_gene = []

        self.essential_gene = essential_gene


    #############################################################################################################
    def get_AA(self, pdb_file):

        """
        Get the PDB resid to AA mapping for the provided PDB
        """
        three_to_one_letter = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'MSE': 'M', 'PHE': 'F', 
        'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 
        'VAL': 'V'}

        resid2AA = {}
        # Define the path to your PDB file

        # Create a PDB parser
        parser = PDBParser(QUIET=True)

        # Parse the PDB file
        structure = parser.get_structure("protein", pdb_file)

        # Initialize an empty list to store amino acid codes
        amino_acid_codes = []

        # Iterate through the structure and extract amino acid codes
        for model in structure:
            for chain in model:
                for residue in chain:
                    if is_aa(residue):
                        resname = residue.get_resname()
                        resid = residue.get_id()[1]
                        if resname in three_to_one_letter:
                            AA = three_to_one_letter[resname]
                        else:
                            AA = 'NC'
                        #print(resname, resid, AA)
                        resid2AA[resid] = AA
        self.resid2AA = resid2AA

#############################################################################################################
class Analyzer:
    """
    Handles sequence and structural analysis for proteins.
    """
    #############################################################################################################
    def __init__(self, pdb_path, uniprot_id, pdb_id, chain):
        self.pdb_path = pdb_path
        self.traj = md.load(pdb_path) 
        print(f'traj: {self.traj}')
        self.gene = uniprot_id
        self.pdb = pdb_id
        self.chain = chain

    #############################################################################################################
    def get_contacts_type1(self, traj, outpath, gene, pdb, chain, contact_potentials):
        """
        Get and save the residues with CA within 8A of eachother (type1 contacts).
        """
        print('Get and save the residues with CA within 8A of eachother (type1 contacts).')
        ## TEST ####################################################
        # Compute the contact map
        contacts = md.compute_contacts(traj, scheme='ca')

        # Get the contact map and rescale distances from nm to Å (optional)
        contact_map = contacts[0][0]
        #contact_map_distances = contact_map * 10  # convert nm to Å if needed

        # Create a square contact matrix
        n_residues = traj.n_residues
        print(f'n_residues: {n_residues}')
        contact_matrix = np.zeros((n_residues, n_residues))
        for (i, j), distance in zip(contacts[1], contact_map):
            contact_matrix[i, j] = distance
            contact_matrix[j, i] = distance  # Ensure symmetry
        contact_matrix = np.where(contact_matrix == 0, 9999, contact_matrix)
        contact_matrix = np.where(contact_matrix <= 0.8, 1, 0)

        contacts = np.nonzero(contact_matrix)
        contacts = np.stack(contacts).T

        res_idx_2_resid = {res_idx:resid.resSeq for res_idx, resid in enumerate(traj.top.residues)}
        res_idx_2_resname = {res_idx:resname.name for res_idx, resname in enumerate(traj.top.residues)}

        # correct them to resid
        vectorized_func = np.vectorize(lambda x: res_idx_2_resid.get(x, x))
        try:
            contact_resids = vectorized_func(contacts)
            #print(f'contact_resids: {contact_resids} {contact_resids.shape}')
        except:
            print(f'Contact vectorization failed. Maybe due to no contacts present')
            return


        # do the same to get resname
        vectorized_func = np.vectorize(lambda x: res_idx_2_resname.get(x, x))
        contact_resnames = vectorized_func(contacts)
        ## TEST ####################################################

        #print(f'contact_resnames: {contact_resnames} {contact_resnames.shape}')
        # make final contact dataframe and save it
        contact_out_df = {'gene': [], 'pdb': [], 'chain': [], 'pdb_resid_i': [], 'pdb_resid_j': [], 'pdb_resname_i':[], 'pdb_resname_j':[]}
        for contact_idx in range(len(contacts)):

            pdb_resid_i = contact_resids[contact_idx][0]
            pdb_resid_j = contact_resids[contact_idx][1]
            diff = abs(pdb_resid_j - pdb_resid_i)
            if diff < 5:
                continue

            contact_out_df['gene'] += [gene]
            contact_out_df['pdb'] += [pdb]
            contact_out_df['chain'] += [chain]

            pdb_resname_i = contact_resnames[contact_idx][0]
            pdb_resname_j = contact_resnames[contact_idx][1]

            contact_out_df['pdb_resid_i'] += [pdb_resid_i]
            contact_out_df['pdb_resid_j'] += [pdb_resid_j]
            contact_out_df['pdb_resname_i'] += [pdb_resname_i]
            contact_out_df['pdb_resname_j'] += [pdb_resname_j]
        
            for pot, cmap_df in contact_potentials.items():
                if pdb_resname_i in cmap_df.columns and pdb_resname_j in cmap_df.columns:
                    pdb_contact_pot_ij = cmap_df.loc[cmap_df['AA'] == pdb_resname_i, pdb_resname_j].iloc[0]
                    #print(pdb_resid_i, pdb_resname_i, pdb_resid_j, pdb_resname_j, pot, pdb_contact_pot_ij)
                else:
                    pdb_contact_pot_ij = np.nan

                if f'{pot}_cont_pot' not in contact_out_df: 
                    contact_out_df[f'{pot}_cont_pot'] = [pdb_contact_pot_ij]
                else:
                    contact_out_df[f'{pot}_cont_pot'] += [pdb_contact_pot_ij]

        contact_out_df = pd.DataFrame(contact_out_df)
        contact_outfile = f'{outpath}contact_type1_lib/{gene}_{pdb}_{chain}_native_contacts.csv'
        contact_out_df.to_csv(contact_outfile, sep='|', index=False)
        print(f'SAVED: {contact_outfile}')
        logging.info(f'SAVED: {contact_outfile}')


    #############################################################################################################
    def get_contacts_type2(self, traj, outpath, gene, pdb, chain, contact_potentials):
        """
        Get and save the residues with heavy atoms within 4.5A of eachother (type2 contacts).
        """
        print('Get and save the residues with heavy atoms within 4.5A of eachother (type2 contacts).')

        # Compute the contact map
        contacts = md.compute_contacts(traj)

        # Get the contact map and rescale distances from nm to Å (optional)
        contact_map = contacts[0][0]
        #contact_map_distances = contact_map * 10  # convert nm to Å if needed

        # Create a square contact matrix
        n_residues = traj.n_residues
        print(f'n_residues: {n_residues}')
        contact_matrix = np.zeros((n_residues, n_residues))
        for (i, j), distance in zip(contacts[1], contact_map):
            contact_matrix[i, j] = distance
            contact_matrix[j, i] = distance  # Ensure symmetry
        contact_matrix = np.where(contact_matrix == 0, 9999, contact_matrix)
        contact_matrix = np.where(contact_matrix <= 0.45, 1, 0)

        contacts = np.nonzero(contact_matrix)
        contacts = np.stack(contacts).T

        res_idx_2_resid = {res_idx:resid.resSeq for res_idx, resid in enumerate(traj.top.residues)}
        res_idx_2_resname = {res_idx:resname.name for res_idx, resname in enumerate(traj.top.residues)}

        # correct them to resid
        vectorized_func = np.vectorize(lambda x: res_idx_2_resid.get(x, x))
        try:
            contact_resids = vectorized_func(contacts)
            #print(f'contact_resids: {contact_resids} {contact_resids.shape}')
        except:
            print(f'Contact vectorization failed. Maybe due to no contacts present')
            return

        # do the same to get resname
        vectorized_func = np.vectorize(lambda x: res_idx_2_resname.get(x, x))
        contact_resnames = vectorized_func(contacts)
        #print(f'contact_resnames: {contact_resnames} {contact_resnames.shape}')

        # make final contact dataframe and save it
        contact_out_df = {'gene': [], 'pdb': [], 'chain': [], 'pdb_resid_i': [], 'pdb_resid_j': [], 'pdb_resname_i':[], 'pdb_resname_j':[]}
        for contact_idx in range(len(contacts)):
            pdb_resid_i = contact_resids[contact_idx][0]
            pdb_resid_j = contact_resids[contact_idx][1]
            diff = abs(pdb_resid_j - pdb_resid_i)
            if diff < 5:
                continue

            contact_out_df['gene'] += [gene]
            contact_out_df['pdb'] += [pdb]
            contact_out_df['chain'] += [chain]

            pdb_resname_i = contact_resnames[contact_idx][0]
            pdb_resname_j = contact_resnames[contact_idx][1]

            contact_out_df['pdb_resid_i'] += [pdb_resid_i]
            contact_out_df['pdb_resid_j'] += [pdb_resid_j]
            contact_out_df['pdb_resname_i'] += [pdb_resname_i]
            contact_out_df['pdb_resname_j'] += [pdb_resname_j]
        
            for pot, cmap_df in contact_potentials.items():
                if pdb_resname_i in cmap_df.columns and pdb_resname_j in cmap_df.columns:
                    pdb_contact_pot_ij = cmap_df.loc[cmap_df['AA'] == pdb_resname_i, pdb_resname_j].iloc[0]
                    #print(pdb_resid_i, pdb_resname_i, pdb_resid_j, pdb_resname_j, pot, pdb_contact_pot_ij)
                else:
                    pdb_contact_pot_ij = np.nan

                if f'{pot}_cont_pot' not in contact_out_df: 
                    contact_out_df[f'{pot}_cont_pot'] = [pdb_contact_pot_ij]
                else:
                    contact_out_df[f'{pot}_cont_pot'] += [pdb_contact_pot_ij]

        contact_out_df = pd.DataFrame(contact_out_df)
        contact_outfile = f'{outpath}contact_type2_lib/{gene}_{pdb}_{chain}_native_contacts.csv'
        contact_out_df.to_csv(contact_outfile, sep='|', index=False)
        print(f'SAVED: {contact_outfile}')
        logging.info(f'SAVED: {contact_outfile}')

    def get_contacts_type3(self, outpath, contact_potentials):
        """
        Get and save the residues with sidechains with heavy atoms within 4.5A of eachother (type3 contacts).
        This excludes residues with backbone heavy atoms within 8A of eachother. 
        """
        print(f'pdb_path: {self.pdb_path}')

        """Find residues with sidechain heavy atoms within a given distance of each other."""
        pdb_parser = PDB.PDBParser(QUIET=True)
        structure = pdb_parser.get_structure('', self.pdb_path)
        neighbor_search = PDB.NeighborSearch([atom for atom in structure.get_atoms() if atom.element != 'H'])

        close_residues = set()
        backbone_atoms = {"N", "CA", "C", "O"}
        distance_cutoff=4.5
        for model in structure:
            for chain in model:
                for residue in chain:
                    resid = residue.get_id()[1]
                    resname = residue.get_resname()
                    if PDB.is_aa(residue, standard=True):
                        sidechain_atoms = []
                        for atom in residue:
                            atom_name = atom.get_name()
                            #print(atom, atom_name)
                            if atom_name not in backbone_atoms and atom.element != 'H':
                                sidechain_atoms += [atom]
                        #print(f'\nresidue: {residue} | sidechain_atoms: {sidechain_atoms}')
                        #sidechain_atoms = [atom for atom in residue if atom.get_name() not in backbone_atoms and atom.element != 'H']
                        for atom in sidechain_atoms:
                            close_atoms = neighbor_search.search(atom.coord, distance_cutoff, level='A')
                            #print(f'\nclose_atoms: {close_atoms}')
                            for close_atom in close_atoms:
                                #print(f'close_atom: {close_atom} | close_atom.get_parent(): {close_atom.get_parent()} | close_atom.get_name(): {close_atom.get_name()}')
                                if close_atom.get_parent() != residue and close_atom.get_name() not in {"N", "CA", "C", "O"}:
                                    #print(residue, close_atom)
                                    close_atom_resid = close_atom.get_parent().get_id()[1]
                                    close_atom_resname = close_atom.get_parent().get_resname()
                                    diff = abs(resid - close_atom_resid)
                                    #if diff >= 4:
                                        #print(resid, resname, close_atom_resid, close_atom_resname, diff)
                                    #    close_residues.add((resid, resname, close_atom_resid, close_atom_resname))
                                    close_residues.add((resid, resname, close_atom_resid, close_atom_resname))

        contact_out_df = {'gene': [], 'pdb': [], 'chain': [], 'pdb_resid_i': [], 'pdb_resid_j': [], 'pdb_resname_i':[], 'pdb_resname_j':[]}
        for contact in close_residues:
            pdb_resid_i = contact[0]
            pdb_resid_j = contact[2]
            diff = abs(pdb_resid_j - pdb_resid_i)
            if diff < 5:
                continue

            contact_out_df['gene'] += [self.gene]
            contact_out_df['pdb'] += [self.pdb]
            contact_out_df['chain'] += [self.chain]

            pdb_resname_i = contact[1]
            pdb_resname_j = contact[3]

            contact_out_df['pdb_resid_i'] += [pdb_resid_i]
            contact_out_df['pdb_resid_j'] += [pdb_resid_j]
            contact_out_df['pdb_resname_i'] += [pdb_resname_i]
            contact_out_df['pdb_resname_j'] += [pdb_resname_j]
            for pot, cmap_df in contact_potentials.items():
                if pdb_resname_i in cmap_df.columns and pdb_resname_j in cmap_df.columns:
                    pdb_contact_pot_ij = cmap_df.loc[cmap_df['AA'] == pdb_resname_i, pdb_resname_j].iloc[0]
                    #print(pdb_resid_i, pdb_resname_i, pdb_resid_j, pdb_resname_j, pot, pdb_contact_pot_ij)
                else:
                    pdb_contact_pot_ij = np.nan

                if f'{pot}_cont_pot' not in contact_out_df: 
                    contact_out_df[f'{pot}_cont_pot'] = [pdb_contact_pot_ij]
                else:
                    contact_out_df[f'{pot}_cont_pot'] += [pdb_contact_pot_ij]

        contact_out_df = pd.DataFrame(contact_out_df)
        contact_out_df = contact_out_df.sort_values(by=['pdb_resid_i'])
        contact_outfile = f'{outpath}contact_type3_lib/{self.gene}_{self.pdb}_{self.chain}_native_contacts.csv'
        contact_out_df.to_csv(contact_outfile, sep='|', index=False)
        print(f'SAVED: {contact_outfile}')
        logging.info(f'SAVED: {contact_outfile}')

    #############################################################################################################
    def get_uent_features(self, traj, outpath, gene, pdb, chain, cntrl_data, prot_size, mapping_pdb2uniprot):
        """
        Get the features for each unique entanglement provided in the clustered_unampped_GE file
        """

        uent_df = {'PDB':[], 
                    'chain':[], 
                    'ENT-ID':[],
                    'Gn':[],
                    'N_term_thread':[],
                    'Gc':[],
                    'C_term_thread':[],
                    'unmapped-NC':[],
                    'unmapped-NC_wbuff':[],
                    'unmapped-crossings':[], 
                    'unmapped-crossings_wbuff':[], 
                    'loopsize': [], 
                    'num_zipper_nc':[], 
                    'perc_bb_loop':[],
                    'num_loop_contacting_res':[],
                    'num_cross_nearest_neighbors':[],
                    'ent_coverage':[],
                    'min_N_prot_depth_left':[],
                    'min_N_thread_depth_left':[],
                    'min_N_thread_slippage_left':[],
                    'min_C_prot_depth_right':[],
                    'min_C_thread_depth_right':[],
                    'min_C_thread_slippage_right':[], 
                    'prot_size':[], 
                    'ACO':[],
                    'RCO':[]}

        #############################################################################################################################################################################
        ### Load entanglement information if present
        topology = traj.topology

        # Get alpha carbon indices from PDB so I can find contacting residues
        haystack_alpha_carbon_indices = traj.top.select('name CA')
        print(f'haystack_alpha_carbon_indices: {haystack_alpha_carbon_indices}')

        ## load clustered entanglements from summary file
        #../failure_to_form/bioinformatics/get_entanglements/entanglements_0.6g/clustered_unmapped_GE/P0AD61_clustered_GE.txt
        ent_file = f'{cntrl_data["clustered_ent_files"]}{cntrl_data["gene"]}_clustered_GE.txt'
        print(f'ent_file: {ent_file}')

        ## parse lines to get native contacts, crossings,
        ent_present = True
        rbuffer = 3
        pdb_NC_list = [] # list of PDB native contact residues +/- rbuffer
        pdb_NC_core_list = [] # list of PDB natvie contact residues
        pdb_crossing_list = [] # list of PDB crossing residues +/- rbuffer
        pdb_crossing_core_list = [] # list of PDB crossing residues
        total_ent_res = set()
        resid_dict = {}
        if os.path.exists(ent_file):
            ent_data = np.asarray([x.strip('\n') for x in open(ent_file, 'r').readlines()])

            for index, line in enumerate(ent_data):
                if index == 0:
                    continue
                logging.info(f'#######: ENT-ID: {index}')
                ent_core = []

                #['Chain C ', " (23, 315, '+7') ", ' 0.7708439147210744 ', ' -0.18254831113041398', 2]
                #['Chain C', "(23, 315, '+7')", '0.77084', '-0.18255', '10', '23-315;26-322;26-325;26-326;26-329;27-325;28-325;22-322;23-321;23-325']
                line = line.split("|")
                print(line)
                logging.info(line)

                ## check if there was a CCbond. if so skip this ent
                if line[-1] == 'True':
                    print('CCBond detected and ent will be skipped')
                    continue

                ## check that the entanglement isnt in a non-mapped area. if so skip it
                #line = line[1].split(',')
                pdb_NCi_core = line[1].split(',')[0]
                pdb_NCj_core = line[1].split(',')[1]
                pdb_NCi_core = int(float(pdb_NCi_core.replace("(", "").strip()))
                pdb_NCj_core = int(float(pdb_NCj_core.strip()))
                pdb_crossing_res_core = [abs(int(float(re.findall(r'\d+', cross)[0]))) for cross in line[1].split(',')[2:]]
                total_core_ent_res = [pdb_NCi_core, pdb_NCj_core] + pdb_crossing_res_core
                print(pdb_NCi_core, pdb_NCj_core, pdb_crossing_res_core, total_core_ent_res)

                uent_df['PDB'] += [pdb]
                uent_df['chain'] += [chain]
                uent_df['ENT-ID'] += [index]

                #########################################################################
                ## get Gn and Gc and if it is present the cluster size
                if len(line) == 7:
                    num_zipper_nc = int(line[-3])
                else:
                    num_zipper_nc = np.nan
                Gn = float(line[2])
                Gc = float(line[3])
                #print(f'Gn: {Gn} | Gc: {Gc} | num_zipper_nc: {num_zipper_nc}')

                # Calcualte the absolute and relative contact orders
                range_strings = line[-2].split(';')
                
                loops = []
                loop_sizes = []
                # handel weird cases when PDB numbering is negative and regular cases
                for l in range_strings:
                    #print(l)
                    if l.count('-') == 1:
                        loops += [(int(i), int(j)) for i,j in [l.split('-')]]

                    elif l.count('-') == 2 and l.startswith('-'):
                        loops += [(int(i), int(j)) for i,j in [l.rsplit('-', 1)]]

                    elif l.count('-') == 2 and not l.startswith('-'):
                        loops += [(int(i), int(j)) for i,j in [l.split('-', 1)]]

                    elif l.count('-') == 3:
                        hyphen_positions = [i for i, char in enumerate(l) if char == '-']
                        second_hyphen_index = hyphen_positions[1]
                        i = l[:second_hyphen_index]
                        j = l[second_hyphen_index+1:]
                        loops += [(int(i), int(j))]

                for i,j in loops:
                    loop_sizes += [j-i]

                #print(f'loop_sizes: {loop_sizes}')
                num_zipper_nc = len(loop_sizes)
                ACO = np.sum(loop_sizes)/num_zipper_nc
                RCO = ACO/prot_size
                #print(f'Gn: {Gn} | Gc: {Gc} | num_zipper_nc: {num_zipper_nc} | ACO: {ACO} | RCO: {RCO}')
                logging.info(f'Gn: {Gn} | Gc: {Gc} | num_zipper_nc: {num_zipper_nc} | ACO: {ACO} | RCO: {RCO}')
                uent_df['Gn'] += [Gn]
                uent_df['Gc'] += [Gc]
                uent_df['num_zipper_nc'] += [num_zipper_nc]
                uent_df['ACO'] += [ACO]
                uent_df['RCO'] += [RCO]


                #########################################################################
                #get PDB native contact and those +/- rbuffer along the primary structure
                line = line[1].split(',')
                pdb_NCi_core = line[0]
                pdb_NCj_core = line[1]
                pdb_NCi_core = int(float(pdb_NCi_core.replace("(", "").strip()))
                pdb_NCj_core = int(float(pdb_NCj_core.strip()))
                pdb_NC_core = [pdb_NCi_core, pdb_NCj_core]
                pdb_NC_core_list += pdb_NC_core

                pdb_NCi = np.arange(pdb_NCi_core - rbuffer, pdb_NCi_core + rbuffer + 1)
                pdb_NCj = np.arange(pdb_NCj_core - rbuffer, pdb_NCj_core + rbuffer + 1)
                pdb_NC = np.hstack([pdb_NCi, pdb_NCj]).tolist()
                pdb_NC_list += pdb_NC

                logging.info(f'pdb_NC: {pdb_NC}')
                logging.info(f'pdb_NC_core: {pdb_NC_core}')
                uent_df['unmapped-NC'] += [",".join([str(r) for r in pdb_NC_core])]
                uent_df['unmapped-NC_wbuff'] += [",".join([str(r) for r in pdb_NC])]

                loopsize = pdb_NCj_core - pdb_NCi_core
                loop_resids = np.arange(pdb_NCi_core, pdb_NCj_core + 1)
                loop_alpha_carbon_indices = [atom.index for atom in traj.top.atoms if atom.name == 'CA' if atom.residue.resSeq in loop_resids]
                loop_nearest_neighbors_indices = md.compute_neighbors(traj, 0.8, query_indices=loop_alpha_carbon_indices, haystack_indices=haystack_alpha_carbon_indices)[0]
                num_loop_contacting_res = len(loop_nearest_neighbors_indices)

                logging.info(f'num_loop_contacting_res: {num_loop_contacting_res}')
                uent_df['loopsize'] += [loopsize]
                uent_df['perc_bb_loop'] += [loopsize/prot_size]
                uent_df['num_loop_contacting_res'] += [num_loop_contacting_res]
                #########################################################################


                #########################################################################
                #get PDB crossings and those +/- rbuffer along the primary structure
                pdb_crossing_res_core = [abs(int(float(re.findall(r'\d+', cross)[0]))) for cross in line[2:]]
                pdb_crossing_res = np.hstack([np.arange(int(x) - rbuffer, int(x) + rbuffer + 1) for x in pdb_crossing_res_core]).tolist()
                logging.info(f'pdb_crossing_res: {pdb_crossing_res}')
                logging.info(f'pdb_crossing_res_core: {pdb_crossing_res_core}')

                pdb_crossing_list += pdb_crossing_res
                pdb_crossing_core_list += pdb_crossing_res_core
                uent_df['unmapped-crossings'] += [",".join([str(c) for c in pdb_crossing_res_core])]
                uent_df['unmapped-crossings_wbuff'] += [",".join([str(c) for c in pdb_crossing_res])]

                ### Get residues in contact with crossing residues +/- cbuff
                cross_alpha_carbon_indices = [atom.index for atom in traj.top.atoms if atom.name == 'CA' if atom.residue.resSeq in pdb_crossing_res]
                cross_nearest_neighbors_indices = md.compute_neighbors(traj, 0.8, query_indices=cross_alpha_carbon_indices, haystack_indices=haystack_alpha_carbon_indices)[0]
                num_cross_nearest_neighbors = len(cross_nearest_neighbors_indices)

                logging.info(f'num_cross_nearest_neighbors: {num_cross_nearest_neighbors}')
                uent_df['num_cross_nearest_neighbors'] += [num_cross_nearest_neighbors]
                #########################################################################


                #########################################################################
                ## Get number of threads in each termini and depth
                N_term_thread = [c for c in pdb_crossing_res_core if c < pdb_NCi_core]            
                num_N_term_thread = len(N_term_thread)
                C_term_thread = [c for c in pdb_crossing_res_core if c > pdb_NCj_core]            
                num_C_term_thread = len(C_term_thread)
                logging.info(f'N_term_thread: {N_term_thread}')

                logging.info(f'C_term_thread: {C_term_thread}')
                uent_df['N_term_thread'] += [num_N_term_thread]
                uent_df['C_term_thread'] += [num_C_term_thread]

                if num_N_term_thread != 0:
                    min_N_thread_slippage_left = min(N_term_thread)
                    min_N_thread_depth_left = min_N_thread_slippage_left / pdb_NCi_core
                    min_N_prot_depth_left = min_N_thread_slippage_left / prot_size
                else:
                    min_N_thread_slippage_left = np.nan
                    min_N_thread_depth_left = np.nan
                    min_N_prot_depth_left = np.nan
                uent_df['min_N_thread_slippage_left'] += [min_N_thread_slippage_left]
                uent_df['min_N_thread_depth_left'] += [min_N_thread_depth_left]
                uent_df['min_N_prot_depth_left'] += [min_N_prot_depth_left]

                if num_C_term_thread != 0:
                    min_C_thread_slippage_right = prot_size - max(C_term_thread)
                    min_C_thread_depth_right = min_C_thread_slippage_right / (prot_size - pdb_NCj_core)
                    min_C_prot_depth_right = min_C_thread_slippage_right / prot_size
                else:
                    min_C_thread_slippage_right = np.nan
                    min_C_thread_depth_right = np.nan
                    min_C_prot_depth_right = np.nan
                uent_df['min_C_thread_slippage_right'] += [min_C_thread_slippage_right]
                uent_df['min_C_thread_depth_right'] += [min_C_thread_depth_right]
                uent_df['min_C_prot_depth_right'] += [min_C_prot_depth_right]
                #########################################################################
                

                #########################################################################
                ### Get entangled residues. Those that are within 8A of the core residues
                print('Get entangled residues. Those that are within 8A of the core residues')
                ent_core = set(pdb_NC).union(set(pdb_crossing_res))
                ent_core = list(ent_core)
                logging.info(f'ent_core: {ent_core}')
                
                ## Get alpha carbon indices
                haystack_alpha_carbon_indices = traj.top.select('name CA')
                #print(haystack_alpha_carbon_indices)

                ent_res = []
                for ent_core_res in ent_core:
                    alpha_carbon_indices = [atom.index for atom in traj.top.atoms if atom.name == 'CA' and atom.residue.resSeq in [ent_core_res]]

                    # Calculate nearest neighbors within 8 angstroms based on alpha carbon coordinates
                    nearest_neighbors_indices = md.compute_neighbors(traj, 0.8, query_indices=alpha_carbon_indices, haystack_indices=haystack_alpha_carbon_indices, periodic=False)[0]
                    #nearest_neighbors_resid = [traj.top.atom(atom_index).residue.resSeq for atom_index in nearest_neighbors_indices] + ent_core
                    for atom_index in nearest_neighbors_indices:
                        dist = math.sqrt(sum((x - y) ** 2 for x, y in zip(traj.xyz[0][alpha_carbon_indices[0]], traj.xyz[0][atom_index])))
                        #dist = euclidean_distance(traj.xyz[0][alpha_carbon_indices[0]], traj.xyz[0][atom_index])
                        #print(ent_core_res, alpha_carbon_indices, traj.xyz[0][alpha_carbon_indices[0]], traj.top.atom(atom_index).residue.resSeq, atom_index, traj.xyz[0][atom_index], dist)
                        logging.info(f'DIST_check: {ent_core_res}, {alpha_carbon_indices}, {traj.xyz[0][alpha_carbon_indices[0]]}, {traj.top.atom(atom_index).residue.resSeq}, {atom_index}, {traj.xyz[0][atom_index]}, {dist}')
                        if float(f'{dist:.4f}') > 0.8:
                            print(f'ERROR distance greater than 8A.. See logs')
                            quit()
                    nearest_neighbors_resid = [traj.top.atom(atom_index).residue.resSeq for atom_index in nearest_neighbors_indices]
                    nearest_neighbors_resid = list(set(nearest_neighbors_resid))
                    nearest_neighbors_res_idx = [traj.top.atom(atom_index).residue.index for atom_index in nearest_neighbors_indices]
                    num_nearest_neighbors = len(nearest_neighbors_indices)

                    ent_res += nearest_neighbors_resid

        
                ent_res = set(ent_res).union(set(ent_core))
                #print(f'ent_res: {ent_res} {len(ent_res)}')
                logging.info(f'ent_res: {ent_res} {len(ent_res)}')
                uent_df['ent_coverage'] += [len(ent_res)/prot_size]
                uent_df['prot_size'] += [prot_size]
                total_ent_res = set(ent_res).union(total_ent_res)
                for res in ent_res:

                    if res in resid_dict:
                        if isinstance(resid_dict[res], list):
                            resid_dict[res] += [str(index)]

                        else:
                            #print(res, resid_dict[res]) 
                            resid_dict[res] = [str(index)]
                    else:
                        resid_dict[res] = [str(index)]
                
                nonent_res = [traj.top.atom(atom_index).residue.resSeq for atom_index in haystack_alpha_carbon_indices if atom_index not in nearest_neighbors_indices]
                #print(f'nonent_res: {nonent_res}')
                for res in nonent_res:
                    if res in resid_dict:
                        continue
                    else:
                        resid_dict[res] = np.nan
                #########################################################################

        else:
            ent_present = False
            print(f'No entanglement file found')
        self.ent_present = ent_present
        self.resid_dict = resid_dict
        self.pdb_NC_list = pdb_NC_list
        self.pdb_NC_core_list = pdb_NC_core_list
        self.pdb_crossing_list = pdb_crossing_list
        self.pdb_crossing_core_list = pdb_crossing_core_list
        logging.info(f'pdb_NC_list: {pdb_NC_list}')
        logging.info(f'pdb_NC_core_list: {pdb_NC_core_list}')
        logging.info(f'pdb_crossing_list: {pdb_crossing_list}')
        logging.info(f'pdb_crossing_core_list: {pdb_crossing_core_list}')
        logging.info(f'total_ent_res: {total_ent_res} {len(total_ent_res)}')

        ### save file for unique entanglement features
        uent_df = pd.DataFrame(uent_df)
        print(f'uent_df:\n{uent_df}')
        uent_outfile = f'{outpath}uent_features_lib/{gene}_{pdb}_{chain}_native_uents.csv'
        uent_df.to_csv(uent_outfile, sep='|', index=False)
        print(f'SAVED: {uent_outfile}')
        logging.info(f'SAVED: {uent_outfile}')
        ########################################################################################################################

    def parse_ranges(self, range_strings):
        result = []
        for range_str in range_strings:
            # Use regular expression to find the integers in the string
            numbers = re.findall(r'-?\d+', range_str)
            # Convert the found numbers from strings to integers
            int_numbers = [int(num) for num in numbers]
            result.append(int_numbers)
        return result

    #############################################################################################################
    def get_residue_features(self, traj, outpath, gene, pdb, chain, scop_resid_data, cntrl_data, prot_size, mapping_pdb2uniprot, resid2AA, stride_data, IDR_dict, sig_LipMS_df, essential_gene):

        """
        Get the per residue feature file for the desired gene, pdb, and chain.
        The resulting file that is saved has the following columns

        'gene': the uniprot assession ID
        'pdb': the PDB ID
        'chain': the chain ID for the PDB
        'essential': essential gene or not judged by Deg
        'pdb_resid': PDB resid
        'resname': residue name
        'AA': amino acid short name
        'nearest_neighbors': residues with alpha carbons within 8A
        'num_nearest_neighbors': number of residues with alpha carbons within 8A
        'region': entangled region 1 or not 0
        'ent_idx': index for the unique entanglement starting at 1
        'res_sasa': residue solvent accessible surface area
        'median_sasa': residue median solvent accessible surface area
        'NC': loop closing native contact 1 or in buffer 2 or not 0
        'crossing': crossing residue 1 or in buffer 2 or not 0
        'mapped_resid': uniprot mapped resid
        'secondary_struct': STRIDE secondary structure
        'SCOP_class': SCOP class
        'IDR': Disprot IDR True or False
        'cut_str': LiPMS cut sites
        'cut_C_Rall': boolean for any cuts at this residue in Cyto-serum only buffer
        'cut_CD_Rall': boolean for any cuts at this residue in Cyto-serum +DnaK buffer
        'cut_CG_Rall': boolean for any cuts at this residue in Cyto-serum +GroEL buffer

        """

        resFeat_df = {'gene': [],
                'pdb': [],
                'chain': [],
                'uniprot_length':[],
                'essential':[],
                'ent_present':[],
                'pdb_resid': [],
                'resname':[],
                'AA':[],
                'nearest_neighbors': [],
                'num_nearest_neighbors': [],
                'region': [],
                'ent_idx':[],
                'res_sasa':[],
                'median_sasa':[],
                'NC':[],
                'crossing':[],
                'mapped_resid':[],
                'secondary_struct': [],
                'SCOP_class':[],
                'IDR':[],
                'cut_str':[],
                'cut_C_Rall':[],
                'cut_CD_Rall':[],
                'cut_CG_Rall':[]}

        ### Generate residue feature file for PDB
        topology = traj.topology

        # Calculate solvent accessible surface area (SASA) for each residue
        sasa = md.shrake_rupley(traj, mode='residue')[0]

        # Get alpha carbon indices
        haystack_alpha_carbon_indices = traj.top.select('name CA')
        #print(f'haystack_alpha_carbon_indices: {haystack_alpha_carbon_indices}')
        pdb_num_res = len(haystack_alpha_carbon_indices)

        # Iterate over each frame in the trajectory
        for frame in traj:

            # Initialize a list to store nearest neighbors for each residue in the current frame
            frame_neighbors = []

            # Iterate over each residue in the topology
            for residue in frame.top.residues:

                # Skip non-protein residues
                if residue.is_protein:
                    # Get the residue index
                    res_idx = residue.index
                    resid = residue.resSeq
                    #print(f'\n{"#"*50}\n{residue} {res_idx} {resid}')

                    if self.ent_present == True:
                        if resid in self.resid_dict:
                            ent_ids = self.resid_dict[resid]
                        else:
                            ent_ids = np.nan
                            print(f'WARNING: resid not found')

                    else:
                        ent_ids = np.nan

                    # Get the SASA
                    res_sasa = sasa[res_idx]

                    # Get the alpha carbon atom indices for the current residue
                    alpha_carbon_indices = [atom.index for atom in residue.atoms if atom.name == 'CA']

                    # Calculate nearest neighbors within 8 angstroms based on alpha carbon coordinates
                    nearest_neighbors_indices = md.compute_neighbors(traj, 0.8, query_indices=alpha_carbon_indices, haystack_indices=haystack_alpha_carbon_indices, periodic=False)[0]
                    nearest_neighbors_resid = [traj.top.atom(atom_index).residue.resSeq for atom_index in nearest_neighbors_indices]
                    nearest_neighbors_res_idx = [traj.top.atom(atom_index).residue.index for atom_index in nearest_neighbors_indices]
                    num_nearest_neighbors = len(nearest_neighbors_indices)
                    median_sasa = np.median(sasa[nearest_neighbors_res_idx])
                    nearest_neighbors_resid_str = ','.join([str(r) for r in nearest_neighbors_resid])

                    #check for mapping to uniprot sequence
                    if resid in mapping_pdb2uniprot:
                        mapped_resid = int(mapping_pdb2uniprot[resid])
                    else:
                        mapped_resid = 'None'
                    #print(f'mapped_resid: {mapped_resid}')


                    ### check if resid is a native contact NC
                    if resid in self.pdb_NC_list:
                        if resid in self.pdb_NC_core_list:
                            NC = 1
                        else:
                            NC = 2
                    else:
                        NC = 0
                    #print(f'NC status: {NC}')


                    ### check if resid is a crossing residue
                    if resid in self.pdb_crossing_list:
                        if resid in self.pdb_crossing_core_list:
                            crossing = 1
                        else:
                            crossing = 2
                    else:
                        crossing = 0
                    #print(f'crossing status: {crossing}')


                    ### Check what the amino acid is
                    if resid in resid2AA:
                        AA = resid2AA[resid]
                    else:
                        AA = 'None'
                    #print(f'resid: {resid} | AA: {AA}')


                    ### Check if the residue is in a secondary structure
                    if resid in stride_data:
                        res_secondary_struct = stride_data[resid]
                    else:
                        res_secondary_struct = 'None'
                    #print(f'res_secondary_struct: {res_secondary_struct}')


                    ### Check if residue is in an IDR region
                    if gene in IDR_dict:
                        if mapped_resid in IDR_dict[gene]:
                            IDR = True
                        else:
                            IDR = False
                    else:
                        IDR = False
                    #print(f'IDR: {IDR}')


                    ### Check if matpped residue is in scope dicitonary
                    if mapped_resid in scop_resid_data:
                        SCOP_class = scop_resid_data[mapped_resid]
                    else:
                        SCOP_class = 'None'
                    #print(f'SCOP_class: {SCOP_class}')


                    ### Check for cuts
                    #print(sig_LipMS_df)
                    if mapped_resid in sig_LipMS_df:
                        cut_C_Rall = False
                        cut_CD_Rall = False
                        cut_CG_Rall = False
                        cut_str = sig_LipMS_df[mapped_resid]
                        if 'C,' in cut_str:
                            cut_C_Rall = True
                        if 'CD,' in cut_str:
                            cut_CD_Rall = True
                        if 'CG,' in cut_str:
                            cut_CG_Rall = True
                    else:
                        cut_str = 'None'
                        cut_C_Rall = False
                        cut_CD_Rall = False
                        cut_CG_Rall = False
                    #print(f'cut_str: {cut_str}')
                    #print(f'cut_C_Rall: {cut_C_Rall}')
                    #print(f'cut_CD_Rall: {cut_CD_Rall}')
                    #print(f'cut_CG_Rall: {cut_CG_Rall}')


                    ### check resname
                    resname = residue.name
                    #print(f'resid: {resid} |  resname: {resname}')

                    ### add residue entry to output dictionary.
                    resFeat_df['gene'] += [gene]
                    resFeat_df['pdb'] += [pdb]
                    resFeat_df['chain'] += [chain]
                    resFeat_df['uniprot_length'] += [prot_size]
                    #resFeat_df['pdb_coverage'] += [pdb_num_res/prot_size]
                    resFeat_df['essential'] += [essential_gene]
                    resFeat_df['ent_present'] += [self.ent_present]
                    resFeat_df['pdb_resid'] += [resid]
                    resFeat_df['AA'] += [AA]
                    resFeat_df['resname'] += [resname]
                    resFeat_df['nearest_neighbors'] += [nearest_neighbors_resid_str]
                    resFeat_df['num_nearest_neighbors'] += [num_nearest_neighbors]
                    resFeat_df['res_sasa'] += [res_sasa]
                    resFeat_df['median_sasa'] += [median_sasa]
                    resFeat_df['mapped_resid'] += [mapped_resid]
                    resFeat_df['NC'] += [NC]
                    resFeat_df['crossing'] += [crossing]
                    resFeat_df['secondary_struct'] += [res_secondary_struct]
                    resFeat_df['IDR'] += [IDR]
                    resFeat_df['SCOP_class'] += [SCOP_class]
                    resFeat_df['cut_str'] += [cut_str]
                    resFeat_df['cut_C_Rall'] += [cut_C_Rall]
                    resFeat_df['cut_CD_Rall'] += [cut_CD_Rall]
                    resFeat_df['cut_CG_Rall'] += [cut_CG_Rall]

                    if isinstance(ent_ids, list):
                        resFeat_df['region'] += [1]
                        resFeat_df['ent_idx'] += [ent_ids]
                    else:
                        resFeat_df['region'] += [0]
                        resFeat_df['ent_idx'] += ['None']

        #############################################################################################################################################################################

        ### save residue feature file
        resFeat_df = pd.DataFrame(resFeat_df)
        resFeat_df['pdb_coverage'] = len(mapping_pdb2uniprot)/prot_size
        if gene in IDR_dict:
            resFeat_df['unresolved_IDR'] = 1 - (resFeat_df['IDR'].sum()/len(IDR_dict[gene]))
        else:
            resFeat_df['unresolved_IDR'] = np.nan
        #print(resFeat_df)
        #print(resFeat_df[['gene', 'region', 'ent_idx', 'pdb_resid', 'mapped_resid', 'AA', 'NC', 'crossing', 'unresolved_IDR']])
        resFeat_outfile = f'{outpath}res_features_lib/{gene}_{pdb}_{chain}_resfeatures.csv'
        resFeat_df['buried'] = resFeat_df['res_sasa'] < 0.1
        resFeat_df.to_csv(resFeat_outfile, sep='|', index=False)
        print(f'SAVED: {resFeat_outfile}')
        logging.info(f'SAVED: {resFeat_outfile}')


#############################################################################################################
def main():

    # Parse user provided arguments
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-f", "--control_file", required=True, help="Path to control file.")
    parser.add_argument("-o", "--outpath", type=str, required=True, help=f"path to outdir")
    parser.add_argument("-u", "--uniprot_id", type=str, required=True, help='Uniprot Ascession ID')
    parser.add_argument("-p", "--pdb_id", type=str, required=True, help="PDB ID to process")
    parser.add_argument("-c", "--chain", type=str, required=True, help="Chain of the PDB to use")
    parser.add_argument("-l", "--log_file", type=str, required=True, help="Path to logging file")
    parser.add_argument("-m", "--lipms_files", type=str, required=True, help="Path to lipms files to use")
    parser.add_argument("-s", "--skip_contact_lib", type=str, required=True, help="True: skip contact lib generation | False: do not skip")
    args = parser.parse_args()

    # Initiate a Logging file
    logging.basicConfig(filename=args.log_file, level=logging.INFO, format='%(asctime)s %(message)s')
    logging.info("Started processing")

    # Make BioDataProcessor object and get contact potential data, stride data, LiPMS data, canonical sequence, mapping
    data_processor = BioDataProcessor(
            control_file = args.control_file,
            uniprot_id = args.uniprot_id,
            pdb_id = args.pdb_id,
            chain = args.chain,
            outpath = args.outpath)
    #print(data_processor)
    #print(data_processor.cntrl_data)

    # Load in contact potential data
    print(f'\nLoad in contact potential data')
    data_processor.load_contact_potentials(data_processor.cntrl_data['path_to_contact_pt'])

    # Load in STRIDE data for protein
    print(f'\nLoad in STRIDE data for protein')
    data_processor.get_stride(data_processor.cntrl_data)
    #print(data_processor.stride_data)

    # Load in LiPMS data 
    print(f'\nLoad in LiPMS data')
    data_processor.load_LipMS_data(glob(args.lipms_files), data_processor.cntrl_data['gene'])
    #print(data_processor.lip_data)

    # Get uniprot canonical sequence for uniprotID
    print(f'\nGet uniprot canonical sequence for uniprotID')
    data_processor.get_canonical_sequence_by_accession(args.uniprot_id)

    # Get mapping information for gene
    print(f'\nGet mapping information for gene')
    data_processor.get_mapping(data_processor.cntrl_data, data_processor.length)
    #print(data_processor.mapping_pdb2uniprot)

    # Get Disprot data
    print(f'\nGet Disprot data')
    data_processor.get_Disprot(data_processor.cntrl_data)
    #print(data_processor.IDR)

    # Get SCOP data
    print(f'\nGet SCOP data')
    data_processor.get_scop_data(data_processor.cntrl_data)
    #print(data_processor.scop_data)

    # Get if gene is essential or not
    print(f'\nGet if gene is essential or not')
    data_processor.get_gene_essentiality(data_processor.cntrl_data)
    #print(data_processor.essential_gene)

    # Get PDB resid to AA mapping
    print(f'\nGet PDB resid to AA mapping')
    data_processor.get_AA(data_processor.pdb_file)
    #print(data_processor.resid2AA)

    # Initiate the Analyzer object
    print(f'\nInitiate the Analyzer object')
    analyzer = Analyzer(data_processor.pdb_file, args.uniprot_id, args.pdb_id, args.chain)

    # Get unique entanglement features
    print(f'\nGet unique entanglement features')
    analyzer.get_uent_features(analyzer.traj, args.outpath, args.uniprot_id, args.pdb_id, args.chain, data_processor.cntrl_data, data_processor.length, data_processor.mapping_pdb2uniprot)

    # Get residue features
    print(f'\nGet residue features')
    analyzer.get_residue_features(analyzer.traj, 
            args.outpath, 
            args.uniprot_id, 
            args.pdb_id, 
            args.chain, 
            data_processor.scop_resid_data, 
            data_processor.cntrl_data, 
            data_processor.length, 
            data_processor.mapping_pdb2uniprot, 
            data_processor.resid2AA, 
            data_processor.stride_data, 
            data_processor.IDR,
            data_processor.lip_data,
            data_processor.essential_gene)

    if args.skip_contact_lib == 'False':
        # Get and save the native contacts of the PDB file
        print(f'\nGet and save the native contacts of the PDB file')
        analyzer.get_contacts_type1(analyzer.traj, args.outpath, args.uniprot_id, args.pdb_id, args.chain, data_processor.contact_potentials)
        analyzer.get_contacts_type2(analyzer.traj, args.outpath, args.uniprot_id, args.pdb_id, args.chain, data_processor.contact_potentials)
        analyzer.get_contacts_type3(args.outpath, data_processor.contact_potentials)

    logging.info("Finished processing")
    print("Processing completed successfully.")

if __name__ == "__main__":
    main()

