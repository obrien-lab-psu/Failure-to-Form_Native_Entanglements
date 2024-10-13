#!/usr/bin/env python3
import requests, logging, os, sys
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.PDB import PDBParser, PDBIO, Select, Polypeptide
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.PDB.PDBList import PDBList
from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time
import argparse
import pandas as pd
from modeller import *
from modeller.automodel import *
from modeller import restraints, physical
import numpy as np
import glob

pd.set_option('display.max_rows', 500)

class Grabber:
    """
    A class to handel downloading PDBs, extracting the specified chains, finding missing resiudes, and downloading the FASTA file
    """

    def __init__(self, args):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
            ("--candidates", type=str, required=True, help="Path to  candidates file")
            ("--outpath", type=str, required=True, help="Path to output directory")
            ("--log", type=str, required=True, help="Path to logging file")
            ("--cath", type=str, requried=True, help="Path to latest cath-domain-boundaries-seqreschopping.txt file")
        """

        # parse the parameters 
        self.outpath = args.outpath
        self.candidates = pd.read_csv(args.candidates)
        self.log = args.log

        # parse the cath file
        self.cath_file = args.cath
        self.cath_data = np.loadtxt(self.cath_file, dtype=str)
        print(self.cath_data)
        df = {'PDB':[], 'CHAIN':[], 'CLASS':[], 'CLASS_ID':[], 'PDB_resid_start':[], 'PDB_resid_end':[]}
        for idstr, rangstr in self.cath_data:
            pdb = idstr[:4]
            chain = idstr[4]
            domain = idstr[6]
            #print(idstr, pdb, chain, domain)
            
            if len(df['PDB']) == 0:
                domainID = 0
            elif pdb != df['PDB'][-1]:
                domainID = 0
            else:
                domainID += 1

            for ri,r in enumerate(rangstr.split(',')):
                start = r.split('-')[0]
                end = r.split('-')[1]
                start, end = int(start), int(end)
                #print(idstr, ri, r, start, end)

                df['PDB'] += [pdb]
                df['CHAIN'] += [chain]
                if domain == '0':
                    df['CLASS'] += ['c']
                elif domain == '1':
                    df['CLASS'] += ['b']
                elif domain == '2':
                    df['CLASS'] += ['a']
                else:
                    df['CLASS'] += ['n']
                df['CLASS_ID'] += [domainID]
                df['PDB_resid_start'] += [start]
                df['PDB_resid_end'] += [end]

        self.cath_data = pd.DataFrame(df)
        print(self.cath_data)

        logging.info(f'candidates:\n{self.candidates}')

        # make PDBs dir, FASTA dir, and MISSING dir
        self.PDBspath = os.path.join(self.outpath, 'PDBs')
        if not os.path.exists(self.PDBspath):
            os.makedirs(self.PDBspath)
            print(f'Made directory: {self.PDBspath}')

        self.FASTApath = os.path.join(self.outpath, 'FASTA')
        if not os.path.exists(self.FASTApath):
            os.makedirs(self.FASTApath)
            print(f'Made directory: {self.FASTApath}')

        self.MISSINGpath = os.path.join(self.outpath, 'MISSING')
        if not os.path.exists(self.MISSINGpath):
            os.makedirs(self.MISSINGpath)
            print(f'Made directory: {self.MISSINGpath}')

        self.DOMAINSpath = os.path.join(self.outpath, 'DOMAINS')
        if not os.path.exists(self.DOMAINSpath):
            os.makedirs(self.DOMAINSpath)
            print(f'Made directory: {self.DOMAINSpath}')
        
        self.Modellerpath = os.path.join(self.outpath, 'MODLLER')
        if not os.path.exists(self.Modellerpath):
            os.makedirs(self.Modellerpath)
            print(f'Made directory: {self.Modellerpath}')

    ### Gather the data from the interwebs
    def run(self,):
        arr = np.asarray(self.candidates[['gene', 'pdb', 'chain']].values, dtype=str)
        for rowi, (gene, pdb, chain) in enumerate(np.unique(arr, axis=0)):
            logging.info(f'{"#"*20} {rowi} {gene} {pdb} {chain} {"#"*20}')
            print(f'{"#"*20} {rowi} {gene} {pdb} {chain} {"#"*20}')
            output_pdb_file = os.path.join(self.PDBspath, f'{gene}_{pdb}_{chain}.pdb')
            output_fasta_file = os.path.join(self.FASTApath, f'{gene}_{pdb}_{chain}.fasta')
            output_missing_file = os.path.join(self.MISSINGpath, f'{gene}_{pdb}_{chain}.missing')
            output_rebuilt_pdb_file = os.path.join(self.PDBspath, f'{gene}_{pdb}_{chain}_rebuilt.pdb')
            output_domain_file = os.path.join(self.DOMAINSpath, f'{gene}_{pdb}_{chain}.txt')
            working_dir = os.path.join(self.Modellerpath, f'{gene}_{pdb}_{chain}/')
            if not os.path.exists(working_dir):
                os.makedirs(working_dir)
                print(f'Made directory: {working_dir}')

            ## processes the PDB and find out informatio 
            fasta_sequence, align_residues_df = self.process_pdb_and_sequence(gene, pdb, chain, output_pdb_file=output_pdb_file, output_fasta_file=output_fasta_file, output_missing_file=output_missing_file, output_domain_file=output_domain_file)
            residues_to_remodell = align_residues_df[(align_residues_df['MISSING'] == True) | (align_residues_df['MUTATION'] == True)]
            logging.info(f'residues_to_remodell:\n{residues_to_remodell}')

            if len(residues_to_remodell) != 0:
                rebuild_missing_residues(output_pdb_file, fasta_sequence, residues_to_remodell, output_rebuilt_pdb_file, working_dir)
            else:
                logging.info(f'No missing residues found. No need to rebuild {gene} {pdb} {chain}')
                # Create a PDB parser
                parser = PDBParser(QUIET=True)
                
                # Parse the structure
                structure = parser.get_structure('structure', output_pdb_file)

                # Save the modified structure using PDBIO
                io = PDBIO()
                io.set_structure(structure)
                io.save(output_rebuilt_pdb_file)

                logging.info(f"Saved rebuilt PDB {output_pdb_file}")   

    ### Extract PDB chain and save it, get fasta sequence and save it, find missing residues and save it
    def process_pdb_and_sequence(self, uniprot_id, pdb_id, chain_id, output_pdb_file='extracted_chain.pdb', output_fasta_file='extracted_chain.fasta', output_missing_file='extracted_chain.missing', output_domain_file='domain.txt'):
        ##################################################################################################
        # Step 1: Download the PDB file
        pdbl = PDBList()
        pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir=self.PDBspath, file_format='pdb')
        
        # remove obsolete directory made by retreival
        os.rmdir('obsolete')
        ##################################################################################################

        ##################################################################################################      
        # Step 2: Parse the PDB file and extract the specified chain
        parser = PDBParser()
        structure = parser.get_structure(pdb_id, pdb_file)
        
        # Find the specified chain
        chain = None
        for model in structure:
            if chain_id in model:
                chain = model[chain_id]
                break
        
        if not chain:
            raise ValueError(f"Chain {chain_id} not found in PDB {pdb_id}.")
        
        io = PDBIO()
        io.set_structure(structure)
        ##################################################################################################   

        ##################################################################################################
        # Step 3: Extract the amino acid sequence of the specified chain
        ppb = Polypeptide.PPBuilder()
        polypeptides = ppb.build_peptides(chain)
        PDB_seq = []
        residues_info = []
        for pp in polypeptides:
            #PDB_seq += [str(pp.get_sequence())]
            for residue in pp:
                # Get the amino acid sequence (1-letter code) for the residue
                aa = Polypeptide.three_to_one(residue.get_resname())
                
                # Get the residue ID (resid)
                resid = residue.id[1]  # The numerical part of the residue ID
                
                # Append to lists
                PDB_seq.append(aa)
                residues_info.append((aa, resid))

        #PDB_seq = polypeptides[0].get_sequence()
        PDB_seq = ''.join(PDB_seq)
        logging.info(polypeptides)
        logging.info(f'{uniprot_id} {pdb_id} {chain_id } PDB_seq: {PDB_seq} {len(PDB_seq)}')
        print(f'{uniprot_id} {pdb_id} {chain_id } PDB_seq: {PDB_seq} {len(PDB_seq)}')

        ## get map of PDB sequence index to amino acid and resid 
        pdb_map = {i:r for i,r in enumerate(residues_info)}
        logging.info(f'pdb_map: {pdb_map}')
        print(f'pdb_map: {pdb_map}')
        ################################################################################################## 

        ##################################################################################################
        # Step 4: Download the canonical FASTA sequence from UniProt
        fasta_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
        response = requests.get(fasta_url)
        
        if response.status_code == 200:
            fasta_file = f"{uniprot_id}_canonical.fasta"
            with open(output_fasta_file, "w") as file:
                file.write(response.text)
            logging.info(f"Canonical FASTA sequence saved to {output_fasta_file}")
            
            # Optionally, parse the sequence and print it
            fasta_record = SeqIO.read(output_fasta_file, "fasta")
            logging.info(f"FASTA Protein sequence for {uniprot_id}:\n{fasta_record.seq}")
            print(f"FASTA Protein sequence for {uniprot_id}:\n{fasta_record.seq}")
        else:
            logging.info(f"Failed to download FASTA file for UniProt ID {uniprot_id}. Status code: {response.status_code}")
        ##################################################################################################

        ##################################################################################################
        # Step 5: align the two sequences and find the missing residues
        align_residues_df = self.find_missing_residues(PDB_seq, fasta_record.seq, pdb_map)
        align_residues_df['uniprot_id'] = uniprot_id
        align_residues_df['pdb_id'] = pdb_id
        align_residues_df['chain_id'] = chain_id
        align_residues_df['CLASS'] = ''
        align_residues_df['CLASS_ID'] = ''
        print(align_residues_df)
        ##################################################################################################

        ##################################################################################################
        # Step 6: make resid map of current PDB resid to updated PDB resid that should match the FAST indexing starting at 1 since Modller hates resides <= 0
        non_missing_mutant_residues_df = align_residues_df[(align_residues_df['MISSING'] == False) & (align_residues_df['MUTATION'] == False)]
        #print(f'non_missing_mutant_residues_df:\n{non_missing_mutant_residues_df}')
        resid_mapping = {int(p):int(f) for f,p in non_missing_mutant_residues_df[['FASTA_resid', 'PDB_resid']].values}
        logging.info(f'resid_mapping: {resid_mapping}')
        print(f'resid_mapping: {resid_mapping}')
        ##################################################################################################

        ##################################################################################################
        # Step 7: get those resid that were mutated to remove them from the PDB before saving 
        # also get those that did not align and remove them as they are likely fringe sequences we do not care about
        #extra_resids = np.asarray(align_residues_df[align_residues_df['MUTATION'] == True]['PDB_resid'].values, dtype=int)
        extra_resids = [p for k,(a,p) in pdb_map.items() if int(p) not in resid_mapping.keys()]
        logging.info(f'extra_resids to remove before rebuild: {extra_resids}')
        print(f'extra_resids to remove before rebuild: {extra_resids}')
        ##################################################################################################

        ##################################################################################################
        # Step 8: check if there is CATH information and save the domain file
        if pdb_id.lower() in self.cath_data['PDB'].values:
            logging.info(f'CATH DATA FOUND')
            print(f'CATH DATA FOUND')

            # get the cath data for this specific pdb
            pdb_cath_data = self.cath_data[(self.cath_data['PDB'] == pdb_id.lower()) & (self.cath_data['CHAIN'] == chain_id)]
            logging.info(f'\n{pdb_cath_data}')
            print(f'{pdb_cath_data}')
  
            class_dict = {}
            for rowi, row in pdb_cath_data.iterrows():
                logging.info(f'Getting residue ranges for cath domain row {rowi}')

                # get the start of the domain. If it is not inside the mapped area make the first residue the start
                start = align_residues_df[align_residues_df['PDB_resid'] == row["PDB_resid_start"]]['PDB_resid'].values
                if len(start) != 0:
                    start = start[0]
                else:
                    logging.info(f'Start residue for this cath domain was not found and is outside the aligned region. The first residue in the protein will be used')
                    start = align_residues_df[~align_residues_df['PDB_resid'].isna()]['PDB_resid'].values[0]

                if row["PDB_resid_end"] in align_residues_df['PDB_resid'].values:
                    end = align_residues_df[align_residues_df['PDB_resid'] == row["PDB_resid_end"]]['PDB_resid'].values[0]
                else:
                    end = max(align_residues_df['PDB_resid'].values)
                
                start, end = int(start), int(end)
                c = row['CLASS']
                class_idx  = row['CLASS_ID']

                logging.info(f'start: {start}, end: {end} {c} {class_idx}')

                # add residues to class dictionary to track all classes associated with this protein
                if class_idx not in class_dict:
                    class_dict[class_idx] = {}
                
                for r in np.arange(start, end + 1):
                    if r not in class_dict[class_idx]:
                        class_dict[class_idx][r] = c
                    else:
                        raise ValueError(f'Residue: {r} already has a class! {class_dict[r]}')

                # update class index counter
                class_idx += 1
        
        else:
            # if no cath data found check secondary structure content via stride and classify it that way
            logging.info(f'NO CATH DATA FOUND! GUESSING for {output_domain_file}')
            class_dict = {0:{}}

            ## launch stride analysis and parse the results
            stride_cmd = f'stride -o {pdb_file} -r{chain_id} -c{chain_id}'
            logging.info(f'stride_cmd: {stride_cmd}')
            stride_results = os.popen(stride_cmd)
            stride_df = {'secondary_class':[], 'start':[], 'end':[], 'delta':[]}
            stride_results = [x.split()[:7] for x in stride_results if 'LOC' in x]
            total_res = 0
            for _, sec, _, start, _, _, end, in stride_results:
                delta = int(end) - int(start)
                total_res += delta#
                stride_df['secondary_class'] += [sec]
                stride_df['start'] += [start]
                stride_df['end'] += [end]
                stride_df['delta'] += [delta]
            stride_df = pd.DataFrame(stride_df)
            stride_df['perc'] = stride_df['delta']/total_res
            logging.info(f'stride_df:\n{stride_df}')
            perc_alpha = np.sum(stride_df[stride_df['secondary_class'].str.contains('AlphaHelix')]['perc'].values)
            perc_strand = np.sum(stride_df[stride_df['secondary_class'].str.contains('Strand')]['perc'].values)
            logging.info(f'perc_alpha: {perc_alpha}')
            logging.info(f'perc_strand: {perc_strand}')

            if perc_alpha >= 0.2 and perc_strand >= 0.2:
                c = 'c'
            elif perc_alpha >= 0.2 and perc_strand < 0.2:
                c = 'a'
            elif perc_alpha < 0.2 and perc_strand >= 0.2:
                c = 'b'
            else:
                raise ValueError(f'There was no dominate alpha helical or beta strand or mixed behaviour!')
            print(f'Secondary structure class from stride analysis: {c}')
           
        ##################################################################################################

        ##################################################################################################
        # Step 9: update align_residues_df with class data
        # if only 1 class is found make the whole PDB that class
        # if more than 1 separate into the classes and find any residues that did not make it into a class
        if len(class_dict) == 1:
            align_residues_df['CLASS_ID'] = 0
            align_residues_df['CLASS'] = c
        else:
            for classID, classData in class_dict.items():
                for r,c in classData.items():
                    align_residues_df.loc[(align_residues_df['PDB_resid'] == float(r)), 'CLASS_ID'] = classID
                    align_residues_df.loc[(align_residues_df['PDB_resid'] == float(r)), 'CLASS'] = c
            
            # check if there are any missing residues without a class
            if len(align_residues_df[align_residues_df['CLASS_ID'] == '']) != 0:
                print(align_residues_df[align_residues_df['CLASS_ID'] == ''])
                print(f'ERROR: some residues are not in a class')

                # First for each row without an assignement try scanning forward to determine the first assignment and then use that 
                for f_res in align_residues_df[align_residues_df['CLASS_ID'] == '']['FASTA_resid'].values:
                    temp = align_residues_df[(align_residues_df['FASTA_resid'] >= f_res) & (align_residues_df['CLASS_ID'] != '')]
                    if len(temp) != 0:
                        new_class = temp['CLASS'].values[0]
                        new_classID = temp['CLASS_ID'].values[0]
                        print(f'FASTA_resid: {f_res} new_class: {new_class} new_classID: {new_classID}')
                        align_residues_df.loc[(align_residues_df['FASTA_resid'] == f_res), 'CLASS'] = new_class
                        align_residues_df.loc[(align_residues_df['FASTA_resid'] == f_res), 'CLASS_ID'] = new_classID
                print(align_residues_df[align_residues_df['CLASS_ID'] == '']) 

                # last for each row without an assignement left over try scanning backwards to determine the last previous assignemnt
                for f_res in align_residues_df[align_residues_df['CLASS_ID'] == '']['FASTA_resid'].values:
                    temp = align_residues_df[(align_residues_df['FASTA_resid'] <= f_res) & (align_residues_df['CLASS_ID'] != '')]
                    if len(temp) != 0:
                        new_class = temp['CLASS'].values[-1]
                        new_classID = temp['CLASS_ID'].values[-1]
                        print(f'FASTA_resid: {f_res} new_class: {new_class} new_classID: {new_classID}')
                        align_residues_df.loc[(align_residues_df['FASTA_resid'] == f_res), 'CLASS'] = new_class
                        align_residues_df.loc[(align_residues_df['FASTA_resid'] == f_res), 'CLASS_ID'] = new_classID
                
                # final quality check to ensure there are no unclassed residues
                if len(align_residues_df[align_residues_df['CLASS_ID'] == '']) != 0:
                    raise ValueError(f'After attempting to readjust domain mapping there is still unmapped residues:\n{align_residues_df[align_residues_df["CLASS_ID"] == ""]}')

        
        # Step 10: get the domain data for output
        outstrs = []
        for domainID, domain in align_residues_df.groupby('CLASS_ID'):
            c = domain['CLASS'].values[0]
            resids = domain['FASTA_resid'].values
            out = find_continuous_ranges(resids)
            outstrs += [out + f' {c}']
        logging.info(outstrs)
        np.savetxt(output_domain_file, outstrs, fmt='%s')
        logging.info(f'SAVED: {output_domain_file}')
        ##################################################################################################

        ##################################################################################################
        # Step 11: save align_residues_df that has the PDB aligned to the FASTA seq and all its info
        logging.info(f'align_residues_df:\n{align_residues_df}')
        print(f'align_residues_df:\n{align_residues_df}')
        align_residues_df.to_csv(output_missing_file, index=False)
        logging.info(f'SAVED: {output_missing_file}')
        ##################################################################################################

        ##################################################################################################
        # Step 12: define class to save only those portions of the PDB I want to modell
        # Define the ChainSelect class with renumbering support
        class ChainSelect(Select):
            def __init__(self, chain_id, exclude_resids, resid_mapping):
                self.chain_id = chain_id
                self.exclude_resids = exclude_resids
                self.resid_mapping = resid_mapping
                self.selected_resids = set()  # To track selected residues

            # Accept only the specified chain
            def accept_chain(self, chain):
                return chain.id == self.chain_id

            # Exclude specified residues and handle duplicates
            def accept_residue(self, residue):
                resid = residue.get_id()[1]  # Get the current residue ID (resid)

                # Exclude residues that should be ignored
                if resid in self.exclude_resids:
                    return False

                # Check if this resid has already been selected (to handle duplicates)
                if resid in self.selected_resids:
                    return False  # Skip duplicates

                # Mark the resid as selected
                self.selected_resids.add(resid)

                return True

            # Renumber residues in reverse order to prevent conflicts
            def renumber_residues_reverse(self, structure):
                for chain in structure.get_chains():
                    # Get list of residues sorted in reverse order of current resid
                    residues = list(chain.get_residues())
                    residues.sort(key=lambda r: r.get_id()[1], reverse=True)

                    # Apply renumbering from highest to lowest
                    for residue in residues:
                        resid = residue.get_id()[1]
                        if resid in self.resid_mapping:
                            new_resid = self.resid_mapping[resid]
                            residue.id = (residue.id[0], new_resid, residue.id[2])

            # Accept only ATOM records (exclude HETATM)
            def accept_atom(self, atom):
                return atom.get_parent().get_id()[0] == ' '

        io.save(output_pdb_file, select=ChainSelect(chain_id=chain_id, exclude_resids=extra_resids, resid_mapping=resid_mapping))
        logging.info(f"Extracted chain {chain_id} saved to {output_pdb_file}")
        ##################################################################################################

        return fasta_record.seq, align_residues_df

    ### Align two PDB sequenes and find missing residues
    def find_missing_residues(self, pdb_seq, canonical_seq, pdb_map):

        # Step 1: Create FASTA files for sequence1 and sequence2
        with open("seq1.fasta", "w") as seq1_file:
            SeqIO.write(SeqRecord(Seq(pdb_seq), id="seq1", description="Query sequence"), seq1_file, "fasta")
        
        with open("seq2.fasta", "w") as seq2_file:
            SeqIO.write(SeqRecord(Seq(canonical_seq), id="seq2", description="Subject sequence"), seq2_file, "fasta")
        
        # Step 2: Create a BLAST database from sequence2 (the subject sequence)
        logging.info("Creating BLAST database from sequence 2...")
        os.system("makeblastdb -in seq2.fasta -dbtype prot")

        # Step 3: Run BLAST (sequence1 against sequence2)
        logging.info("Running BLAST...")
        blastp_cline = NcbiblastpCommandline(query="seq1.fasta", db="seq2.fasta", outfmt=5, out="blast_results.xml")
        stdout, stderr = blastp_cline()

        # Step 4: Parse the BLAST results and track the best alignment
        best_hsp = None
        best_alignment = None
        best_score = -float("inf")  # Initialize with the lowest possible score

        with open("blast_results.xml") as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            
            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        if hsp.score > best_score:
                            best_hsp = hsp
                            best_alignment = alignment
                            best_score = hsp.score

        # Step 5: Output the best alignment with missing residues adjustment
        if best_hsp and best_alignment:
            logging.info("****Best Alignment****")
            logging.info(f"Score: {best_hsp.score}")
            logging.info(f"E-value: {best_hsp.expect}")

            # Adjust the output to show missing residues
            subject_start = best_hsp.sbjct_start
            logging.info(f'subject_start: {subject_start}')
            subject_end = best_hsp.sbjct_end
            logging.info(f'subject_end: {subject_end}')
            query_start = best_hsp.query_start
            logging.info(f'query_start: {query_start}')
            query_end = best_hsp.query_end
            logging.info(f'query_end: {query_end}')
        

            # Find the number of missing residues at the beginning of the query sequence
            front_missing_residues = subject_start - 1
            logging.info(f'front_missing_residues: {front_missing_residues}')
            back_missing_residues = len(canonical_seq) - subject_end
            logging.info(f'back_missing_residues: {back_missing_residues}')
            
            # Create the adjusted alignment with missing residues
            query_alignment = '-' * front_missing_residues + best_hsp.query + '-' * back_missing_residues
            match_alignment = '-' * front_missing_residues + best_hsp.match + '-' * back_missing_residues

            #subject_alignment = best_hsp.sbjct
            subject_alignment = canonical_seq

            logging.info(f"Query: {query_alignment}")
            logging.info(f"Match: {match_alignment}")
            logging.info(f"Subje: {subject_alignment} {len(subject_alignment)}")
         
        
        pdb_idx = query_start - 1 # start index in 0 for the PDB residue in the matched sequence
        ## determine if there is a missing tail offset
        align_df = {'FASTA_resid':[], 'FASTA':[], 'PDB':[], 'PDB_resid':[], 'MISSING':[], 'MUTATION':[]}
        for f_i, (f_res, p_res) in enumerate(zip(subject_alignment, query_alignment)):
            # determine if there was a missing residue
            if p_res == '-':
                missing = True
            else:
                missing = False
            
            # determine if there was a mutation
            if f_res != p_res and p_res != '-':
                mutation = True
            else:
                mutation = False

            # get PDB idx
            if p_res != '-':
                p_i = pdb_map[pdb_idx][1]
                pdb_idx += 1
            else:
                p_i = np.nan

            #print(f_i + 1, f_res, p_i, p_res, missing, mutation)
            align_df['FASTA_resid'] += [f_i + 1]
            align_df['FASTA'] += [f_res]
            align_df['PDB'] += [p_res]
            align_df['PDB_resid'] += [p_i]
            align_df['MISSING'] += [missing]
            align_df['MUTATION'] += [mutation]
        align_df = pd.DataFrame(align_df)

        ## clean up files from blast analysis so they dont affect next analysis and clutter dir
        os.remove('blast_results.xml')
        seq_files = glob.glob('./seq*')
        logging.info(seq_files)
        for f in seq_files:
            os.remove(f)
        logging.info(f'Removed junk seq files and blast results')
        return align_df

### Rebuild missing residues with modeller 
def rebuild_missing_residues(pdb_path, fasta_sequence, missing_residues_df, output_rebuilt_pdb_file, working_dir):
    logging.info(f'Rebuilding missing residues')

    # Ensure the working directory exists
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    
    # change to working dir
    os.chdir(working_dir)


    # Initialize the MODELLER environment
    env = Environ()
    env.io.output_directory = working_dir

    # Align the target sequence (with missing residues) to the existing structure
    aln = Alignment(env)
    mdl = Model(env, file=pdb_path)  # Load the existing PDB structure
    aln.append_model(mdl, align_codes='existing', atom_files=pdb_path)
    
    # Append the complete sequence to the alignment
    fasta_sequence = str(fasta_sequence)
    logging.info(fasta_sequence)
    aln.append_sequence(fasta_sequence)

    target_code = 'target'
    aln[-1].code = target_code  # Assign the alignment code for the target sequence
    aln.align2d()
    
    # Save the alignment file
    #alignment_file = 'alignment.ali'
    alignment_file = os.path.join(working_dir, 'alignment.ali')
    aln.write(file=alignment_file)
    
    # Define a custom model class to restrain the non-missing residues
    class MyModel(automodel):
        def select_atoms(self):
            # Select only missing residues for refinement
            # Select only the regions around existing residues for refinement
            selection = []
            selected_residues = []
            for _, row in missing_residues_df.iterrows():
                res_id = int(row['FASTA_resid'])
                chain = 'A'
                # Add the missing residue position as a selection target
                try:
                    #selection.append(self.residue_range(f"{res_id}:{chain}", f"{res_id}:{chain}"))
                    if res_id - 1 < 1:
                        start = res_id
                    else:
                        start = res_id - 1

                    if res_id + 1 > len(fasta_sequence):
                        end = res_id
                    else:
                        end = res_id+1
                    selected_residues += [np.arange(start, end + 1)]
                    selection.append(self.residue_range(f"{start}:{chain}", f"{end}:{chain}"))
                except KeyError:
                    print(f"Missing residue {res_id} in chain {chain} not found in alignment. Skipping.")
            print(f'selected_residues: {set(np.hstack(selected_residues))}')
            return Selection(*selection)

        ### This is a special restraint for P77754_6BIE_A and P31142 to prevent odd rebuilding of the very long temrinal tails. 
        # it is not required for other proteins and is commented out for the rest
        def special_restraints(self, aln):
            rsr = self.restraints
            
            # Add excluded volume restraints to prevent clashes
            rsr.make(self, restraint_type='stereo', spline_on_site=True)

            # Add custom distance restraints to keep the rebuilt tail away from the core
            #at1 = self.atoms['CA:153:A']  # Example: alpha carbon of the first residue (adjust accordingly)

            # Selection for P77754_6BIE_A
            #at1_range = self.residue_range('150:A', '161:A') # terminal tail to be rebuilt that i do not want interacting with protein
            #at2_range = self.residue_range('55:A', '140:A')

            # Selection for P31142_1URH_A
            #at1_range = self.residue_range('272:A', '281:A') # terminal tail to be rebuilt that i do not want interacting with protein
            #at2_range = self.residue_range('5:A', '260:A')
            
            # Loop through all atoms in both ranges and apply lower-bound distance restraints
            for at1_res in at1_range:
                for at2_res in at2_range:
                    # Apply the distance restraint between atoms in at1_res and at2_res
                    for at1_atom in at1_res.atoms:
                        for at2_atom in at2_res.atoms:
                            rsr.add(forms.LowerBound(group=physical.xy_distance, feature=features.Distance(at1_atom, at2_atom), mean=20.0, stdev=0.1))

    # Set up the model and perform the rebuilding
    mdl = MyModel(env, alnfile=alignment_file,
                  knowns='existing', sequence=target_code)
    mdl.starting_model = 1
    mdl.ending_model = 1
    mdl.make()
    
    # Save the rebuilt model
    mdl.write(file=output_rebuilt_pdb_file)
    logging.info(f"Rebuilt PDB saved as {output_rebuilt_pdb_file}")

    ## correct the PDB chain
    change_chain_id_biopython(output_rebuilt_pdb_file, output_rebuilt_pdb_file, new_chain=missing_residues_df['chain_id'].values[0])

def change_chain_id_biopython(input_pdb_file, output_pdb_file, new_chain='B'):
    """
    Read a PDB file, change the chain ID to the specified value, and save the modified structure.
    
    Args:
    - input_pdb_file: Path to the input PDB file.
    - output_pdb_file: Path to save the modified PDB file.
    - new_chain: The new chain ID to set (default: 'B').
    """
    # Create a PDB parser
    parser = PDBParser(QUIET=True)
    
    # Parse the structure
    structure = parser.get_structure('structure', input_pdb_file)
    
    # Loop over all chains in all models and residues, and change the chain ID
    for model in structure:
        for chain in model:
            chain.id = new_chain  # Change the chain ID to 'B'

    # Save the modified structure using PDBIO
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb_file)

    logging.info(f"Saved modified PDB with chain ID '{new_chain}' to {output_pdb_file}")   

def find_continuous_ranges(arr):
    # Ensure the array is sorted and a numpy array
    arr = np.sort(arr)
    
    # Initialize variables to store ranges
    ranges = []
    
    # Start the first range
    start = arr[0]
    prev = arr[0]
    
    # Iterate over the array to find continuous ranges
    for i in arr[1:]:
        if i != prev + 1:
            # End the current range and start a new one
            ranges.append(f"{start}:{prev}")
            start = i
        prev = i
    
    # Append the last range
    ranges.append(f"{start}:{prev}")
    
    # Return ranges as a space-separated string
    return " ".join(ranges)
############## MAIN #################
def main():
    
    script_name = f'Get_PDBs'
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--candidates", type=str, required=True, help="Path to candidates file")
    parser.add_argument("--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("--log", type=str, required=True, help="Path to logging file")
    parser.add_argument("--cath", type=str, required=True, help="Path to latest cath-domain-boundaries-seqreschopping.txt file")
    args = parser.parse_args()

    ## make output folder
    if not os.path.exists(args.outpath):
        os.makedirs(args.outpath)
        print(f'Made directory: {args.outpath}')

    # Setup logging configuration
    logging.basicConfig(filename=args.log, level=logging.INFO, format='%(asctime)s %(message)s')
    logging.info(f'{"#"*100}\nNEW RUN {script_name}')

    # initialize the simulation object 
    grab = Grabber(args)

    # Start the simulation
    grab.run()


if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
print(f'NORMAL TERMINATION: {end_time - start_time}')
logging.info(f'NORMAL TERMINATION: {end_time - start_time}')