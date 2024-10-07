#!/usr/bin/env python3
import requests, logging, os, sys
from Bio.PDB import PDBParser, PDBIO, Select, Polypeptide
from Bio.PDB.PDBList import PDBList
from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
import time
import argparse
import pandas as pd
from modeller import *
from modeller.automodel import *

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
        """

        # parse the parameters 
        self.outpath = args.outpath
        self.candidates = pd.read_csv(args.candidates)
        self.log = args.log

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
        
        self.Modellerpath = os.path.join(self.outpath, 'MODLLER')
        if not os.path.exists(self.Modellerpath):
            os.makedirs(self.Modellerpath)
            print(f'Made directory: {self.Modellerpath}')

    ### Gather the data from the interwebs
    def run(self,):
        for setID, set_df in self.candidates.groupby('set'):
            print(set_df)
            for gene, pdb, chain in set_df[['gene', 'pdb', 'chain']].values:
                logging.info(f'{"#"*20} {gene} {pdb} {chain} {"#"*20}')
                output_pdb_file = os.path.join(self.PDBspath, f'{gene}_{pdb}_{chain}.pdb')
                output_fasta_file = os.path.join(self.FASTApath, f'{gene}_{pdb}_{chain}.fasta')
                output_missing_file = os.path.join(self.MISSINGpath, f'{gene}_{pdb}_{chain}.missing')
                output_rebuilt_pdb_file = os.path.join(self.PDBspath, f'{gene}_{pdb}_{chain}_rebuilt.pdb')
                working_dir = os.path.join(self.Modellerpath, f'{gene}_{pdb}_{chain}/')
                if not os.path.exists(working_dir):
                    os.makedirs(working_dir)
                    print(f'Made directory: {working_dir}')

                fasta_sequence, missing_residues_df = self.process_pdb_and_sequence(gene, pdb, chain, output_pdb_file=output_pdb_file, output_fasta_file=output_fasta_file, output_missing_file=output_missing_file)
                if len(missing_residues_df) != 0:
                    rebuild_missing_residues(output_pdb_file, fasta_sequence, missing_residues_df, output_rebuilt_pdb_file, working_dir)
                else:
                    logging.info(f'No missing residues found. No need to rebuild.')
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
    def process_pdb_and_sequence(self, uniprot_id, pdb_id, chain_id, output_pdb_file='extracted_chain.pdb', output_fasta_file='extracted_chain.fasta', output_missing_file='extracted_chain.missing'):
        # Step 1: Download the PDB file
        pdbl = PDBList()
        pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir=self.PDBspath, file_format='pdb')
        
        # remove obsolete directory made by retreival
        os.rmdir('obsolete')
       
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
        
        # Save only the specified chain
        class ChainSelect(Select):
            def accept_chain(self, c):
                return c.id == chain_id
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb_file, select=ChainSelect())
        logging.info(f"Extracted chain {chain_id} saved to {output_pdb_file}")

        # Step 3: Extract the amino acid sequence of the specified chain
        ppb = Polypeptide.PPBuilder()
        polypeptides = ppb.build_peptides(chain)
        PDB_seq = []
        for pp in polypeptides:
            PDB_seq += [str(pp.get_sequence())]
        #PDB_seq = polypeptides[0].get_sequence()
        PDB_seq = ''.join(PDB_seq)
        logging.info(polypeptides)
        logging.info(f'{uniprot_id} {pdb_id} {chain_id } PDB_seq: {PDB_seq} {len(PDB_seq)}')


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
        else:
            logging.info(f"Failed to download FASTA file for UniProt ID {uniprot_id}. Status code: {response.status_code}")

        # Step 5: align the two sequences and find the missing residues
        missing_residues = self.find_missing_residues(fasta_record.seq, PDB_seq)
        missing_residues['uniprot_id'] = uniprot_id
        missing_residues['pdb_id'] = pdb_id
        missing_residues['chain_id'] = chain_id
        logging.info(f'missing_residues:\n{missing_residues}')
        missing_residues.to_csv(output_missing_file, index=False)
        logging.info(f'SAVED: {output_missing_file}')

        return fasta_record.seq, missing_residues

    ### Align two PDB sequenes and find missing residues
    def find_missing_residues(self, canonical_seq, pdb_seq):
        # Align the two sequences
        alignments = pairwise2.align.globalxx(canonical_seq, pdb_seq)
        
        # Use the first (best) alignment
        alignment = alignments[0]
        fasta_aligned, pdb_aligned, score, begin, end = alignment
        
        # Print the alignment (optional)
        alin_str = format_alignment(*alignment)
        logging.info(f'{"#"*20} Alignment {"#"*20}\n{alin_str}')
        
        # Identify missing residues in the PDB sequence relative to the canonical sequence
        missing_residues = {'AA':[], 'resid':[]}
        for i, (f_res, p_res) in enumerate(zip(fasta_aligned, pdb_aligned)):
            # Check for gaps in the PDB sequence where the FASTA has a residue
            if f_res != '-' and p_res == '-':
                #missing_residues.append((i + 1, f_res))  # Use 1-based indexing
                missing_residues['AA'] += [f_res]
                missing_residues['resid'] += [i + 1]
        missing_residues = pd.DataFrame(missing_residues)

        return missing_residues

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
    logging.info(fasta_sequence, type(fasta_sequence))
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
            for _, row in missing_residues_df.iterrows():
                res_id = int(row['resid'])
                chain = 'A'
                # Add the missing residue position as a selection target
                try:
                    selection.append(self.residue_range(f"{res_id}:{chain}", f"{res_id}:{chain}"))
                except KeyError:
                    print(f"Missing residue {res_id} in chain {chain} not found in alignment. Skipping.")
            return Selection(*selection)
        

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

############## MAIN #################
def main():
    
    script_name = f'Get_PDBs'
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--candidates", type=str, required=True, help="Path to candidates file")
    parser.add_argument("--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("--log", type=str, required=True, help="Path to logging file")
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