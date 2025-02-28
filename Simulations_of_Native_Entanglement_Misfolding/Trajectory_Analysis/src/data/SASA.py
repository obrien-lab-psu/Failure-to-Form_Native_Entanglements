#!/usr/bin/env python3
import requests, logging, os, sys
import time
import argparse
import pandas as pd
import numpy as np
import glob
import MDAnalysis as mda
from scipy.spatial.distance import pdist, squareform
from topoly import lasso_type  # used pip
import itertools
import concurrent.futures
import mdtraj as md


pd.set_option('display.max_rows', 5000)

class Analysis:
    """
    A class to handel the analyssis of a C-alpha CG trajectory. 
    Current analysis available:
    (1) - Fraction of native contacts (Q)
    (2) - Change of entanglement (G)
    """
    #######################################################################################
    def __init__(self, args):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
        parser.add_argument("--outpath", type=str, required=True, help="Path to output directory")
        parser.add_argument("--outname", type=str, required=True, help="base name for output files")
        parser.add_argument("--psf", type=str, required=True, help="Path to CA protein structure file")
        parser.add_argument("--cor", type=str, required=True, help="Path to CA native coordinates file")
        parser.add_argument("--quench_dcds", type=str, required=True, help="Path to trajectory to analyze")
        parser.add_argument("--native_dcds", type=str, required=True, help="Path to native trajectory to analyze")
        parser.add_argument("--FASTApath", type=str, required=True, help="Path to directory containing FASTA files")
        parser.add_argument("--start", type=int, required=False, help="First frame to analyze 0 indexed", default=0)
        parser.add_argument("--end", type=int, required=False, help="Last frame to analyze 0 indexed", default=None)
        parser.add_argument("--stride", type=int, required=False, help="Frame stride", default=1)
        """

        # parse the parameters 
        self.outpath = args.outpath
        logging.info(f'outpath: {self.outpath}')

        self.outname = args.outname
        logging.info(f'outname: {self.outname}')

        self.psf = args.psf
        logging.info(f' psf: {self. psf}')

        self.cor = args.cor
        logging.info(f'cor: {self.cor}')

        self.FASTApath = args.FASTApath 
        logging.info(f'FASTApath: {self.FASTApath}')
        self.FASTAseq = ''.join([l.strip('\n') for l in open(self.FASTApath).readlines()[1:]])
        print(f'FASTAseq: {self.FASTAseq}')

        self.quench_dcds = args.quench_dcds 
        logging.info(f'quench_dcds: {self.quench_dcds}')

        self.native_dcds = args.native_dcds 
        logging.info(f'native_dcds: {self.native_dcds}')

        #self.ref_universe = mda.Universe(self.psf, self.cor, format='CRD')
        #print(f'ref_universe: {self.ref_universe}')

        #self.traj_universe = mda.Universe(self.psf, self.dcd, format='DCD')
        #print(f'traj_universe: {self.traj_universe}')

        self.start = args.start
        self.end = args.end
        self.stride = args.stride
        print(f'START: {self.start} | END: {self.end} | STRIDE: {self.stride}')

    #######################################################################################

    #######################################################################################
    def Zeta(self, ):
        print(f'Calculating the relative difference in SASA')
        logging.info(f'Calculating the relative difference in SASA')
        """
        calulate the average hydrophobic SASA in the reference simulations = <A(t)>ref
        for each frame then calculate the relative difference = (A(t)/<A(t)>ref)*100
        """
        # make directory for Q data if it doesnt exist
        self.ZetaPath = os.path.join(self.outpath, 'Zeta')
        if not os.path.exists(self.ZetaPath):
            os.makedirs(self.ZetaPath)
            print(f'Made directory: {self.ZetaPath}')

        ## get the native state dcd files
        native_SASA_outfile = os.path.join(self.ZetaPath, f'{self.outname}_Native.SASA')
        if not os.path.exists(native_SASA_outfile):
            native_SASA_df = []
            native_dcd_files = glob.glob(self.native_dcds)
            for f in native_dcd_files:
                traj = f.split('/')[-1].split('.')[0]
                logging.info(f, traj)
                df = calculate_sasa_timeseries(f, self.psf, self.FASTAseq)
                df['traj'] = traj
                #print(df)
                native_SASA_df += [df]
            native_SASA_df = pd.concat(native_SASA_df)
            native_SASA_df.to_csv(native_SASA_outfile, index=False)
            logging.info(f'SAVED: {native_SASA_outfile}')
        else:
            native_SASA_df = pd.read_csv(native_SASA_outfile)
            logging.info(f'LOADED: {native_SASA_outfile}')
        print(f'native_SASA_df:\n{native_SASA_df}')


        ## get the quench state dcd files
        quench_SASA_outfile = os.path.join(self.ZetaPath, f'{self.outname}_quench.SASA')
        if not os.path.exists(quench_SASA_outfile):
            quench_SASA_df = []
            quench_dcd_files = glob.glob(self.quench_dcds)
            for f in quench_dcd_files:
                traj = f.split('/')[-1].split('_')[-2].replace('t','')
                logging.info(f, traj)
                df = calculate_sasa_timeseries(f, self.psf, self.FASTAseq)
                df['traj'] = traj
                #print(df)
                quench_SASA_df += [df]
            quench_SASA_df = pd.concat(quench_SASA_df)
            quench_SASA_df.to_csv(quench_SASA_outfile, index=False)
            logging.info(f'SAVED: {quench_SASA_outfile}')
        else:
            quench_SASA_df = pd.read_csv(quench_SASA_outfile)
            logging.info(f'LOADED: {quench_SASA_outfile}')
        print(f'quench_SASA_df:\n{quench_SASA_df}')

        ## Calculate zeta for the quench data
        quench_Zeta_outfile = os.path.join(self.ZetaPath, f'{self.outname}_quench.Zeta')
        meanAnative = np.mean(native_SASA_df['SASA'].values)
        print(f'meanAnative: {meanAnative}')

        quench_SASA_df['Zeta'] = ((quench_SASA_df['SASA']/meanAnative)-1)*100
        print(f'quench_SASA_df:\n{quench_SASA_df}')
        quench_SASA_df.to_csv(quench_Zeta_outfile, index=False)
        logging.info(f'SAVED: quench_Zeta_outfile: {quench_Zeta_outfile}')
        
    #######################################################################################  

    #######################################################################################

#######################################################################################
# Function to calculate SASA timeseries
def calculate_sasa_timeseries(dcd_file, pdb_file, FASTA, amino_acids=['I', 'V', 'L', 'F', 'C', 'M', 'A', 'G', 'T', 'S', 'W', 'Y', 'P']):
    """
    Calculate the total SASA timeseries for a protein trajectory.

    Parameters:
    - dcd_file: str, path to the .dcd trajectory file.
    - pdb_file: str, path to the .pdb topology file.

    Returns:
    - DataFrame with columns ['frame', 'SASA'].
    """
    # Load the trajectory
    traj = md.load(dcd_file, top=pdb_file)

    selected_indices = [i for i,A in enumerate(FASTA) if A in amino_acids]
    #print(f'# selected_indices: {len(selected_indices)}')

    # Calculate SASA for each frame
    sasa_per_frame = md.shrake_rupley(traj, mode='residue', probe_radius=4.0)
    sasa_per_frame = sasa_per_frame[:, selected_indices]

    # Sum SASA per frame to get total SASA
    total_sasa = sasa_per_frame.sum(axis=1)

    # Create DataFrame
    df = pd.DataFrame({'frame': np.arange(len(total_sasa)), 'SASA': total_sasa})

    return df
#######################################################################################

############## MAIN #################
def main():
    
    script_name = f'SASA'
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("--outname", type=str, required=True, help="base name for output files")
    parser.add_argument("--psf", type=str, required=True, help="Path to CA protein structure file")
    parser.add_argument("--cor", type=str, required=True, help="Path to CA native coordinates file")
    parser.add_argument("--quench_dcds", type=str, required=True, help="Path to trajectory to analyze")
    parser.add_argument("--native_dcds", type=str, required=True, help="Path to native trajectory to analyze")
    parser.add_argument("--FASTApath", type=str, required=True, help="Path to directory containing FASTA files")
    parser.add_argument("--start", type=int, required=False, help="First frame to analyze 0 indexed", default=0)
    parser.add_argument("--end", type=int, required=False, help="Last frame to analyze 0 indexed", default=None)
    parser.add_argument("--stride", type=int, required=False, help="Frame stride", default=1)
    args = parser.parse_args()

    ## make output folder
    if not os.path.exists(args.outpath):
        os.makedirs(args.outpath)
        print(f'Made directory: {args.outpath}')

    ## make logging dir
    logs = os.path.join(args.outpath, 'logs')
    if not os.path.exists(logs):
        os.makedirs(logs)
        print(f'Made directory: {logs}')    

    # Setup logging configuration
    logfile = os.path.join(logs, f'{args.outname}.log')
    print(f'logfile: {logfile}')
    logging.basicConfig(filename=logfile, level=logging.INFO, format='%(asctime)s %(message)s')
    logging.info(f'{"#"*50}NEW RUN {script_name}{"#"*50}')

    # initialize the simulation object 
    anal = Analysis(args)

    # Get the fraction of native contacts Q
    anal.Zeta()



if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    
print(f'NORMAL TERMINATION: {end_time - start_time}')
logging.info(f'NORMAL TERMINATION: {end_time - start_time}')