#!/usr/bin/env python3
import logging, os, sys
import time
import argparse
import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.colors as mcolors
import seaborn as sns
from scipy.stats import mode, permutation_test, ttest_1samp, ttest_ind
import pickle
#pd.set_option('display.max_rows', 5000)

class Analysis:
    """
    A class to collect the G, Q, K, and clustered EntInfo data into a single binary dictionary that can then be read into other analysis
    ALong with this I remove annotated mirror image trajectories from the final set of data and calcualte the Q, G, K stats for the metastable states in the MSM model
    We also identify native frames either by MSM or Ref
    """
    #######################################################################################
    def __init__(self, args):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
        ("--candidates", type=str, required=True, help="A file containing two columns. The candidate tag and the groupID")
        ("--toplevel", type=str, required=True, help="file containing relative paths to either native state GQ files or the MSM file for various proteins")
        ("--outname", type=str, required=True, help="base name for output files")
        ("--MirrorOutfile", type=str, required=True, help="file containing trajectories identified as a mirror")
        """

        # parse the parameters 
        self.candidates = pd.read_csv(args.candidates)
        logging.info(f'candidates:\n{self.candidates}')

        self.outname = args.outname
        logging.info(f'outname: {self.outname}')

        self.toplevel = args.toplevel
        logging.info(f'toplevel: {self.toplevel}')
        print(f'toplevel: {self.toplevel}')

        self.MirrorOutfile = args.MirrorOutfile
        print(f'MirrorOutfile: {self.MirrorOutfile}')

        ## Get the set of clustered entanglement files
        self.ClusterEntInfoFiles = glob.glob(os.path.join(self.toplevel, '*/Cluster_ChangesInEnt/*_clustered.EntInfo'))
        print(f'Number of clustered EntInfo files found: {len(self.ClusterEntInfoFiles)}')


        ## get the set of G, Q, and K files
        self.QFiles = glob.glob(os.path.join(self.toplevel, '*/Q/*.Q'))
        print(f'Number of Q files found: {len(self.QFiles)}')

        self.GFiles = glob.glob(os.path.join(self.toplevel, '*/G/*.G'))
        print(f'Number of G files found: {len(self.GFiles)}')

        self.EntInfoFiles = glob.glob(os.path.join(self.toplevel, '*/G/*.EntInfo'))
        print(f'Number of EntInfo files found: {len(self.EntInfoFiles)}')

        self.KFiles = glob.glob(os.path.join(self.toplevel, '*/Mirror/Quench/*.dat'))
        print(f'Number of K files found: {len(self.KFiles)}')


        ## get the set of native G, Q, and K files
        self.nQFiles = glob.glob(os.path.join(self.toplevel, '*/Native/Q/*.Q'))
        print(f'Number of Native Q files found: {len(self.QFiles)}')

        self.nGFiles = glob.glob(os.path.join(self.toplevel, '*/Native/G/*.G'))
        print(f'Number of Native G files found: {len(self.GFiles)}')

        self.nKFiles = glob.glob(os.path.join(self.toplevel, '*/Mirror/Native/*.dat'))
        print(f'Number of Native K files found: {len(self.KFiles)}')


        ## get the MSM files
        self.MSMFiles = glob.glob(os.path.join(self.toplevel, '*/BuildKineticModel/*_MSMmapping.csv'))
        print(f'Number of MSM files found: {len(self.KFiles)}')

    #######################################################################################

    #######################################################################################
    def load_OP(self,):
        """
        Loops through the file paths in the GQfiles file and determines if it is either 
        1. GQ data from the native sims
        2. MSM data file 
        3. the tag is in the caidate list
        and loads the data as appropriate
        """

        data = {}
        print(self.candidates)
        for gene, pdb, chain in self.candidates[['gene', 'pdb', 'chain']].values:
            tag = f'{gene}_{pdb}_{chain}'
            print(f'Loading files for {gene} {pdb} {chain} {tag}')

            if tag in data:
                continue

            if tag not in data:
                data[tag] = {'df':pd.DataFrame(), 'K':{}}


            ###########################################################################
            ## Load Quench Q, G, K data
            KFiles = [f for f in self.KFiles if tag in f]
            print(f'KFiles: {len(KFiles)}')
            for f in KFiles:
                traj = f.split('_')[-2].replace('t', '')
                #print(f, traj)
                data[tag]['K'][traj] = np.loadtxt(f, usecols=(-1), dtype=str)[1:].astype(float)

            ###########################################################################
            ## Load MSM files
            MSMFiles = [f for f in self.MSMFiles if tag in f]
            if len(MSMFiles) > 1:
                raise ValueError(f"There should only be 1 or 0 MSM file but there are {len(MSMFiles)}")
                
            if MSMFiles:
                df = pd.read_csv(MSMFiles[0])
                df['K'] = np.nan

                ## add the K data to a column for easy manipulation later 
                new_dfs = []
                for traj, traj_df in df.groupby('traj'):
                    traj_K = data[tag]['K'][str(traj)]
                    if len(traj_df) == len(traj_K):
                        temp_df = traj_df.copy()
                        temp_df['K'] = traj_K
                        new_dfs += [temp_df]
                        #df.loc[traj_df.index,'K'] = traj_K

                    elif len(traj_K) > len(traj_df):
                        #print(traj_K, traj_K.shape)
                        temp_df = traj_df.copy()
                        temp_df['K'] = traj_K[-len(traj_df):]
                        new_dfs += [temp_df]
                        #df.loc[traj_df.index,'K'] = traj_K[-len(traj_df):]

                    else:
                        print(traj_df)
                        print(traj_K, traj_K.shape)
                        raise ValueError(f'The length of the MSM traj {len(traj_df)} != the K values {len(traj_K)}')

                df = pd.concat(new_dfs)
                ## remove any traj that has been identified to have a mirror
                data[tag]['df'] = df
                print(df)


        self.data = data
        logging.info(f'Number of genes in data: {len(self.data)}')    
        print(f'Number of genes in data: {len(self.data)}')            
    #######################################################################################  

    #######################################################################################
    def Find_potential_Mirrors(self,):
        print(f'Finding potential mirrors')
        outdf = {'gene':[], 'pdb':[], 'chain':[], 'traj':[], '<Q>':[], '<G>':[], '<K>':[], 'Mirror':[]}
        for tag in self.data:
            gene, pdb, chain = tag.split('_')
            print(f'{tag} {gene} {pdb} {chain}')

            df = self.data[tag]['df']
            for traj, traj_df in df.groupby('traj'):
                last_10_percent = traj_df.tail(int(len(traj_df) * 0.1))
                #print(last_10_percent)
                meanQ = np.mean(last_10_percent['Q'].values)    
                meanG = np.mean(last_10_percent['G'].values)
                meanK = np.mean(last_10_percent['K'].values)

                Mirror = 'False'
                if meanK < 0.3:
                    Mirror = 'True'
                
                if meanG == 0 and meanK < 0.5:
                    Mirror = 'True'
                
                if meanK < 0.5:
                    outdf['gene'] += [gene]
                    outdf['pdb'] += [pdb]
                    outdf['chain'] += [chain]
                    outdf['traj'] += [traj]
                    outdf['<Q>'] += [meanQ]
                    outdf['<G>'] += [meanG]
                    outdf['<K>'] += [meanK]
                    outdf['Mirror'] += [Mirror]

        outdf = pd.DataFrame(outdf)
        print(f'outdf:\n{outdf}')
        outdf.to_csv(self.MirrorOutfile)
        print(f'SAVED: {self.MirrorOutfile}')
    #######################################################################################



#######################################################################################
def bootstrap(data):
    boot_means = []
    for b in range(10):
        boot_samp = np.random.choice(data, size=len(data))
        boot_mean = np.mean(boot_samp)
        #print(b, boot_mean)
        boot_means += [boot_mean]

    lb = np.percentile(boot_means, 2.5)
    ub = np.percentile(boot_means, 97.5)
    return (lb, ub)
#######################################################################################

#######################################################################################
def get_stats(arr):
    """
    Get the <> and 95% ci for a stats array
    """
    mean = np.mean(arr)
    std = np.std(arr)
    lb, ub = bootstrap(arr)
    median = np.median(arr)
    return (median, mean, std, lb, ub)
#######################################################################################

#######################################################################################
def statistic(x, y, axis):
    return np.mean(x, axis=axis) - np.mean(y, axis=axis)
#######################################################################################

############## MAIN #################
def main():
    
    script_name = f'FindMirrorStructures'
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--candidates", type=str, required=True, help="A file containing two columns. The candidate tag and the groupID")
    parser.add_argument("--toplevel", type=str, required=True, help="file containing relative paths to either native state GQ files or the MSM file for various proteins")
    parser.add_argument("--MirrorOutfile", type=str, required=True, help="file containing trajectories identified as a mirror")
    parser.add_argument("--outname", type=str, required=True, help="base name for output files")
    args = parser.parse_args() 
   
    # Setup logging configuration
    logfile = os.path.join(logs, f'{args.outname}.log')
    print(f'logfile: {logfile}')
    logging.basicConfig(filename=logfile, level=logging.INFO, format='%(asctime)s %(message)s')
    logging.info(f'{"#"*50}NEW RUN {script_name}{"#"*50}')

    
    # Step 0: initialize the simulation object 
    anal = Analysis(args)

    # Define a file to save the intermeditate dictionary too for quicker loading
    DATA_outfile = os.path.join(data_path, f'DATA.pkl')
    logging.info(f'DATA_outfile: {DATA_outfile}')

    if not os.path.exists(DATA_outfile):
        logging.info(f'No DATA.pkl file found. Generating from scratch.')
        # Step 1: load the G and Q data
        anal.load_OP()

        # Step 2: for each candidate get the native G and Q states
        anal.Find_potential_Mirrors()

    print(f'logfile: {logfile}')

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()

print(f'NORMAL TERMINATION: {time.time() - start_time}')
logging.info(f'NORMAL TERMINATION: {time.time() - start_time}')