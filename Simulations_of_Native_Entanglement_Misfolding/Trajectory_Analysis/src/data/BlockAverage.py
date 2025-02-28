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
from rpy2.robjects import r, FloatVector
from rpy2.robjects.packages import importr
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
        ("--outpath", type=str, required=True, help="Path to output directory")
        ("--candidates", type=str, required=True, help="A file containing two columns. The candidate tag and the groupID")
        ("--toplevel", type=str, required=True, help="file containing relative paths to either native state GQ files or the MSM file for various proteins")
        ("--outname", type=str, required=True, help="base name for output files")
        ("--Mirrorfile", type=str, required=True, help="file containing trajectories identified as a mirror")
        ("--setID", type=int, required=True, help="setID to use when collecting OP files")
        """

        # parse the parameters 
        self.candidates = pd.read_csv(args.candidates)
        logging.info(f'candidates:\n{self.candidates}')

        self.outpath = args.outpath
        logging.info(f'outpath: {self.outpath}')

        self.outname = args.outname
        logging.info(f'outname: {self.outname}')

        self.toplevel = args.toplevel
        logging.info(f'toplevel: {self.toplevel}')
        print(f'toplevel: {self.toplevel}')

        self.Mirrorfile = args.Mirrorfile
        print(f'Mirrorfile: {self.Mirrorfile}')
        self.Mirror_df = pd.read_csv(self.Mirrorfile)
        self.Mirror_df = self.Mirror_df[self.Mirror_df['Mirror'] == True]
        print(self.Mirror_df)

        self.setID = args.setID
        print(f'setID: {self.setID} {type(self.setID)}')

        ## Get the set of clustered entanglement files
        self.ClusterEntInfoFiles = glob.glob(os.path.join(self.toplevel, '*/Cluster_ChangesInEnt/*_clustered.EntInfo'))
        print(f'Number of clustered EntInfo files found: {len(self.ClusterEntInfoFiles)}')


        ## get the set of G, Q, and K files
        self.QFiles = glob.glob(os.path.join(self.toplevel, '*/Q/*.Q'))
        print(f'Number of Q files found: {len(self.QFiles)}')

        self.GFiles = glob.glob(os.path.join(self.toplevel, '*/Cluster_ChangesInEnt/*_clustered.G'))
        print(f'Number of G files found: {len(self.GFiles)}')

        self.KFiles = glob.glob(os.path.join(self.toplevel, '*/Mirror/Quench/*.dat'))
        print(f'Number of K files found: {len(self.KFiles)}')


        ## get the set of native G, Q, and K files
        self.nQFiles = glob.glob(os.path.join(self.toplevel, '*/Native/Q/*.Q'))
        print(f'Number of Native Q files found: {len(self.QFiles)}')

        self.nGFiles = glob.glob(os.path.join(self.toplevel, '*/Native/Cluster_ChangesInEnt/*_clustered.G'))
        print(f'Number of Native G files found: {len(self.nGFiles)}')

        self.nKFiles = glob.glob(os.path.join(self.toplevel, '*/Mirror/Native/*.dat'))
        print(f'Number of Native K files found: {len(self.KFiles)}')


        ## get the MSM files
        #print(os.path.join(self.toplevel, f'*/BuildKineticModel/setID{self.setID}/*_MSMmapping.csv'))
        self.MSMFiles = glob.glob(os.path.join(self.toplevel, f'*/BuildKineticModel/setID{self.setID}/*_MSMmapping.csv'))
        print(f'Number of MSM files found: {len(self.MSMFiles)}')

        ## make logging dir
        self.data_path = os.path.join(self.outpath, 'DATA')
        if not os.path.exists(self.data_path):
            os.makedirs(self.data_path)
            print(f'Made directory: {self.data_path}')   
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

        Quench_dfs = []
        Native_dfs = []
        self.candidates = self.candidates[self.candidates['set'] == self.setID] ## only get those candidates in this set
        print(self.candidates)
        for gene, pdb, chain in self.candidates[['gene', 'pdb', 'chain']].values:
            tag = f'{gene}_{pdb}_{chain}'
            print(f'Loading files for {gene} {pdb} {chain} {tag}')

            ###########################################################################
            ### find the MSM file and load it
            MSMFiles = [f for f in self.MSMFiles if tag in f]
            #print(f'MSMFiles: {MSMFiles}')
            if len(MSMFiles) > 1:
                raise ValueError(f"There should only be 1 or 0 MSM file but there are {len(MSMFiles)}")
            df = pd.read_csv(MSMFiles[0])
            df.loc[:, 'gene'] = gene
            df.loc[:, 'pdb'] = pdb
            df.loc[:, 'chain'] = chain            
            df.loc[:, 'Mirror'] = False
            #print(df)
            ###########################################################################


            ###########################################################################
            ## Load Quench K data and add it to the df
            KFiles = [f for f in self.KFiles if tag in f]
            #print(f'KFiles: {len(KFiles)}')
            df['K'] = np.nan #initialize the column with nan values
            for f in KFiles:
                traj = f.split('_')[-2].replace('t', '')
                #print(f, traj)
                K = np.loadtxt(f, usecols=(-1), dtype=str)[1:].astype(float)
                frames = df.loc[(df['traj'] == int(traj)),'frame'].values
                K = K[frames] ## only get those frames of K that were found in the MSM datafile. This depends on the setID 
                df.loc[(df['traj'] == int(traj)), 'K'] = K
            #print(df)
            ###########################################################################


            ###########################################################################
            ### Get the mirror trajs and mark them
            Mirrors = self.Mirror_df[(self.Mirror_df['gene'] == gene) & (self.Mirror_df['pdb'] == pdb) & (self.Mirror_df['chain'] == chain)]
            Mirror_trajs = Mirrors['traj'].values
            #print(f'Mirror_trajs: {Mirror_trajs}')
            if len(Mirror_trajs) > 0:
                for traj in Mirror_trajs:
                    df.loc[(df['traj'] == traj), 'Mirror'] = True
            #print(df, df['Mirror'].sum())

            Quench_dfs += [df]
            ###########################################################################


            ###########################################################################
            ## Load Native Q, G, K data
            ndf = {'traj':[], 'frame':[], 'Q':[], 'G':[], 'K':[]}
            QFiles = [f for f in self.nQFiles if tag in f]
            #print(f'nQFiles: {len(QFiles)}')
            for f in QFiles:
                traj = f.split('_')[-1].replace('t', '').replace('.Q', '')
                #print(f, traj)
                Q = pd.read_csv(f)['Q'].values

                ## get corresponding Native G file
                Gf = [f for f in self.nGFiles if f'{tag}_t{traj}_' in f][0]
                #print(Gf)
                G = pd.read_csv(Gf)['G'].values

                ## get corresponding Native K file
                Kf = [f for f in self.nKFiles if f'{tag}' in f]
                Kf = [f for f in Kf if f'K_{traj}.dat' in f][0]
                #print(Kf)
                K = np.loadtxt(Kf, usecols=(-1), dtype=str)[1:].astype(float)           

                if Q.shape != G.shape:
                    raise ValueError(f'The shape of the native Q data {Q.shape} != the G {G.shape}')
            
                for frame, (Q, G, K) in enumerate(zip(Q, G, K)):
                    ndf['traj'] += [traj]
                    ndf['frame'] += [frame]
                    ndf['Q'] += [Q]
                    ndf['G'] += [G]
                    ndf['K'] += [K]

            ndf = pd.DataFrame(ndf)
            ndf.loc[:, 'gene'] = gene
            ndf.loc[:, 'pdb'] = pdb
            ndf.loc[:, 'chain'] = chain
            #print(f'Native GQ df:\n{ndf}')
            Native_dfs += [ndf]
            ###########################################################################


        ###########################################################################
        ### Combine the Quench_dfs and the Native_dfs and save them
        Quench_dfs = pd.concat(Quench_dfs, ignore_index=True)
        #print(f'Quench_dfs:\n{Quench_dfs}')
        #Quench_outfile = os.path.join(self.data_path, f'Quench_Collected_GQK.csv')
        #Quench_dfs.to_csv(Quench_outfile, index=False)
        #print(f'SAVED: {Quench_outfile}')

        Native_dfs = pd.concat(Native_dfs, ignore_index=True)
        #print(f'Native_dfs:\n{Native_dfs}')
        Native_outfile = os.path.join(self.data_path, f'Native_Collected_GQK.csv')
        Native_dfs.to_csv(Native_outfile, index=False)
        #print(f'SAVED: {Native_outfile}')

        return Quench_dfs, Native_dfs
    #######################################################################################  

    #######################################################################################
    def BlockAverage(self, df, label):
        print(f'Block Averaging')
        """
        1. for each block size split the data into consecutive chunks M of size n
        2. for each chunk calculate the mean <X(m)> where m is a chunk of the set M
        3. calculate the overall mean of the traj with the blocks <X>
        4. calculate the standard deviation of the block means sigma and its associated error
        """
        outfile = os.path.join(self.data_path, f'{label}_BlockAverage.csv')
        if not os.path.exists(outfile):
            block_res = {'gene':[], 'traj':[], 'blocking_opt':[], 'num_blocks':[], 'overall_mean':[], 'overall_std':[], 'var_est':[], 'overall_std_err':[]}
            for gene, gene_df in df.groupby('gene'):
                for traj, traj_df in gene_df.groupby('traj'):
                    #print(traj_df)
                    print(label, gene, traj)
                    Qdata = traj_df['Q'].values
                    num_blocks = len(Qdata)
                    blocked_array = Qdata.copy()
                    blocking_opt = 1
                    while num_blocks >= 2:
                        blocked_array = blocking_transformation(blocked_array)
                        num_blocks = len(blocked_array)

                        overall_mean = np.mean(blocked_array)
                        overall_std = np.std(blocked_array)
                        var_est = (overall_std**2)/(num_blocks - 1)
                        overall_std_err = np.sqrt((2*overall_std**4)/(num_blocks-1)**3)

                        block_res['gene'] += [gene]
                        block_res['traj'] += [traj]
                        block_res['blocking_opt'] += [blocking_opt]
                        block_res['num_blocks'] += [num_blocks]
                        block_res['overall_mean'] += [overall_mean]
                        block_res['overall_std'] += [overall_std]
                        block_res['var_est'] += [var_est]
                        block_res['overall_std_err'] += [overall_std_err]
                        blocking_opt += 1

            block_res = pd.DataFrame(block_res)
            print(f'block_res:\n{block_res}')
            outfile = os.path.join(self.data_path, f'{label}_BlockAverage.csv')
            block_res.to_csv(outfile, index=False)
            print(f'SAVED: {outfile}')
        else:
            block_res = pd.read_csv(outfile)
            print(f'LOADED: {outfile}') 
            print(f'block_res:\n{block_res}')

        ##################################################################
        ## plot data
        for gene, gene_df in block_res.groupby('gene'):
            print(gene_df)

            # Initialize the plot
            plt.figure(figsize=(10, 7))

            # Plot 'overall_std' vs 'block_size' for each 'traj'
            for traj, traj_df in gene_df.groupby('traj'):
                #print(traj_df)
                # Plot
                plt.errorbar(
                    traj_df['blocking_opt'], 
                    traj_df['var_est'], 
                    yerr=traj_df['overall_std_err'], 
                    fmt='o-', 
                    label=f'Traj {traj}',
                    capsize=5
                )

            # Plot customization
            plt.title('Overall Standard Deviation vs Block Size')
            plt.xlabel('blocking_opt')
            plt.ylabel('Overall Standard Deviation')
            plt.ylim(bottom=0.00000001)
            #plt.yscale('log')
            plt.grid(True)
            #plt.legend()

            # Show the plot
            plt.tight_layout()
            outfile = os.path.join(self.data_path, f'{label}_BlockAverage_{gene}.png')
            plt.savefig(outfile)
            print(f'SAVED: {outfile}')
            plt.close()
 

#######################################################################################
def bootstrap(data, reps=10000):
    boot_means = []
    for b in range(reps):
        boot_samp = np.random.choice(data, size=len(data))
        boot_mean = np.mean(boot_samp)
        #print(b, boot_mean)
        boot_means += [boot_mean]

    lb = np.percentile(boot_means, 2.5)
    ub = np.percentile(boot_means, 97.5)
    return (lb, ub)
#######################################################################################

#######################################################################################
def get_stats(arr, reps=10000):
    """
    Get the <> and 95% ci for a stats array
    """
    means = []
    modes = []
    window_size = 200 # ~15ns
    #print(arr)
    for a in arr:
        means += [np.mean(a[i:i + window_size]) for i in range(len(a) - window_size + 1)]
        modes += [mode(a[i:i + window_size]).mode for i in range(len(a) - window_size + 1)]
    #print(means)
    #print(modes)
    means = np.hstack(means)
    modes = np.hstack(modes)
    #print(means, means.shape)
    #print(modes, modes.shape)

    mean = np.mean(means)
    std = np.std(means)
    lb, ub = bootstrap(means, reps=reps)
   
    mean_mode = np.mean(modes)
    std_mode = np.std(modes)
    mode_lb, mode_ub = bootstrap(modes, reps=reps)

    return (mean, std, lb, ub, mean_mode, std_mode, mode_lb, mode_ub)
#######################################################################################

#######################################################################################
def blocking_transformation(arr):
    """
    Function to take an array and perform the blocking operation stated
        equations 20 and 21 from https://pubs.aip.org/aip/jcp/article/91/1/461/91554/Error-estimates-on-averages-of-correlated
    """
    blocked_array = [np.mean(arr[i:i+2]) for i in range(0, len(arr) - 2, 2)]
    return blocked_array

#######################################################################################

#######################################################################################
def statistic(x, y, axis):
    return np.mean(x, axis=axis) - np.mean(y, axis=axis)
#######################################################################################


############## MAIN #################
def main():
    
    script_name = f'CollectAndProcessOP'
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("--candidates", type=str, required=True, help="A file containing two columns. The candidate tag and the groupID")
    parser.add_argument("--toplevel", type=str, required=True, help="file containing relative paths to either native state GQ files or the MSM file for various proteins")
    parser.add_argument("--outname", type=str, required=True, help="base name for output files")
    parser.add_argument("--Mirrorfile", type=str, required=True, help="file containing trajectories identified as a mirror")
    parser.add_argument("--setID", type=int, required=True, help="setID to use when collecting OP files")

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
    
    ## make DATA dir
    data_path = os.path.join(args.outpath, 'DATA')
    if not os.path.exists(data_path):
        os.makedirs(data_path)
        print(f'Made directory: {data_path}')  
   

    # Setup logging configuration
    logfile = os.path.join(logs, f'{args.outname}.log')
    print(f'logfile: {logfile}')
    logging.basicConfig(filename=logfile, level=logging.INFO, format='%(asctime)s %(message)s')
    logging.info(f'{"#"*50}NEW RUN {script_name}{"#"*50}')

    
    # Step 0: initialize the simulation object 
    anal = Analysis(args)
    
    logging.info(f'No DATA.pkl file found. Generating from scratch.')
    # Step 1: load the G and Q data
    Quench_df, Native_df = anal.load_OP()

    # Step 2: for each candidate get the native G and Q states
    anal.BlockAverage(Native_df, 'Native')
    anal.BlockAverage(Quench_df, 'Quench')


    print(f'logfile: {logfile}')

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()

print(f'NORMAL TERMINATION: {time.time() - start_time}')
logging.info(f'NORMAL TERMINATION: {time.time() - start_time}')