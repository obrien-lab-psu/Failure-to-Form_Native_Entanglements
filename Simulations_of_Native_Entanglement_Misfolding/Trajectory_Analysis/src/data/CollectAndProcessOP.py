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
            print(f'MSMFiles: {MSMFiles}')
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
            print(f'KFiles: {len(KFiles)}')
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
            print(f'Mirror_trajs: {Mirror_trajs}')
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
            print(f'nQFiles: {len(QFiles)}')
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
            print(f'Native GQ df:\n{ndf}')
            Native_dfs += [ndf]
            ###########################################################################


        ###########################################################################
        ### Combine the Quench_dfs and the Native_dfs and save them
        Quench_dfs = pd.concat(Quench_dfs, ignore_index=True)
        print(f'Quench_dfs:\n{Quench_dfs}')
        #Quench_outfile = os.path.join(self.data_path, f'Quench_Collected_GQK.csv')
        #Quench_dfs.to_csv(Quench_outfile, index=False)
        #print(f'SAVED: {Quench_outfile}')

        Native_dfs = pd.concat(Native_dfs, ignore_index=True)
        print(f'Native_dfs:\n{Native_dfs}')
        Native_outfile = os.path.join(self.data_path, f'Native_Collected_GQK.csv')
        Native_dfs.to_csv(Native_outfile, index=False)
        print(f'SAVED: {Native_outfile}')

        return Quench_dfs, Native_dfs
    #######################################################################################  

    #######################################################################################
    def get_nativeSIM_stats(self, Native_df):
        """
        Get the <> and 95% ci for G and Q in the native state sims
        """
        #print(self.candidates)
        self.stats_df = {'gene':[], 'pdb':[], 'chain':[], 'Qmean':[], 'Qstd':[], 'Q_lb':[], 'Q_ub':[], 'Qthreshold':[],
                        'Gmean':[], 'Gstd':[], 'G_lb':[], 'G_ub':[], 'Gthreshold':[], 
                        'Kmean':[], 'Kstd':[], 'K_lb':[], 'K_ub':[], 'n':[]}
        
        for gene, pdb, chain in self.candidates[['gene', 'pdb', 'chain']].values:
            tag = f'{gene}_{pdb}_{chain}'
            print(tag)
            gene_df = Native_df[Native_df['gene'] == gene]
            #print(gene_df)
            Qdata = gene_df['Q'].values
            Gdata = gene_df['G'].values
            Kdata = gene_df['K'].values       

            Qmedian, Qmean, Qstd, Q_lb, Q_ub = get_stats(Qdata)
            Gmedian, Gmean, Gstd, G_lb, G_ub = get_stats(Gdata)
            Kmedian, Kmean, Kstd, K_lb, K_ub = get_stats(Kdata)
            print(f'<Q> = {Qmean} | std(Q) = {Qstd} | 95%ci = ({Q_lb}, {Q_ub})')
            print(f'<G> = {Gmean} | std(G) = {Gstd} | 95%ci = ({G_lb}, {G_ub})')
            print(f'<K> = {Kmean} | std(K) = {Kstd} | 95%ci = ({K_lb}, {K_ub})')

            self.stats_df['gene'] += [gene]
            self.stats_df['pdb'] += [pdb]
            self.stats_df['chain'] += [chain]
            self.stats_df['Qmean'] += [Qmean]
            self.stats_df['Qstd'] += [Qstd]
            self.stats_df['Q_lb'] += [Q_lb]
            self.stats_df['Q_ub'] += [Q_ub]
            self.stats_df['Qthreshold'] += [Qmean - 3*Qstd]
            self.stats_df['Gmean'] += [Gmean]
            self.stats_df['Gstd'] += [Gstd]
            self.stats_df['G_lb'] += [G_lb]
            self.stats_df['G_ub'] += [G_ub]
            self.stats_df['Gthreshold'] += [Gmean + 3*Gstd]
            self.stats_df['Kmean'] += [Kmean]
            self.stats_df['Kstd'] += [Kstd]
            self.stats_df['K_lb'] += [K_lb]
            self.stats_df['K_ub'] += [K_ub]
            self.stats_df['n'] += [len(Qdata)]

        self.stats_df = pd.DataFrame(self.stats_df)
        print(f'Native QGK stats:\n{self.stats_df}')
        Native_stats_outfile = os.path.join(self.data_path, f'Native_Collected_GQK_STATS.csv')
        self.stats_df.to_csv(Native_stats_outfile, index=False)
        print(f'SAVED: {Native_stats_outfile}')
    #######################################################################################

    #######################################################################################
    def find_native_state(self, df):
        """
        For each candidate find the native state amongst the MSM data defined by the highest <Q> and lowest <G>
        Tag this state in a new column called NativeByMSM in the self.data[tag]['MSM'] dataframe
        Also find native frames by using the Q >= <Q>ns - 3*sigma and G <= <G>ns + 3*sigma and tag a column named NativeByRef in the self.data[tag]['MSM'] dataframe
        """
        print(f'Find the native states by MSM and Ref')
        ###################################################################################################
        df['NativeByRef'] = False
        df['NativeByMSS'] = False
        new_dfs = []
        local_mss_dfs = []
        for gene, pdb, chain in self.candidates[['gene', 'pdb', 'chain']].values:
            tag = f'{gene}_{pdb}_{chain}'
            print(tag)

            gene_df = df[df['gene'] == gene]

            #######################################################
            ### Identify native frames by the reference simulations
            Qthreshold = self.stats_df.loc[(self.stats_df['gene'] == gene), 'Qthreshold'].values[0]
            Gthreshold = self.stats_df.loc[(self.stats_df['gene'] == gene), 'Gthreshold'].values[0]
            print(f'Qthreshold: {Qthreshold} | Gthreshold: {Gthreshold}')
            gene_df.loc[(gene_df['Q'] >= Qthreshold) & (gene_df['G'] <= Gthreshold), 'NativeByRef'] = True
            #print(gene_df)
            #######################################################


            #######################################################
            ### Identify native frames by the MSM
            local_mss_GQmeans = {'gene':[], 'pdb':[], 'chain':[], 'mssID':[], 'Qmean':[], 'Qstd':[], 'Q_lb':[], 'Q_ub':[], 
                        'Gmean':[], 'Gstd':[], 'G_lb':[], 'G_ub':[], 
                        'Kmean':[], 'Kstd':[], 'K_lb':[], 'K_ub':[], 'n':[]}
            for mssID, mss_df in df.groupby(['metastablestate']):
                mssID = mssID[0]
                #print(mssID, mss_df)
                
                Qmedian, Qmean, Qstd, Q_lb, Q_ub = get_stats(mss_df['Q'].values)
                mssID_Qstats = (Qmean, Qstd, Q_lb, Q_ub)

                Gmedian, Gmean, Gstd, G_lb, G_ub = get_stats(mss_df['G'].values)
                mssID_Gstats = (Gmean, Gstd, G_lb, G_ub)

                Kmedian, Kmean, Kstd, K_lb, K_ub = get_stats(mss_df['K'].values)
                mssID_Kstats = (Kmean, Kstd, K_lb, K_ub)

                ## update dictionary traccking the mean G and Q for this protein only
                local_mss_GQmeans['gene'] += [gene]
                local_mss_GQmeans['pdb'] += [pdb]
                local_mss_GQmeans['chain'] += [chain]
                local_mss_GQmeans['mssID'] += [mssID]
                local_mss_GQmeans['Qmean'] += [Qmean]
                local_mss_GQmeans['Qstd'] += [Qstd]
                local_mss_GQmeans['Q_lb'] += [Q_lb]
                local_mss_GQmeans['Q_ub'] += [Q_ub]
                local_mss_GQmeans['Gmean'] += [Gmean]
                local_mss_GQmeans['Gstd'] += [Gstd]
                local_mss_GQmeans['G_lb'] += [G_lb]
                local_mss_GQmeans['G_ub'] += [G_ub]
                local_mss_GQmeans['Kmean'] += [Kmean]
                local_mss_GQmeans['Kstd'] += [Kstd]
                local_mss_GQmeans['K_lb'] += [K_lb]
                local_mss_GQmeans['K_ub'] += [K_ub]
                local_mss_GQmeans['n'] += [len(mss_df['Q'].values)]

            # determine the native MSM state by the highest <Q> and lowest <G>
            logging.info(f'Determine the native MSM state by the highest <Q> and lowest <G>')
            local_mss_GQmeans = pd.DataFrame(local_mss_GQmeans)
            local_mss_GQmeans = local_mss_GQmeans.sort_values(by=['Qmean', 'Gmean'], ascending=False, ignore_index=True)
            print(local_mss_GQmeans.to_string())
            nativeMSS = local_mss_GQmeans.iloc[0]['mssID']
            logging.info(f'{tag} | nativeMSS: {nativeMSS}')
            print(f'{tag} | nativeMSS: {nativeMSS}')
            gene_df.loc[:, 'NativeByMSS'] = gene_df['metastablestate'] == nativeMSS
            #print(gene_df)
            #######################################################

            new_dfs += [gene_df]
            local_mss_dfs += [local_mss_GQmeans]


        ## save the collected QGK data that has native state frames Identified
        Quench_Collected_GQK = pd.concat(new_dfs, ignore_index=True)
        print(f'Quench_Collected_GQK:\n{Quench_Collected_GQK}')
        Quench_Collected_GQK_outfile = os.path.join(self.data_path, f'Quench_Collected_GQK.csv')
        Quench_Collected_GQK.to_csv(Quench_Collected_GQK_outfile, index=False)
        print(f'SAVED: {Quench_Collected_GQK_outfile}')


        ## save the metastable state stats
        local_mss_dfs = pd.concat(local_mss_dfs, ignore_index=True)
        print(f'local_mss_dfs:\n{local_mss_dfs}')
        MSS_stats_outfile = os.path.join(self.data_path, f'MSS_Collected_GQK_STATS.csv')
        local_mss_dfs.to_csv(MSS_stats_outfile, index=False)
        print(f'SAVED: {MSS_stats_outfile}')

        return Quench_Collected_GQK
        ###################################################################################################

    ###################################################################################################
    def get_ClusterEntInfoFiles(self, df):
        """
        Collect the change in entanglement info files and grab those frames and traj that have already been controlled for mirros
        """
        collected_EntInfo_dfs = []
        for gene, gene_df in df.groupby('gene'):
            for traj, traj_df in gene_df.groupby('traj'):
                #print(traj_df)
                frames = traj_df['frame'].values
                #print(f'frames: {frames}')

                EntInfo_f = [f for f in self.ClusterEntInfoFiles if f'{gene}' in f]
                EntInfo_f = [f for f in EntInfo_f if f'_t{traj}_clustered.EntInfo' in f]
                print(EntInfo_f)

                if len(EntInfo_f) != 1:
                    print(f'Failed to find clustered EntInfo file.')
                    continue
                
                EntInfo = pd.read_csv(EntInfo_f[0], low_memory=False)
                EntInfo = EntInfo[EntInfo['Frame'].isin(frames)]
                EntInfo['gene'] = gene
                EntInfo['traj'] = traj
                
                collected_EntInfo_dfs += [EntInfo]

        ## save the collected EntInfo data that has native state frames Identified
        collected_EntInfo = pd.concat(collected_EntInfo_dfs, ignore_index=True)
        print(f'collected_EntInfo:\n{collected_EntInfo}')
        collected_EntInfo_outfile = os.path.join(self.data_path, f'Collected_EntInfo.csv')
        collected_EntInfo.to_csv(collected_EntInfo_outfile, index=False)
        print(f'SAVED: {collected_EntInfo_outfile}')
    ###################################################################################################

    ###################################################################################################

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

    # Define a file to save the intermeditate dictionary too for quicker loading
    DATA_outfile = os.path.join(data_path, f'DATA.pkl')
    logging.info(f'DATA_outfile: {DATA_outfile}')
    
    logging.info(f'No DATA.pkl file found. Generating from scratch.')
    # Step 1: load the G and Q data
    Quench_df, Native_df = anal.load_OP()

    # Step 2: for each candidate get the native G and Q states
    anal.get_nativeSIM_stats(Native_df)
    
    # Step 3: Identify the native state in the MSM data
    Quench_df = anal.find_native_state(Quench_df)

    # Step 4: Get the clustered EntInfo files
    anal.get_ClusterEntInfoFiles(Quench_df)

    print(f'logfile: {logfile}')

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()

print(f'NORMAL TERMINATION: {time.time() - start_time}')
logging.info(f'NORMAL TERMINATION: {time.time() - start_time}')