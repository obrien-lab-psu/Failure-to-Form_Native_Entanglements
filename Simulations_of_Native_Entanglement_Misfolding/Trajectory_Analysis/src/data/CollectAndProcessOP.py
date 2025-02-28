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
from multiprocessing import Pool, cpu_count
from scipy.ndimage import label
from types import SimpleNamespace
#pd.set_option('display.max_rows', 5000)

class Analysis:
    """
    A class to collect the G, Q, K, Zeta, and clustered EntInfo data into a single binary dictionary that can then be read into other analysis
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
        ("--gene", type=str, required=True, help="gene.")
        """

        # parse the parameters 
        self.candidates = pd.read_csv(args.candidates)
        logging.info(f'candidates:\n{self.candidates}')

        self.gene = args.gene
        logging.info(f'gene: {self.gene}')

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

        self.downsample = args.downsample
        print(f'downsample: {self.downsample}')

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

        self.ZetaFiles = glob.glob(os.path.join(self.toplevel, '*/Zeta/*_quench.Zeta'))
        print(f'Number of Zeta files found: {len(self.ZetaFiles)}')


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
        self.candidates = self.candidates[self.candidates['gene'] == self.gene]
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
            ## Load Quench Zeta data and add it to the df
            ZFiles = [f for f in self.ZetaFiles if tag in f]
            print(f'ZFiles: {len(ZFiles)}')
            df['Z'] = np.nan #initialize the column with nan values
            for f in ZFiles:
                print(f)
                Zeta = pd.read_csv(f)
                #print(Zeta)
                for traj, traj_df in Zeta.groupby('traj'):
                    frames = df.loc[(df['traj'] == int(traj)),'frame'].values
                    Zeta_frames = traj_df[traj_df['frame'].isin(frames)] ## only get those frames of Zeta that were found in the MSM datafile. This depends on the setID 
                    #print(Zeta_frames)
                    df.loc[(df['traj'] == int(traj)), 'Z'] = Zeta_frames['Zeta'].values
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
        Native_outfile = os.path.join(self.data_path, f'{self.outname}_Native_Collected_GQK.csv')
        Native_dfs.to_csv(Native_outfile, index=False)
        print(f'SAVED: {Native_outfile}')

        return Quench_dfs, Native_dfs
    #######################################################################################  

    #######################################################################################
    def get_nativeSIM_stats(self, Native_df):
        """
        Get the following parameters for Q, G, and K in the native state simulations
        median 
        mean (ci about the mean)
        std in Q
        <Qmode> across all 15ns sliding windows (ci of this distribution of modes)
        threshold = <Qmode> - 3*sigma
        """
        #print(self.candidates)
        self.stats_df = {'gene':[], 'pdb':[], 'chain':[], 'n':[],
                         'Qmean':[], 'Qstd':[], 'Q_lb':[], 'Q_ub':[], 'meanQmode':[], 'stdQmode':[], 'Qmode_lb':[], 'Qmode_ub':[], 'Qthreshold':[],
                        'Gmean':[], 'Gstd':[], 'G_lb':[], 'G_ub':[], 'meanGmode':[], 'stdGmode':[], 'Gmode_lb':[], 'Gmode_ub':[], 'Gthreshold':[], 
                        'Kmean':[], 'Kstd':[], 'K_lb':[], 'K_ub':[], 'meanKmode':[], 'stdKmode':[], 'Kmode_lb':[], 'Kmode_ub':[], 'Kthreshold':[]}
        Native_stats_outfile = os.path.join(self.data_path, f'{self.outname}_Native_Collected_GQK_STATS.csv')
        if not os.path.exists(Native_stats_outfile):
            for gene, pdb, chain in self.candidates[['gene', 'pdb', 'chain']].values:
                tag = f'{gene}_{pdb}_{chain}'
                print(tag)
                gene_df = Native_df[Native_df['gene'] == gene]
                print(gene_df)

                Qdata = [traj_df['Q'].values[::self.downsample] for traj, traj_df in gene_df.groupby('traj')]
                Gdata = [traj_df['G'].values[::self.downsample] for traj, traj_df in gene_df.groupby('traj')]
                Kdata = [traj_df['K'].values[::self.downsample] for traj, traj_df in gene_df.groupby('traj')]       

                Qmean, Qstd, Q_lb, Q_ub, meanQmode, stdQmode, Qmode_lb, Qmode_ub = get_stats(Qdata, reps=10000)
                Qthreshold = meanQmode - 3*stdQmode
                Gmean, Gstd, G_lb, G_ub, meanGmode, stdGmode, Gmode_lb, Gmode_ub  = get_stats(Gdata, reps=10000)
                Gthreshold = meanGmode + 3*stdGmode
                Kmean, Kstd, K_lb, K_ub, meanKmode, stdKmode, Kmode_lb, Kmode_ub  = get_stats(Kdata, reps=10000)
                Kthreshold = meanKmode - 3*stdKmode
                #print(f'<<Q>> = {Qmean} | std(<Q>) = {Qstd} | 95%ci = ({Q_lb}, {Q_ub}) | <Qmode> = {meanQmode} | std(Qmode) = {stdQmode} | 95%ci = ({Qmode_lb}, {Qmode_ub}) | Qthreshold = {Qthreshold}')
                #print(f'<<G>> = {Gmean} | std(<G>) = {Gstd} | 95%ci = ({G_lb}, {G_ub}) | <Gmode> = {meanGmode} | std(Gmode) = {stdGmode} | 95%ci = ({Gmode_lb}, {Gmode_ub}) | Gthreshold = {Gthreshold}')
                #print(f'<<K>> = {Kmean} | std(<K>) = {Kstd} | 95%ci = ({K_lb}, {K_ub}) | <Kmode> = {meanKmode} | std(Kmode) = {stdKmode} | 95%ci = ({Kmode_lb}, {Kmode_ub}) | Kthreshold = {Kthreshold}')

                self.stats_df['gene'] += [gene]
                self.stats_df['pdb'] += [pdb]
                self.stats_df['chain'] += [chain]
                self.stats_df['n'] += [len(Qdata)]
                self.stats_df['Qmean'] += [Qmean]
                self.stats_df['Qstd'] += [Qstd]
                self.stats_df['Q_lb'] += [Q_lb]
                self.stats_df['Q_ub'] += [Q_ub]
                self.stats_df['meanQmode'] += [meanQmode]
                self.stats_df['stdQmode'] += [stdQmode]
                self.stats_df['Qmode_lb'] += [Qmode_lb]
                self.stats_df['Qmode_ub'] += [Qmode_ub]
                self.stats_df['Qthreshold'] += [Qthreshold]
                self.stats_df['Gmean'] += [Gmean]
                self.stats_df['Gstd'] += [Gstd]
                self.stats_df['G_lb'] += [G_lb]
                self.stats_df['G_ub'] += [G_ub]
                self.stats_df['meanGmode'] += [meanGmode]
                self.stats_df['stdGmode'] += [stdGmode]
                self.stats_df['Gmode_lb'] += [Gmode_lb]
                self.stats_df['Gmode_ub'] += [Gmode_ub]
                self.stats_df['Gthreshold'] += [Gthreshold]
                self.stats_df['Kmean'] += [Kmean]
                self.stats_df['Kstd'] += [Kstd]
                self.stats_df['K_lb'] += [K_lb]
                self.stats_df['K_ub'] += [K_ub]
                self.stats_df['meanKmode'] += [meanKmode]
                self.stats_df['stdKmode'] += [stdKmode]
                self.stats_df['Kmode_lb'] += [Kmode_lb]
                self.stats_df['Kmode_ub'] += [Kmode_ub]
                self.stats_df['Kthreshold'] += [Kthreshold]

            self.stats_df = pd.DataFrame(self.stats_df)
            print(f'Native QGK stats:\n{self.stats_df}')
            Native_stats_outfile = os.path.join(self.data_path, f'{self.outname}_Native_Collected_GQK_STATS.csv')
            self.stats_df.to_csv(Native_stats_outfile, index=False)
            print(f'SAVED: {Native_stats_outfile}')
        else:
            self.stats_df = pd.read_csv(Native_stats_outfile)
            print(f'LOADED: {Native_stats_outfile}')
            print(f'Native QGK stats:\n{self.stats_df}')
    #######################################################################################

    #######################################################################################
    def get_quenchSIM_stats(self, df):
        """
        For each candidate find the native state amongst the MSM data defined by the highest <Q> and lowest <G>
        Tag this state in a new column called NativeByMSM in the self.data[tag]['MSM'] dataframe
        Also find native frames by using the Q >= <Q>ns - 3*sigma and G <= <G>ns + 3*sigma and tag a column named NativeByRef in the self.data[tag]['MSM'] dataframe

        using a 15ns sliding window calculate the Qmode and Gmode of all windows in the trajectory (sliding)
        
        """
        print(f'Find the native states by Ref')
        ###################################################################################################
        Quench_stats_outfile = os.path.join(self.data_path, f'{self.outname}_Quench_Collected_Traj_GQK_STATS.csv')
        if not os.path.exists(Quench_stats_outfile):
            self.traj_stats_df = {'gene':[], 'pdb':[], 'chain':[], 'traj':[],'n':[], 'NativeFolded':[],
                        'Qmean':[], 'Qstd':[], 'Q_lb':[], 'Q_ub':[], 'meanQmode':[], 'stdQmode':[], 'Qmode_lb':[], 'Qmode_ub':[], 
                    'Gmean':[], 'Gstd':[], 'G_lb':[], 'G_ub':[], 'meanGmode':[], 'stdGmode':[], 'Gmode_lb':[], 'Gmode_ub':[], 
                    'Kmean':[], 'Kstd':[], 'K_lb':[], 'K_ub':[], 'meanKmode':[], 'stdKmode':[], 'Kmode_lb':[], 'Kmode_ub':[],
                    'Zmean':[], 'Zstd':[], 'Z_lb':[], 'Z_ub':[], 'meanZmode':[], 'stdZmode':[], 'Zmode_lb':[], 'Zmode_ub':[]}
            
            for gene, pdb, chain in self.candidates[['gene', 'pdb', 'chain']].values:
                tag = f'{gene}_{pdb}_{chain}'
                gene_df = df[df['gene'] == gene]

                #######################################################
                ### Identify native frames by the reference simulations
                Qthreshold = self.stats_df.loc[(self.stats_df['gene'] == gene), 'Qthreshold'].values[0]
                Gthreshold = self.stats_df.loc[(self.stats_df['gene'] == gene), 'Gthreshold'].values[0]
                Kthreshold = self.stats_df.loc[(self.stats_df['gene'] == gene), 'Kthreshold'].values[0]
                print(f'{tag} Qthreshold: {Qthreshold} | Gthreshold: {Gthreshold} | Kthreshold: {Kthreshold}')

                #######################################################
                ### Identify native folded trajectories by the reference simulations
                for traj, traj_df in gene_df.groupby('traj'):
                    #print(traj_df)
                    
                    if self.setID == 3:
                        Qdata = traj_df['Q'].values[::self.downsample]
                        Gdata = traj_df['G'].values[::self.downsample]
                        Kdata = traj_df['K'].values[::self.downsample]
                        Zdata = traj_df['Z'].values[::self.downsample]
                    else:
                        Qdata = traj_df.tail(2667)['Q'].values[::self.downsample]
                        Gdata = traj_df.tail(2667)['G'].values[::self.downsample]
                        Kdata = traj_df.tail(2667)['K'].values[::self.downsample]
                        Zdata = traj_df.tail(2667)['Z'].values[::self.downsample]

                    Qmean, Qstd, Q_lb, Q_ub, meanQmode, stdQmode, Qmode_lb, Qmode_ub = get_stats([Qdata], reps=10000)
                    Gmean, Gstd, G_lb, G_ub, meanGmode, stdGmode, Gmode_lb, Gmode_ub  = get_stats([Gdata], reps=10000)
                    Kmean, Kstd, K_lb, K_ub, meanKmode, stdKmode, Kmode_lb, Kmode_ub  = get_stats([Kdata], reps=10000)
                    Zmean, Zstd, Z_lb, Z_ub, meanZmode, stdZmode, Zmode_lb, Zmode_ub  = get_stats([Zdata], reps=10000)
                    #print(f'{gene} {traj} | <<Q>> = {Qmean} | std(<Q>) = {Qstd} | 95%ci = ({Q_lb}, {Q_ub}) | <Qmode> = {meanQmode} | std(Qmode) = {stdQmode} | 95%ci = ({Qmode_lb}, {Qmode_ub}) | Qthreshold = {Qthreshold}')
                    #print(f'{gene} {traj} | <<G>> = {Gmean} | std(<G>) = {Gstd} | 95%ci = ({G_lb}, {G_ub}) | <Gmode> = {meanGmode} | std(Gmode) = {stdGmode} | 95%ci = ({Gmode_lb}, {Gmode_ub}) | Gthreshold = {Gthreshold}')
                    #print(f'{gene} {traj} | <<K>> = {Kmean} | std(<K>) = {Kstd} | 95%ci = ({K_lb}, {K_ub}) | <Kmode> = {meanKmode} | std(Kmode) = {stdKmode} | 95%ci = ({Kmode_lb}, {Kmode_ub}) | Kthreshold = {Kthreshold}')

                    self.traj_stats_df['gene'] += [gene]
                    self.traj_stats_df['pdb'] += [pdb]
                    self.traj_stats_df['chain'] += [chain]
                    self.traj_stats_df['traj'] += [traj]
                    self.traj_stats_df['n'] += [len(Qdata)]
                    self.traj_stats_df['Qmean'] += [Qmean]
                    self.traj_stats_df['Qstd'] += [Qstd]
                    self.traj_stats_df['Q_lb'] += [Q_lb]
                    self.traj_stats_df['Q_ub'] += [Q_ub]
                    self.traj_stats_df['meanQmode'] += [meanQmode]
                    self.traj_stats_df['stdQmode'] += [stdQmode]
                    self.traj_stats_df['Qmode_lb'] += [Qmode_lb]
                    self.traj_stats_df['Qmode_ub'] += [Qmode_ub]
                    self.traj_stats_df['Gmean'] += [Gmean]
                    self.traj_stats_df['Gstd'] += [Gstd]
                    self.traj_stats_df['G_lb'] += [G_lb]
                    self.traj_stats_df['G_ub'] += [G_ub]
                    self.traj_stats_df['meanGmode'] += [meanGmode]
                    self.traj_stats_df['stdGmode'] += [stdGmode]
                    self.traj_stats_df['Gmode_lb'] += [Gmode_lb]
                    self.traj_stats_df['Gmode_ub'] += [Gmode_ub]
                    self.traj_stats_df['Kmean'] += [Kmean]
                    self.traj_stats_df['Kstd'] += [Kstd]
                    self.traj_stats_df['K_lb'] += [K_lb]
                    self.traj_stats_df['K_ub'] += [K_ub]
                    self.traj_stats_df['meanKmode'] += [meanKmode]
                    self.traj_stats_df['stdKmode'] += [stdKmode]
                    self.traj_stats_df['Kmode_lb'] += [Kmode_lb]
                    self.traj_stats_df['Kmode_ub'] += [Kmode_ub]
                    self.traj_stats_df['Zmean'] += [Zmean]
                    self.traj_stats_df['Zstd'] += [Zstd]
                    self.traj_stats_df['Z_lb'] += [Z_lb]
                    self.traj_stats_df['Z_ub'] += [Z_ub]
                    self.traj_stats_df['meanZmode'] += [meanZmode]
                    self.traj_stats_df['stdZmode'] += [stdZmode]
                    self.traj_stats_df['Zmode_lb'] += [Zmode_lb]
                    self.traj_stats_df['Zmode_ub'] += [Zmode_ub]
                    
                    if meanQmode >= Qthreshold and meanGmode <= Gthreshold:
                        self.traj_stats_df['NativeFolded'] += [True]
                    else:
                        self.traj_stats_df['NativeFolded'] += [False]
                #######################################################

            self.traj_stats_df = pd.DataFrame(self.traj_stats_df)
            print(f'Traj QGK stats:\n{self.traj_stats_df}')
            Quench_stats_outfile = os.path.join(self.data_path, f'{self.outname}_Quench_Collected_Traj_GQK_STATS.csv')
            self.traj_stats_df.to_csv(Quench_stats_outfile, index=False)
            print(f'SAVED: {Quench_stats_outfile}')
        else:
            self.traj_stats_df = pd.read_csv(Quench_stats_outfile)
            print(f'LOADED: {Quench_stats_outfile}')
            print(f'Quench_Collected_GQK_STATS:\n{self.traj_stats_df}')
    ###################################################################################################

    ###################################################################################################
    def find_native_MSS(self, df):
        print(f'Find the native MSS')
        ##################################################################################################

        Quench_mss_stats_outfile = os.path.join(self.data_path, f'{self.outname}_Quench_Collected_MSS_GQK_STATS.csv')
        if not os.path.exists(Quench_mss_stats_outfile):
            self.mss_stats_df = {'gene':[], 'pdb':[], 'chain':[], 'mssID':[],'n':[], 'NativeFolded':[],
                        'Qmean':[], 'Qstd':[], 'Q_lb':[], 'Q_ub':[], 'Qmode':[], 'stdQmode':[], 'Qmode_lb':[], 'Qmode_ub':[], 
                    'Gmean':[], 'Gstd':[], 'G_lb':[], 'G_ub':[], 'Gmode':[], 'stdGmode':[], 'Gmode_lb':[], 'Gmode_ub':[], 
                    'Kmean':[], 'Kstd':[], 'K_lb':[], 'K_ub':[], 'Kmode':[], 'stdKmode':[], 'Kmode_lb':[], 'Kmode_ub':[],
                    'Zmean':[], 'Zstd':[], 'Z_lb':[], 'Z_ub':[], 'Zmode':[], 'stdZmode':[], 'Zmode_lb':[], 'Zmode_ub':[]}
            
            for gene, pdb, chain in self.candidates[['gene', 'pdb', 'chain']].values:
                tag = f'{gene}_{pdb}_{chain}'

                #######################################################
                ### Identify native frames by the reference simulations
                Qthreshold = self.stats_df.loc[(self.stats_df['gene'] == gene), 'Qthreshold'].values[0]
                Gthreshold = self.stats_df.loc[(self.stats_df['gene'] == gene), 'Gthreshold'].values[0]
                Kthreshold = self.stats_df.loc[(self.stats_df['gene'] == gene), 'Kthreshold'].values[0]
                print(f'{tag} Qthreshold: {Qthreshold} | Gthreshold: {Gthreshold} | Kthreshold: {Kthreshold}')

                gene_df = df[df['gene'] == gene]
                #print(gene_df)
  
                for mssID, mss_df in gene_df.groupby(['metastablestate']):
                    mssID = mssID[0]
                    #print(mssID, mss_df)
                    if self.setID == 3:
                        Qdata = mss_df['Q'].values[::self.downsample]
                        Gdata = mss_df['G'].values[::self.downsample]
                        Kdata = mss_df['K'].values[::self.downsample]
                        Zdata = mss_df['Z'].values[::self.downsample]
                    else:
                        Qdata = mss_df.tail(2667)['Q'].values[::self.downsample]
                        Gdata = mss_df.tail(2667)['G'].values[::self.downsample]
                        Kdata = mss_df.tail(2667)['K'].values[::self.downsample]
                        Zdata = mss_df.tail(2667)['Z'].values[::self.downsample]

                    Qmean, Qstd, Q_lb, Q_ub, Qmode, stdQmode, Qmode_lb, Qmode_ub = get_mss_stats(Qdata, reps=10000)
                    Gmean, Gstd, G_lb, G_ub, Gmode, stdGmode, Gmode_lb, Gmode_ub  = get_mss_stats(Gdata, reps=10000)
                    Kmean, Kstd, K_lb, K_ub, Kmode, stdKmode, Kmode_lb, Kmode_ub  = get_mss_stats(Kdata, reps=10000)
                    Zmean, Zstd, Z_lb, Z_ub, Zmode, stdZmode, Zmode_lb, Zmode_ub  = get_mss_stats(Zdata, reps=10000)
                    #print(f'{gene} mss{mssID} | <<Q>> = {Qmean} | std(<Q>) = {Qstd} | 95%ci = ({Q_lb}, {Q_ub}) | <Qmode> = {Qmode} | std(Qmode) = {stdQmode} | 95%ci = ({Qmode_lb}, {Qmode_ub}) | Qthreshold = {Qthreshold}')
                    #print(f'{gene} mss{mssID} | <<G>> = {Gmean} | std(<G>) = {Gstd} | 95%ci = ({G_lb}, {G_ub}) | <Gmode> = {Gmode} | std(Gmode) = {stdGmode} | 95%ci = ({Gmode_lb}, {Gmode_ub}) | Gthreshold = {Gthreshold}')
                    #print(f'{gene} mss{mssID} | <<K>> = {Kmean} | std(<K>) = {Kstd} | 95%ci = ({K_lb}, {K_ub}) | <Kmode> = {Kmode} | std(Kmode) = {stdKmode} | 95%ci = ({Kmode_lb}, {Kmode_ub}) | Kthreshold = {Kthreshold}')

                    self.mss_stats_df['gene'] += [gene]
                    self.mss_stats_df['pdb'] += [pdb]
                    self.mss_stats_df['chain'] += [chain]
                    self.mss_stats_df['mssID'] += [mssID]
                    self.mss_stats_df['n'] += [len(Qdata)]
                    self.mss_stats_df['Qmean'] += [Qmean]
                    self.mss_stats_df['Qstd'] += [Qstd]
                    self.mss_stats_df['Q_lb'] += [Q_lb]
                    self.mss_stats_df['Q_ub'] += [Q_ub]
                    self.mss_stats_df['Qmode'] += [Qmode]
                    self.mss_stats_df['stdQmode'] += [stdQmode]
                    self.mss_stats_df['Qmode_lb'] += [Qmode_lb]
                    self.mss_stats_df['Qmode_ub'] += [Qmode_ub]
                    self.mss_stats_df['Gmean'] += [Gmean]
                    self.mss_stats_df['Gstd'] += [Gstd]
                    self.mss_stats_df['G_lb'] += [G_lb]
                    self.mss_stats_df['G_ub'] += [G_ub]
                    self.mss_stats_df['Gmode'] += [Gmode]
                    self.mss_stats_df['stdGmode'] += [stdGmode]
                    self.mss_stats_df['Gmode_lb'] += [Gmode_lb]
                    self.mss_stats_df['Gmode_ub'] += [Gmode_ub]
                    self.mss_stats_df['Kmean'] += [Kmean]
                    self.mss_stats_df['Kstd'] += [Kstd]
                    self.mss_stats_df['K_lb'] += [K_lb]
                    self.mss_stats_df['K_ub'] += [K_ub]
                    self.mss_stats_df['Kmode'] += [Kmode]
                    self.mss_stats_df['stdKmode'] += [stdKmode]
                    self.mss_stats_df['Kmode_lb'] += [Kmode_lb]
                    self.mss_stats_df['Kmode_ub'] += [Kmode_ub]
                    self.mss_stats_df['Zmean'] += [Zmean]
                    self.mss_stats_df['Zstd'] += [Zstd]
                    self.mss_stats_df['Z_lb'] += [Z_lb]
                    self.mss_stats_df['Z_ub'] += [Z_ub]
                    self.mss_stats_df['Zmode'] += [Zmode]
                    self.mss_stats_df['stdZmode'] += [stdZmode]
                    self.mss_stats_df['Zmode_lb'] += [Zmode_lb]
                    self.mss_stats_df['Zmode_ub'] += [Zmode_ub]
                    
                    if Qmode >= Qthreshold and Gmode <= Gthreshold:
                        self.mss_stats_df['NativeFolded'] += [True]
                    else:
                        self.mss_stats_df['NativeFolded'] += [False]
                #######################################################

            self.mss_stats_df = pd.DataFrame(self.mss_stats_df)
            print(f'MSS QGK stats:\n{self.mss_stats_df}')
            Quench_mss_stats_outfile = os.path.join(self.data_path, f'{self.outname}_Quench_Collected_MSS_GQK_STATS.csv')
            self.mss_stats_df.to_csv(Quench_mss_stats_outfile, index=False)
            print(f'SAVED: {Quench_mss_stats_outfile}')
        else:
            self.mss_stats_df = pd.read_csv(Quench_mss_stats_outfile)
            print(f'LOADED: {Quench_mss_stats_outfile}')
            print(f'Quench_Collected_MSS_GQK_STATS:\n{self.mss_stats_df}')
        ###################################################################################################
    ###################################################################################################

    ###################################################################################################
    def get_ClusterEntInfoFiles(self, df):
        """
        Collect the change in entanglement info files and grab those frames and traj that have already been controlled for mirros
        """
        keys = ['TotalChanges', 'TotalLoss', 'TotalGain', 'TrueLoss', 'PartialLoss', 'TrueGain', 'PartialGain', 'TrueLossOverlap', 'PartialLossOverlap', 'TrueGainOverlap', 'PartialGainOverlap']
        Quench_outfile = os.path.join(self.data_path, f'{self.outname}_Quench_Collected_GQKZ_wEntInfo.csv')
        if not os.path.exists(Quench_outfile):
            for key in keys:
                df[key] = 0

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
                    else:
                        EntInfo_f = EntInfo_f[0]
                    
                    EntInfo = pd.read_csv(EntInfo_f, low_memory=False)
                    EntInfo = EntInfo[EntInfo['Frame'].isin(frames)]
                    EntInfo['gene'] = gene
                    EntInfo['traj'] = traj
                    EntInfo['crossingsN'] = EntInfo['crossingsN'].astype(str)
                    EntInfo['crossingsC'] = EntInfo['crossingsC'].astype(str)
                    print(EntInfo)

                    # quality check if there are no EntInfo for the last 200ns
                    if len(EntInfo) == 0:
                        continue
                    
                    results = parallel_process_genes(EntInfo)   
                    results = pd.concat(results)
                    #print(f'results:\n{results}')
    

                    # Merging DataFrame B with the L column from DataFrame A based on the G column
                    #df = df.merge(results[["gene", "traj", 'frame', "TotalChanges"]], on=["gene", "traj", 'frame'], how="left")
                    for rowi, row in results.iterrows():
                        frame = row['frame']
                        for key in keys:
                            df.loc[(df['gene'] == gene) & (df['traj'] == traj) & (df['frame'] == frame), key] = row[key]

            print(f'Quench_Collected_GQKZ_wEntInfo updated:\n{df}')
            Quench_outfile = os.path.join(self.data_path, f'{self.outname}_Quench_Collected_GQKZ_wEntInfo.csv')
            df.to_csv(Quench_outfile, index=False)
            print(f'SAVED: {Quench_outfile}')
        else:
            df = pd.read_csv(Quench_outfile)
            print(f'LOADED: {Quench_outfile}')
            print(f'Quench_Collected_GQKZ_wEntInfo:\n{df}')

        return df
    ###################################################################################################

    ###################################################################################################
    def update_GQK_df(self, df):
        df['NativeByRef'] = False
        df['NativeByMSS'] = False
        df['NativeFolded'] = False
        print(df)
        ###################################################################################################
        Quench_outfile = os.path.join(self.data_path, f'{self.outname}_Quench_Collected_GQK.csv')
        if not os.path.exists(Quench_outfile):
            
            for gene, pdb, chain in self.candidates[['gene', 'pdb', 'chain']].values:
                tag = f'{gene}_{pdb}_{chain}'
                gene_df = df[df['gene'] == gene]
                print(gene_df)

                #######################################################
                ### Identify native frames by the reference simulations
                Qthreshold = self.stats_df.loc[(self.stats_df['gene'] == gene), 'Qthreshold'].values[0]
                Gthreshold = self.stats_df.loc[(self.stats_df['gene'] == gene), 'Gthreshold'].values[0]
                Kthreshold = self.stats_df.loc[(self.stats_df['gene'] == gene), 'Kthreshold'].values[0]
                print(f'{tag} Qthreshold: {Qthreshold} | Gthreshold: {Gthreshold} | Kthreshold: {Kthreshold}')
                #######################################################

                ## determine NativeByMSS self.mss_stats_df
                mss_stats = self.mss_stats_df[(self.mss_stats_df['gene'] == gene)]
                print(mss_stats)
                native_mssID = mss_stats[mss_stats['NativeFolded'] == True]
                if len(native_mssID) == 1:
                    native_mssID = native_mssID['mssID'].values[0]
                    
                elif len(native_mssID) > 1:
                    native_mssID = native_mssID.sort_values(by=['n'], ascending=False)
                    print(native_mssID)
                    native_mssID = int(native_mssID['mssID'].values[0])

                else:
                    continue
    
                print(f'native_mssID: {native_mssID}')
                df.loc[(df['gene'] == gene) & (df['metastablestate'] == native_mssID), 'NativeByMSS'] = True
                #print(gene_df)
                #######################################################

                #######################################################
                ### Identify native folded trajectories by the reference simulations
                for traj, traj_df in gene_df.groupby('traj'):
                    #print(traj_df)

                    ## determine NativeByRef
                    df.loc[(df['gene'] == gene) & (df['traj'] == traj) & (df['Q'] >= Qthreshold) & (df['G'] <= Gthreshold), 'NativeByRef'] = True
                    #print(traj_df)

                    ## determine NativeFolded
                    traj_stats = self.traj_stats_df[(self.traj_stats_df['gene'] == gene) & (self.traj_stats_df['traj'] == traj)]
                    if len(traj_stats) != 1:
                        raise ValueError(f'NO trajectory stats found')
                    NativeFolded = traj_stats['NativeFolded'].values[0]
                    #print(traj_stats, NativeFolded)
                    df.loc[(df['gene'] == gene) & (df['traj'] == traj), 'NativeFolded'] = NativeFolded


            print(f'QGK updated:\n{df}')
            Quench_outfile = os.path.join(self.data_path, f'{self.outname}_Quench_Collected_GQK.csv')
            df.to_csv(Quench_outfile, index=False)
            print(f'SAVED: {Quench_outfile}')
        else:
            df = pd.read_csv(Quench_outfile)
            print(f'LOADED: {Quench_outfile}')
            print(f'Quench_Collected_GQK:\n{df}')

        return df
    ###################################################################################################

    #######################################################################################
    def Lifetimes(self, df):
        """
        Caluclates the following metrics and confidence intervals 
        Lifetimes of gain or loss as a function of Q bins

        """
        print(f'Lifetimes stats from df')
        #print(f'Mechanism Lifetimes from df:\n{df}')
        #print(f'self.quench_traj_stats_df:\n{self.quench_traj_stats_df}')
        print(f'df: {df}')
        keys = ['TotalChanges', 'TotalLoss', 'TotalGain', 'TrueLoss', 'PartialLoss', 'TrueGain', 'PartialGain', 'TrueLossOverlap', 'PartialLossOverlap', 'TrueGainOverlap', 'PartialGainOverlap']

        Mechanism_Lifetimes_summary_outfile = os.path.join(self.data_path, f'{self.outname}_ChangeType_Lifetimes_summary_setID{self.setID}.csv')
        if not os.path.exists(Mechanism_Lifetimes_summary_outfile):
            # Define some datastructures
            column = 'meanQmode'
            n_resamples = 10000

            outdf = {'gene':[], 'traj':[], 'ChangeType':[], 
                    'meanLifetime':[], 'meanLifetime_lb':[], 'meanLifetime_ub':[],
                    'meanFract':[], 'meanFract_lb':[], 'meanFract_ub':[]}
            #df = df[df['traj'] == 26]

            for traj, traj_df in df.groupby('traj'):
                #print(traj_df)
                #print(traj_df['NativeFolded'])
                
                Fract_df = traj_df[traj_df['TotalChanges'] != 0]
                for key in keys:
                    outdf['gene'] += [self.gene]
                    outdf['traj'] += [traj]
                    outdf['ChangeType'] += [key]
                    
                    ## Make sure to only consider frames with both gain and loss when doing overlap conditional probabilities
                    if 'Overlap' in key:
                        #print(key)
                        both_frames = Fract_df[(Fract_df['TotalLoss'] != 0) & (Fract_df['TotalGain'] != 0)]
                        #print(both_frames)
                        if 'TrueLoss' in key:
                            Fract_arr = both_frames[key]/both_frames['TrueLoss']
                        elif 'PartialLoss' in key:
                            Fract_arr = both_frames[key]/both_frames['PartialLoss']
                        if 'TrueGain' in key:
                            Fract_arr = both_frames[key]/both_frames['TrueGain']
                        elif 'PartialGain' in key:
                            Fract_arr = both_frames[key]/both_frames['PartialGain']
                    else:
                        Fract_arr = Fract_df[key]/Fract_df['TotalChanges']
                    Fract_arr = Fract_arr.values

                    ## get the fraction of the changetype and its stats
                    mean = np.nanmean(Fract_arr)
                    std = np.std(Fract_arr)
                    lb, ub = bootstrap(Fract_arr, reps=n_resamples)
                    
                    #median, mean, std, lb, ub = get_mss_stats(lifetimes, reps=n_resamples)
                    print(f'{key} Fractions: {mean} ({lb}, {ub})')  
                    outdf['meanFract'] += [mean]
                    outdf['meanFract_lb'] += [lb]
                    outdf['meanFract_ub'] += [ub]

                    ## get the mean lifetime of the changetype and its stats
                    array = traj_df[key].values
                    binary_array = np.where(array > 0, 1, 0)

                    labels_ones, num_features_ones = label(binary_array)
                    lifetimes = [np.sum(labels_ones == i)*0.075 for i in range(1, num_features_ones + 1)]
                    #print(f'lifetimes: {lifetimes}')

                    mean = np.mean(lifetimes)
                    std = np.std(lifetimes)
                    lb, ub = bootstrap(lifetimes, reps=n_resamples)
                    
                    print(f'{key} lifetimes: {mean} ({lb}, {ub})')  
                    outdf['meanLifetime'] += [mean]
                    outdf['meanLifetime_lb'] += [lb]
                    outdf['meanLifetime_ub'] += [ub]

            outdf = pd.DataFrame(outdf)
            print(f'outdf:\n{outdf.to_string()}')
            ## save the Mechanism stats summary file
            outdf.to_csv(Mechanism_Lifetimes_summary_outfile, index=False)
            logging.info(f'SAVED: {Mechanism_Lifetimes_summary_outfile}')
            print(f'SAVED: {Mechanism_Lifetimes_summary_outfile}')

        else:
            Mechanism_Lifetimes_summary_outfile = os.path.join(self.data_path, f'{self.outname}_ChangeType_Lifetimes_summary_setID{self.setID}.csv')
            outdf = pd.read_csv(Mechanism_Lifetimes_summary_outfile)
            logging.info(f'LOADED: {Mechanism_Lifetimes_summary_outfile}')
            print(f'LOADED: {Mechanism_Lifetimes_summary_outfile}')
    #######################################################################################

    #######################################################################################
    def StructureCat(self, df):
        """
        Caluclates the following metrics and confidence intervals 
        Lifetimes of gain or loss as a function of Q bins

        """
        print(f'Fraction of structures with onlyLoss, onlyGain, BothLossGain')
        #print(f'Mechanism Lifetimes from df:\n{df}')
        #print(f'self.quench_traj_stats_df:\n{self.quench_traj_stats_df}')
        df = df[['gene', 'traj', 'TotalChanges', 'TotalLoss', 'TotalGain']]
        df = df[df['TotalChanges'] != 0]
        print(f'df:\n{df}')

        Mechanism_StructureCategory_summary_outfile = os.path.join(self.data_path, f'{self.outname}_StructureCategory_summary_setID{self.setID}.csv')
        if not os.path.exists(Mechanism_StructureCategory_summary_outfile):
            # Define some datastructures
            n_resamples = 10000

            outdf = {'gene':[], 'traj':[], 'ChangeType':[], 
                    'mean':[], 'lb':[], 'ub':[]}

            LossOnly = []
            GainOnly = []
            BothLossGain = []
            for traj, traj_df in df.groupby('traj'):
                #print(traj_df)

                N = len(traj_df)
                LossOnly_n = len(traj_df[(traj_df['TotalLoss'] != 0) & (traj_df['TotalGain'] == 0)])
                GainOnly_n = len(traj_df[(traj_df['TotalLoss'] == 0) & (traj_df['TotalGain'] != 0)])
                BothLossGain_n = len(traj_df[(traj_df['TotalLoss'] != 0) & (traj_df['TotalGain'] != 0)])

                LossOnly = [1] * LossOnly_n + [0] * (N - LossOnly_n)
                GainOnly = [1] * GainOnly_n + [0] * (N - GainOnly_n)
                BothLossGain = [1] * BothLossGain_n + [0] * (N - BothLossGain_n)

                for ChangeType, arr in [('LossOnly', LossOnly), ('GainOnly', GainOnly), ('BothLossGain', BothLossGain)]:
                    ## get the fraction of the changetype and its stats
                    mean = np.nanmean(arr)
                    std = np.std(arr)
                    lb, ub = bootstrap(arr, reps=n_resamples)
                    
                    #median, mean, std, lb, ub = get_mss_stats(lifetimes, reps=n_resamples)
                    print(f'{ChangeType} Fraction structures: {mean} ({lb}, {ub})')  
                    outdf['gene'] += [self.gene]
                    outdf['traj'] += [traj]
                    outdf['ChangeType'] += [ChangeType]
                    outdf['mean'] += [mean]
                    outdf['lb'] += [lb]
                    outdf['ub'] += [ub]

            outdf = pd.DataFrame(outdf)
            print(f'outdf:\n{outdf.to_string()}')
            ## save the Mechanism stats summary file
            outdf.to_csv(Mechanism_StructureCategory_summary_outfile, index=False)
            logging.info(f'SAVED: {Mechanism_StructureCategory_summary_outfile}')
            print(f'SAVED: {Mechanism_StructureCategory_summary_outfile}')

        else:
            outdf = pd.read_csv(Mechanism_StructureCategory_summary_outfile)
            logging.info(f'LOADED: {Mechanism_StructureCategory_summary_outfile}')
            print(f'LOADED: {Mechanism_StructureCategory_summary_outfile}')
    #######################################################################################

####################################################
def parallel_process_genes(combined_EntInfo):
    #print(combined_EntInfo)
    # Prepare arguments for parallel processing
    args = [(frame, frame_df) for frame, frame_df in combined_EntInfo.groupby('Frame')]

    # Use multiprocessing to parallelize the processing
    #num_cpu = cpu_count()
    num_cpu = 20
    print(f'cpu_count: {num_cpu}')
    with Pool(num_cpu) as pool:
        results = pool.map(process_gene, args)

    # Flatten the list of results
    #flattened_results = [item for sublist in results for item in sublist]
    return results
#######################################################################################

#######################################################################################
def process_gene(args):
    frame, frame_df = args
    gene = frame_df['gene'].values[0]
    traj = frame_df['traj'].values[0]

    ### track the different types of changes in each frame
    ## total changes will be the number of unique cID in the frame
    ## True loss will be any native entanglement going to a 0 linking value
    ## partial loss will be any loss going to anything other than 0
    ## True gain will be going from a linking value of 0 to anyting else
    ## partial gain will be going from a non zero linking value to anything else
    results = {'gene':[], 'traj':[], 'frame':[], 'TotalChanges':[], 'TotalLoss':[], 'TotalGain':[],
               'TrueLoss':[], 'TrueGain':[], 'TrueLossOverlap':[], 'TrueGainOverlap':[], 'PartialLoss':[], 'PartialGain':[], 'PartialLossOverlap':[], 'PartialGainOverlap':[]}

    ## determine the overlap stats for the frame
    TrueLoss, PartialLoss, TrueLossOverlap, PartialLossOverlap, TrueGain, PartialGain, TrueGainOverlap, PartialGainOverlap = EntanglementDecomp(frame_df)
    results['gene'] += [gene]
    results['traj'] += [traj]
    results['frame'] += [frame]
    results['TotalChanges'] += [TrueLoss + PartialLoss + TrueGain + PartialGain]
    results['TotalLoss'] += [TrueLoss + PartialLoss]
    results['TotalGain'] += [TrueGain + PartialGain]
    results['TrueLoss'] += [TrueLoss]
    results['TrueGain'] += [TrueGain]
    results['TrueLossOverlap'] += [TrueLossOverlap]
    results['TrueGainOverlap'] += [TrueGainOverlap]
    results['PartialLoss'] += [PartialLoss]
    results['PartialGain'] += [PartialGain]
    results['PartialLossOverlap'] += [PartialLossOverlap]
    results['PartialGainOverlap'] += [PartialGainOverlap]
    results = pd.DataFrame(results)
    return results
#######################################################################################

#######################################################################################
def EntanglementDecomp(df):

    ## if there is only a single unique ent then determine if its a loss or gain and return
    TrueLoss = 0
    TrueGain = 0
    PartialLoss = 0
    PartialGain = 0
    TrueLossOverlap = 0
    TrueGainOverlap = 0
    PartialLossOverlap = 0
    PartialGainOverlap = 0
    #print(df)

    ## else make dictionary of residues for each unique ent to make it easier to judge overlap
    cID_res = {}
    for cID, cID_df in df.groupby('cID'):
        #print(cID_df)
        ijres = []
        Nres = []
        Cres = []
        for rowi, row in cID_df.iterrows():
            i, j = row['i'], row['j']
            ijres += [np.arange(i - 3, i + 4)]
            ijres += [np.arange(j - 3, j + 4)]

            NchangeType, CchangeType = row['NchangeType'], row['CchangeType']
            gn, gc = row['gn'], row['gc']
            delta_gn, delta_gc = row['delta_gn'], row['delta_gc']

            if NchangeType != 'NoChange':
                Nres += [np.arange(abs(float(c)) - 3, abs(float(c)) + 4) for c in row['crossingsN'].split(',')]

            if CchangeType != 'NoChange':
                Cres += [np.arange(abs(float(c)) - 3, abs(float(c)) + 4) for c in row['crossingsC'].split(',')]    

        ijres = np.hstack(ijres)
        ijres = set([int(r) for r in ijres if r > 0])
        
        if len(Nres) != 0:
            Nres = np.hstack(Nres)
            Nres = set([int(r) for r in Nres if r > 0])
            Nres = Nres.union(ijres)
        else:
            Nres = set()

        if len(Cres) != 0:
            Cres = np.hstack(Cres)
            Cres = set([int(r) for r in Cres if r > 0])
            Cres = Cres.union(ijres)
        else:
            Cres = set()
        
        cID_res[cID] = {'ijres':ijres, 'Nres':Nres, 'Ntype':NchangeType, 'gn':gn, 'delta_gn':delta_gn, 'Cres':Cres, 'Ctype':CchangeType, 'gc':gc, 'delta_gc':delta_gc}
    #########################################################

    #########################################################
    # determine True versus Partial changes and overlap
    for cID, cID_info in cID_res.items():
        #print('\n', cID, cID_info)

        for termini in ['N', 'C']:
            #############################
            # Determine if there was a true or partial loss or gain in the N or C
            if cID_info[f'{termini}type'] != 'NoChange':

                ### True Loss
                if 'Loss' in cID_info[f'{termini}type'] and int(cID_info[f'g{termini.lower()}']) == 0:
                    TrueLoss += 1
                
                    ### Check for overlap of a TrueLoss in the N or C terminus
                    for Comp_cID, Comp_cID_info in cID_res.items():
                        if Comp_cID == cID:
                            continue
                        if cID_info[f'{termini}res'].intersection(Comp_cID_info[f'{termini}res']):
                            if ('Loss' in cID_info[f'{termini}type'] and 'Gain' in Comp_cID_info[f'{termini}type']):
                                TrueLossOverlap += 1
                                break
                
                ### Partial Loss
                elif 'Loss' in cID_info[f'{termini}type'] and int(cID_info[f'g{termini.lower()}']) != 0:
                    PartialLoss += 1
                
                    ### Check for overlap of a Partial Loss in the N or Cterminus
                    for Comp_cID, Comp_cID_info in cID_res.items():
                        if Comp_cID == cID:
                            continue
                        if cID_info[f'{termini}res'].intersection(Comp_cID_info[f'{termini}res']):
                            if ('Loss' in cID_info[f'{termini}type'] and 'Gain' in Comp_cID_info[f'{termini}type']):
                                PartialLossOverlap += 1
                                break                    
            
                ### True Gain
                if 'Gain' in cID_info[f'{termini}type'] and int(abs(cID_info[f'g{termini.lower()}']) - cID_info[f'delta_g{termini.lower()}']) == 0:
                    TrueGain += 1
                
                    ### Check for overlap of a TrueGain in the N terminus
                    for Comp_cID, Comp_cID_info in cID_res.items():
                        if Comp_cID == cID:
                            continue
                        if cID_info[f'{termini}res'].intersection(Comp_cID_info[f'{termini}res']):
                            if ('Gain' in cID_info[f'{termini}type'] and 'Loss' in Comp_cID_info[f'{termini}type']):
                                TrueGainOverlap += 1   
                                break
                        
                ### Partial Gain
                elif 'Gain' in cID_info[f'{termini}type'] and int(abs(cID_info[f'g{termini.lower()}']) - cID_info[f'delta_g{termini.lower()}']) != 0:
                    PartialGain += 1

                    ### Check for overlap of a TrueGain in the N terminus
                    for Comp_cID, Comp_cID_info in cID_res.items():
                        if Comp_cID == cID:
                            continue
                        if cID_info[f'{termini}res'].intersection(Comp_cID_info[f'{termini}res']):
                            if ('Gain' in cID_info[f'{termini}type'] and 'Loss' in Comp_cID_info[f'{termini}type']):
                                PartialGainOverlap += 1   
                                break                    

        #print(TrueLoss, PartialLoss, TrueLossOverlap, PartialLossOverlap, TrueGain, PartialGain, TrueGainOverlap, PartialGainOverlap)
    return  TrueLoss, PartialLoss, TrueLossOverlap, PartialLossOverlap, TrueGain, PartialGain, TrueGainOverlap, PartialGainOverlap
#######################################################################################

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
def get_mss_stats(arr, reps=10000):
    """
    Get the <> and 95% ci for a stats array
    """
    mean = np.mean(arr)
    std = np.std(arr)
    lb, ub = bootstrap(arr, reps=reps)
   
    arr_mode = mode(arr).mode
    std_mode, (mode_lb, mode_ub) = bootstrap_mss_mode(arr, reps=reps)

    return (mean, std, lb, ub, arr_mode, std_mode, mode_lb, mode_ub)
#######################################################################################

#######################################################################################
def bootstrap_mss_mode(data, reps=10000):
    boot_modes = []
    for b in range(reps):
        boot_samp = np.random.choice(data, size=len(data))
        boot_mode = mode(boot_samp).mode
        #print(b, boot_mean)
        boot_modes += [boot_mode]

    lb = np.percentile(boot_modes, 2.5)
    ub = np.percentile(boot_modes, 97.5)
    std = np.std(boot_modes)
    return std, (lb, ub)
#######################################################################################

#######################################################################################
def meanXmode(gene_df, metric):
    modes = []
    for traj, traj_df in gene_df.groupby('traj'):
        print(traj_df)
        modes += [mode(traj_df[metric].values).mode]
    meanmode = np.mean(modes)
    varmode = np.std(modes)**2

    return meanmode, varmode
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
    parser.add_argument("--downsample", type=int, required=False, help="number of points to downsample by.", default=1)
    parser.add_argument("--gene", type=str, required=True, help="gene.")

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
   
    # Step 3: for each quench traj get the GQK stats
    anal.get_quenchSIM_stats(Quench_df)
  
    anal.find_native_MSS(Quench_df)

    Quench_df = anal.update_GQK_df(Quench_df)

    # Step 4: Get the clustered EntInfo files
    if anal.setID in [0, 3]:
        Quench_df = anal.get_ClusterEntInfoFiles(Quench_df)
        anal.Lifetimes(Quench_df)
        anal.StructureCat(Quench_df)

    print(f'logfile: {logfile}')

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()

print(f'NORMAL TERMINATION: {time.time() - start_time}')
logging.info(f'NORMAL TERMINATION: {time.time() - start_time}')