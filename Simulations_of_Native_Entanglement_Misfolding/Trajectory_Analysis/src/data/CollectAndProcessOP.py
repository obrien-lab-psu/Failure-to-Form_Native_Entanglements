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

        self.nEntInfoFiles = glob.glob(os.path.join(self.toplevel, '*/Native/G/*.EntInfo'))
        print(f'Number of Native G files found: {len(self.nEntInfoFiles)}')

        self.nKFiles = glob.glob(os.path.join(self.toplevel, '*/Mirror/Native/*.dat'))
        print(f'Number of Native K files found: {len(self.KFiles)}')


        ## get the MSM files
        self.MSMFiles = glob.glob(os.path.join(self.toplevel, '*/BuildKineticModel/*_MSMmapping.csv'))
        print(f'Number of MSM files found: {len(self.KFiles)}')


        ## make Plots dir
        self.plot_path = os.path.join(self.outpath, 'Plots')
        if not os.path.exists(self.plot_path):
            os.makedirs(self.plot_path)
            print(f'Made directory: {self.plot_path}')  
        print(f'plot_path: {self.plot_path}') 
        
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

        data = {}
        print(self.candidates)
        for gene, pdb, chain in self.candidates[['gene', 'pdb', 'chain']].values:
            tag = f'{gene}_{pdb}_{chain}'
            print(f'Loading files for {gene} {pdb} {chain} {tag}')

            if tag in data:
                continue

            if tag not in data:
                data[tag] = {'Q':{'Native':{}, 'Quench':{}}, 'G':{'Native':{}, 'Quench':{}}, 'K':{'Native':{}, 'Quench':{}}, 'EntInfo':pd.DataFrame(), 
                            'Qstats':{'Native':[], 'Quench':[]}, 'Gstats':{'Native':[], 'Quench':[]}, 'Kstats':{'Native':[], 'Quench':[]}, 
                            'MSM':{'Qstats':{}, 'Gstats':{}, 'Kstats':{}}, 'df':pd.DataFrame()}

            ### Get the mirror trajs for removal
            ## remove any traj that has been identified to have a mirror
            Mirrors = self.Mirror_df[(self.Mirror_df['gene'] == gene) & (self.Mirror_df['pdb'] == pdb) & (self.Mirror_df['chain'] == chain)]
            Mirror_trajs = Mirrors['traj'].values
            print(f'Mirror_trajs: {Mirror_trajs}')

            ###########################################################################
            ## Load Quench Q, G, K data
            QFiles = [f for f in self.QFiles if tag in f]
            print(f'QFiles: {len(QFiles)}')
            for f in QFiles:
                traj = f.split('_')[-1].replace('t', '').replace('.Q', '')
                #print(f, traj)
                data[tag]['Q']['Quench'][traj] = pd.read_csv(f)['Q'].values
              
            GFiles = [f for f in self.GFiles if tag in f]
            print(f'GFiles: {len(GFiles)}')
            for f in GFiles:
                traj = f.split('_')[-1].replace('t', '').replace('.G', '')
                #print(f, traj)
                data[tag]['G']['Quench'][traj] = pd.read_csv(f)['G'].values

            KFiles = [f for f in self.KFiles if tag in f]
            print(f'KFiles: {len(KFiles)}')
            for f in KFiles:
                traj = f.split('_')[-2].replace('t', '')
                #print(f, traj)
                data[tag]['K']['Quench'][traj] = np.loadtxt(f, usecols=(-1), dtype=str)[1:].astype(float)

            EntInfoFiles = [f for f in self.ClusterEntInfoFiles if tag in f]
            print(f'EntInfoFiles: {len(EntInfoFiles)}')
            EntInfo_dfs = []
            for f in EntInfoFiles:
                traj = f.split('_')[-2].replace('t', '')
                #print(f, traj)
                if int(traj) in Mirror_trajs:
                    #print(f'Mirror traj: {traj}')
                    continue
                EntInfo_df = pd.read_csv(f, low_memory=False)
                EntInfo_df['traj'] = traj
                EntInfo_dfs += [EntInfo_df]
            EntInfo_dfs  = pd.concat(EntInfo_dfs)
            data[tag]['EntInfo'] = EntInfo_dfs

            ###########################################################################
            ## Load Native Q, G, K data
            QFiles = [f for f in self.nQFiles if tag in f]
            print(f'nQFiles: {len(QFiles)}')
            for f in QFiles:
                traj = f.split('_')[-1].replace('t', '').replace('.Q', '')
                #print(f, traj)
                data[tag]['Q']['Native'][traj] = pd.read_csv(f)['Q'].values
              
            GFiles = [f for f in self.nGFiles if tag in f]
            print(f'nGFiles: {len(GFiles)}')
            for f in GFiles:
                traj = f.split('_')[-1].replace('t', '').replace('.G', '')
                #print(f, traj)
                data[tag]['G']['Native'][traj] = pd.read_csv(f)['G'].values

            KFiles = [f for f in self.nKFiles if tag in f]
            print(f'nKFiles: {len(KFiles)}')
            for f in KFiles:
                traj = f.split('_')[-1].replace('.dat', '')
                #print(f, traj)
                data[tag]['K']['Native'][traj] = np.loadtxt(f, usecols=(-1), dtype=str)[1:].astype(float)

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
                    traj_K = data[tag]['K']['Quench'][str(traj)]
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
                df = df[~df['traj'].isin(Mirror_trajs)]
                data[tag]['df'] = df

        self.data = data
        logging.info(f'Number of genes in data: {len(self.data)}')    
        print(f'Number of genes in data: {len(self.data)}')            
    #######################################################################################  

    #######################################################################################
    def get_nativeSIM_stats(self,):
        """
        Get the <> and 95% ci for G and Q in the native state sims
        """
        #print(self.candidates)
        self.stats_df = {'gene_tag':[], 'macrostateID':[], 'desc':[], 'Qmean':[], 'Qstd':[], 'Q_lb':[], 'Q_ub':[], 
                        'Gmean':[], 'Gstd':[], 'G_lb':[], 'G_ub':[], 'Kmean':[], 'Kstd':[], 'K_lb':[], 'K_ub':[], 'n':[], 'native':[]}
        for tag, tag_data in self.data.items():
            Qdata = np.hstack(list(tag_data["Q"]['Native'].values()))
            Gdata = np.hstack(list(tag_data["G"]['Native'].values()))
            Kdata = np.hstack(list(tag_data["K"]['Native'].values()))
            print(f'{tag} with {len(Qdata)} {Qdata.shape} Q arrays | {len(Gdata)} {Gdata.shape} G arrays | {len(Kdata)} {Kdata.shape} K arrays')

            Qmedian, Qmean, Qstd, Q_lb, Q_ub = get_stats(Qdata)
            Gmedian, Gmean, Gstd, G_lb, G_ub = get_stats(Gdata)
            Kmedian, Kmean, Kstd, K_lb, K_ub = get_stats(Kdata)
            print(f'<Q> = {Qmean} | std(Q) = {Qstd} | 95%ci = ({Q_lb}, {Q_ub})')
            print(f'<G> = {Gmean} | std(G) = {Gstd} | 95%ci = ({G_lb}, {G_ub})')
            print(f'<K> = {Kmean} | std(K) = {Kstd} | 95%ci = ({K_lb}, {K_ub})')

            self.data[tag]['Qstats']['Native'] = (Qmean, Qstd, Q_lb, Q_ub)
            self.data[tag]['Gstats']['Native'] = (Gmean, Gstd, G_lb, G_ub)
            self.data[tag]['Kstats']['Native'] = (Kmean, Kstd, K_lb, K_ub)

            self.stats_df['gene_tag'] += [tag]
            self.stats_df['macrostateID'] += [-1]
            self.stats_df['desc'] += ['Ref State Sims']
            self.stats_df['Qmean'] += [Qmean]
            self.stats_df['Qstd'] += [Qstd]
            self.stats_df['Q_lb'] += [Q_lb]
            self.stats_df['Q_ub'] += [Q_ub]
            self.stats_df['Gmean'] += [Gmean]
            self.stats_df['Gstd'] += [Gstd]
            self.stats_df['G_lb'] += [G_lb]
            self.stats_df['G_ub'] += [G_ub]
            self.stats_df['Kmean'] += [Kmean]
            self.stats_df['Kstd'] += [Kstd]
            self.stats_df['K_lb'] += [K_lb]
            self.stats_df['K_ub'] += [K_ub]
            self.stats_df['n'] += [len(Qdata)]
            self.stats_df['native'] += [True]
          
    #######################################################################################

    #######################################################################################
    def find_native_state(self, ):
        """
        For each candidate find the native state amongst the MSM data defined by the highest <Q> and lowest <G>
        Tag this state in a new column called NativeByMSM in the self.data[tag]['MSM'] dataframe
        Also find native frames by using the Q >= <Q>ns - 3*sigma and G <= <G>ns + 3*sigma and tag a column named NativeByRef in the self.data[tag]['MSM'] dataframe
        """
        print(f'Find the native states by MSM and Ref')
        ###################################################################################################
        for tag, tag_data in self.data.items():
            df = tag_data['df']
            print(f'Finding Native states in {tag}')
            if len(df) != 0:
                #print(f'df {tag}\n{df}')

                local_mss_GQmeans = {'mssID':[], 'Qmean':[], 'Gmean':[]}
                for mssID, mss_df in df.groupby(['metastablestate']):
                    mssID = mssID[0]
                    #print(mssID, mss_df)
                    
                    Qmedian, Qmean, Qstd, Q_lb, Q_ub = get_stats(mss_df['Q'].values)
                    mssID_Qstats = (Qmean, Qstd, Q_lb, Q_ub)

                    Gmedian, Gmean, Gstd, G_lb, G_ub = get_stats(mss_df['G'].values)
                    mssID_Gstats = (Gmean, Gstd, G_lb, G_ub)

                    Kmedian, Kmean, Kstd, K_lb, K_ub = get_stats(mss_df['K'].values)
                    mssID_Kstats = (Kmean, Kstd, K_lb, K_ub)

                    #print(f'<Q> = {Qmean} | std(Q) = {Qstd} | 95%ci = ({Q_lb}, {Q_ub})')
                    #print(f'<G> = {Gmean} | std(G) = {Gstd} | 95%ci = ({G_lb}, {G_ub})')
                    #print(f'<K> = {Kmean} | std(K) = {Kstd} | 95%ci = ({K_lb}, {K_ub})')

                    ## update dictionary traccking the mean G and Q for this protein only
                    local_mss_GQmeans['mssID'] += [mssID]
                    local_mss_GQmeans['Qmean'] += [Qmean]
                    local_mss_GQmeans['Gmean'] += [Gmean]

                    ## update the global self.data dictionary MSM level with the Q and G stats
                    self.data[tag]['MSM']['Qstats'][mssID] = mssID_Qstats
                    self.data[tag]['MSM']['Gstats'][mssID] = mssID_Gstats
                    self.data[tag]['MSM']['Kstats'][mssID] = mssID_Kstats

                # determine the native MSM state by the highest <Q> and lowest <G>
                logging.info(f'Determine the native MSM state by the highest <Q> and lowest <G>')
                local_mss_GQmeans = pd.DataFrame(local_mss_GQmeans)
                local_mss_GQmeans = local_mss_GQmeans.sort_values(by=['Qmean', 'Gmean'], ascending=False, ignore_index=True)
                #print(local_mss_GQmeans)
                nativeMSS = local_mss_GQmeans.iloc[0]['mssID']
                logging.info(f'{tag} | nativeMSS: {nativeMSS}')
                df['NativeByMSS'] = df['metastablestate'] == nativeMSS

                # determine the native frames by having a Q and G within 3sigma of the reference sim averages
                logging.info(f'Determine the native frames by having a Q and G within 3sigma of the reference sim averages')
                Native_Qmean = self.data[tag]['Qstats']['Native'][0] 
                Native_Qstd = self.data[tag]['Qstats']['Native'][1]
                Qthreshold = Native_Qmean - 3*Native_Qstd
                logging.info(f'Native_Qmean: {Native_Qmean} | Native_Qstd: {Native_Qstd} | Qthreshold: {Qthreshold}')

                Native_Gmean = self.data[tag]['Gstats']['Native'][0] 
                Native_Gstd = self.data[tag]['Gstats']['Native'][1]
                Gthreshold = Native_Gmean + 3*Native_Gstd
                logging.info(f'Native_Gmean: {Native_Gmean} | Native_Gstd: {Native_Gstd} | Gthreshold: {Gthreshold}')

                df['NativeByRef'] = (df['Q'] >= Qthreshold) & (df['G'] <= Gthreshold)


                ###################################################################################################
                ## loop through again now that we have assigned native frames and MSM stateas and make a summary stats file for output
                for mssID, mss_df in df.groupby(['metastablestate']):
                    mssID = mssID[0]
                    #print(mss_df)

                    mssID_Qstats = self.data[tag]['MSM']['Qstats'][mssID]
                    Qmean, Qstd, Q_lb, Q_ub = mssID_Qstats

                    mssID_Gstats = self.data[tag]['MSM']['Gstats'][mssID]
                    Gmean, Gstd, G_lb, G_ub = mssID_Gstats

                    mssID_Kstats = self.data[tag]['MSM']['Kstats'][mssID]
                    Kmean, Kstd, K_lb, K_ub = mssID_Kstats

                    ## update the stats dataframe containing just high level stats info for each set of Q and G values analuzed
                    self.stats_df['gene_tag'] += [tag]
                    self.stats_df['macrostateID'] += [mssID]
                    self.stats_df['desc'] += ['T-Quench sims']
                    self.stats_df['Qmean'] += [Qmean]
                    self.stats_df['Qstd'] += [Qstd]
                    self.stats_df['Q_lb'] += [Q_lb]
                    self.stats_df['Q_ub'] += [Q_ub]
                    self.stats_df['Gmean'] += [Gmean]
                    self.stats_df['Gstd'] += [Gstd]
                    self.stats_df['G_lb'] += [G_lb]
                    self.stats_df['G_ub'] += [G_ub]
                    self.stats_df['Kmean'] += [Kmean]
                    self.stats_df['Kstd'] += [Kstd]
                    self.stats_df['K_lb'] += [K_lb]
                    self.stats_df['K_ub'] += [K_ub]
                    self.stats_df['n'] += [len(mss_df)]
                    if mssID == nativeMSS:
                        self.stats_df['native'] += [True]
                    else:
                         self.stats_df['native'] += [False]
                

        self.stats_df = pd.DataFrame(self.stats_df)
        logging.info(f'self.stats_df:\n{self.stats_df}')
        stats_outfile = os.path.join(self.data_path, f'GQstats.csv')
        self.stats_df.to_csv(stats_outfile, index=False)
        logging.info(f'SAVED: {stats_outfile}')
        print(self.stats_df)

        #################################################################################
        ## Save the GQK and NativeBy dataframe dataframes 
        mss_dfs = []
        for tag, tag_data in self.data.items():
            gene, pdb, chain = tag.split('_')
            df = tag_data['df']
            df['gene'] = gene
            df['pdb'] = pdb
            df['chain'] = chain
            #print(df)
            mss_dfs += [df]
        mss_dfs = pd.concat(mss_dfs)
        print(f'mss_dfs:\n{mss_dfs}')
        mss_df_outfile = os.path.join(self.data_path, f'CollectedGQK.csv')
        mss_dfs.to_csv(mss_df_outfile)
        print(f'SAVED: {mss_df_outfile}')

        EntInfo_dfs = []
        for i, (tag, tag_data) in enumerate(self.data.items()):
            gene, pdb, chain = tag.split('_')
            df = tag_data['EntInfo']
            df['gene'] = gene
            df['pdb'] = pdb
            df['chain'] = chain
            #print(df)
            EntInfo_dfs_outfile = os.path.join(self.data_path, f'CollectedEntInfo_{tag}.csv')
            df.to_csv(EntInfo_dfs_outfile)
            print(f'SAVED: {EntInfo_dfs_outfile} {i}')
        ###################################################################################################
    
    #######################################################################################

#######################################################################################
def bootstrap(data):
    boot_means = []
    for b in range(10000):
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

    # Step 1: load the G and Q data
    anal.load_OP()

    # Step 2: for each candidate get the native G and Q states
    anal.get_nativeSIM_stats()
    
    # Step 3: Identify the native state in the MSM data
    anal.find_native_state()
    
    print(f'logfile: {logfile}')

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()

print(f'NORMAL TERMINATION: {time.time() - start_time}')
logging.info(f'NORMAL TERMINATION: {time.time() - start_time}')