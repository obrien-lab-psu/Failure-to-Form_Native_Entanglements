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
    A class to collect the G, Q, K and calcualte the running average for each trajectory and plot the results
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


        ## get the set of G, Q, and K files
        self.QFiles = glob.glob(os.path.join(self.toplevel, '*/Q/*.Q'))
        print(f'Number of Q files found: {len(self.QFiles)}')

        self.GFiles = glob.glob(os.path.join(self.toplevel, '*/Cluster_ChangesInEnt/*_clustered.G'))
        print(f'Number of G files found: {len(self.GFiles)}')

        self.KFiles = glob.glob(os.path.join(self.toplevel, '*/Mirror/Quench/*.dat'))
        print(f'Number of K files found: {len(self.KFiles)}')


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
        self.candidates = self.candidates[self.candidates['set'] == self.setID] ## only get those candidates in this set
        print(self.candidates)
        for gene, pdb, chain in self.candidates[['gene', 'pdb', 'chain']].values:
            tag = f'{gene}_{pdb}_{chain}'
            print(f'Loading files for {gene} {pdb} {chain} {tag}')
            file = f'../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CollectAndProcessOP/setID{self.setID}/DATA/{tag}_Quench_Collected_GQK.csv'
            if os.path.exists(file):
                print(f'Loading {file}')
                df = pd.read_csv(file)
                Quench_dfs += [df]
            else:
                raise ValueError(f'File {file} does not exist')
        Quench_dfs = pd.concat(Quench_dfs, ignore_index=True)
        print(f'Quench_dfs:\n{Quench_dfs}')

        return Quench_dfs
    #######################################################################################  

    #######################################################################################
    def RunningAverage(self, df):
        print(f'Test of convergence by plotting Running averages')
        """
        For each order parameter Q, G, K, Z calculate the running average and standard deviation
        """
        outfile = os.path.join(self.data_path, f'RunningAverage.csv')
        if not os.path.exists(outfile):
            outdf = {'gene':[], 'traj':[], 'frame':[], 'RunningAvg_Q':[], 'RunningStd_Q':[], 'RunningAvg_G':[], 'RunningStd_G':[], 'RunningAvg_K':[], 'RunningStd_K':[], 'RunningAvg_Z':[], 'RunningStd_Z':[]}
            for gene, gene_df in df.groupby('gene'):
                for traj, traj_df in gene_df.groupby('traj'):
                    print(traj_df)
                    print(gene, traj)
         
                    Qdata = traj_df['Q'].values
                    Gdata = traj_df['G'].values
                    Kdata = traj_df['K'].values
                    Zdata = traj_df['Z'].values
                    
                    # for each frame calculate the running average at time t' which is the normalized sum of the values from t=0 to t=t'
                    for frame, value in enumerate(Qdata):
                        running_avg_Q = np.mean(Qdata[:frame+1])
                        running_std_Q = np.std(Qdata[:frame+1])

                        running_avg_G = np.mean(Gdata[:frame+1])
                        running_std_G = np.std(Gdata[:frame+1])

                        running_avg_K = np.mean(Kdata[:frame+1])
                        running_std_K = np.std(Kdata[:frame+1])

                        running_avg_Z = np.mean(Zdata[:frame+1])
                        running_std_Z = np.std(Zdata[:frame+1])

                        outdf['gene'] += [gene]
                        outdf['traj'] += [traj]
                        outdf['frame'] += [frame]

                        outdf['RunningAvg_Q'] += [running_avg_Q]
                        outdf['RunningStd_Q'] += [running_std_Q]

                        outdf['RunningAvg_G'] += [running_avg_G]    
                        outdf['RunningStd_G'] += [running_std_G]

                        outdf['RunningAvg_K'] += [running_avg_K]
                        outdf['RunningStd_K'] += [running_std_K]

                        outdf['RunningAvg_Z'] += [running_avg_Z]
                        outdf['RunningStd_Z'] += [running_std_Z]

                        
            outdf = pd.DataFrame(outdf)
            print(f'outdf:\n{outdf}')
            outdf.to_csv(outfile, index=False)
            print(f'SAVED: {outfile}')
        else:
            outdf = pd.read_csv(outfile)
            print(f'LOADED: {outfile}') 
            print(f'outdf:\n{outdf}')
        return outdf
    #######################################################################################

    #######################################################################################
    def plot_convergence_heatmaps(self, df, output_dir='data/heatmap_outputs'):
        """
        Generates and saves vertical heatmaps (1 column x 8 rows) of running averages and standard deviations
        for each order parameter (Q, G, K, Z) per gene from the given dataframe. Adds a vertical dashed red
        line at frame 24000 and sets y-axis ticks every 5 trajectories.

        Parameters:
        - df (pd.DataFrame): Input dataframe with structure as described.
        - output_dir (str): Directory to save the output figures.
        """
        import matplotlib.ticker as ticker

        os.makedirs(output_dir, exist_ok=True)
        order_params = ['Q', 'G', 'K', 'Z']
        line_position = 24000

        for gene, gene_df in df.groupby('gene'):
            fig, axes = plt.subplots(nrows=8, ncols=1, figsize=(12, 32), sharex=True, constrained_layout=True)
            fig.suptitle(f'Running Average and Std Dev Heatmaps for {gene}', fontsize=20)

            for i, param in enumerate(order_params):
                for j, stat in enumerate(['Avg', 'Std']):
                    ax = axes[i * 2 + j]
                    value_col = f'Running{stat}_{param}'
                    pivot_table = gene_df.pivot(index='traj', columns='frame', values=value_col)

                    cbar_label = f"<{param}>" if stat == 'Avg' else f"Std({param})"
                    sns.heatmap(
                        pivot_table,
                        ax=ax,
                        cmap='nipy_spectral',
                        cbar=True,
                        cbar_kws={'label': cbar_label}
                    )

                    ax.set_title(f'{param} - Running {stat}', fontsize=12)
                    ax.set_ylabel('Trajectory')

                    # Set x-axis ticks every 500 frames
                    max_frame = pivot_table.columns.max()
                    ax.set_xticks(range(0, max_frame + 1, 500))
                    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())

                    # Set y-axis ticks every 5 trajectories
                    max_traj = pivot_table.index.max()
                    ax.set_yticks(range(0, max_traj + 1, 5))
                    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())

                    # Draw vertical red dashed line at frame 24000
                    ax.axvline(x=line_position, color='white', linestyle='--', linewidth=1)

            axes[-1].set_xlabel('Frame')

            #output_path = os.path.join(self.data_path, output_dir)
            output_path = os.path.join(output_dir, f'{gene}_convergence_heatmaps.png')
            plt.savefig(output_path, dpi=300)
            plt.close(fig)
            print(f"Saved vertical heatmaps for {gene} in '{output_path}'.")
       

        print(f"Saved vertical heatmaps for {df['gene'].nunique()} proteins in '{output_dir}'.")

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
    
    # Step 1: load the G and Q data
    Quench_df = anal.load_OP()

    # Step 2: for each candidate get the  G, Q, K running average
    RunQuench_df = anal.RunningAverage(Quench_df)


    ## make logging dir
    heatmap_outdir = f'data/heatmap_outputs/setID{args.setID}/'
    if not os.path.exists(heatmap_outdir):
        os.makedirs(heatmap_outdir)
        print(f'Made directory: {heatmap_outdir}')  

    anal.plot_convergence_heatmaps(RunQuench_df, output_dir=heatmap_outdir)

    print(f'logfile: {logfile}')

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()

print(f'NORMAL TERMINATION: {time.time() - start_time}')
logging.info(f'NORMAL TERMINATION: {time.time() - start_time}')