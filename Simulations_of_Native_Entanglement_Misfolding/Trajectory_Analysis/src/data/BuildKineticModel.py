#!/usr/bin/env python3
import logging, os, sys
import time
import argparse
import pandas as pd
import numpy as np
import glob
import pyemma as pem
import parmed as pmd
import mdtraj as mdt
import matplotlib.pyplot as plt
import msmtools
from matplotlib.cm import get_cmap
from matplotlib.colors import ListedColormap, BoundaryNorm

#pd.set_option('display.max_rows', 5000)

class Analysis:
    """
    A class to build higherachical clustering models across a set of trajectories
    """
    #######################################################################################
    def __init__(self, args):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
        ("--outpath", type=str, required=True, help="Path to output directory")
        ("--OPpath", type=str, required=True, help="Path to directory containing G and Q directories created by GQ.py")
        ("--outname", type=str, required=True, help="base name for output files")
        ("--psf", type=str, required=True, help="Path to CA protein structure file")
        ("--dcds", type=str, required=True, help="Path to trajectory to analyze")
        ("--start", type=int, required=False, help="First frame to analyze 0 indexed", default=0)
        ("--end", type=int, required=False, help="Last frame to analyze 0 indexed", default=-1)
        ("--stride", type=int, required=False, help="Frame stride", default=1)
        """

        # parse the parameters 
        self.OPpath = args.OPpath
        logging.info(f'OPpath: {self.OPpath}')

        self.outpath = args.outpath
        logging.info(f'outpath: {self.outpath}')

        self.outname = args.outname
        logging.info(f'outname: {self.outname}')

        self.psf = args.psf
        logging.info(f' psf: {self. psf}')

        self.dcds = args.dcds 
        logging.info(f'dcds: {self.dcds}')

        self.start = args.start
        self.end = args.end
        self.stride = args.stride
        print(f'START: {self.start} | END: {self.end} | STRIDE: {self.stride}')

        self.n_cluster = 50 # Number of k-means clusters to group. Default is 400.
        self.kmean_stride = 10 # Stride of reading trajectory frame when clustring by k-means. 
        self.n_small_states = 2 # Number of clusters for the inactive microstates after MSM modeling to be clustered into
        self.n_large_states = 10  # Adjust based on your system
        self.dt = 0.015/1000 # timestep used in MD simulations in ns

    #######################################################################################

    #######################################################################################
    def load_OP(self,):
        """
        Loads the GQ values of each trajectory into a 2D array and then appends it to a list
        The list should have Nt = number of trajectories and each array should be n x 2 where n is the number of frames
        """

        cor_list = []
        cor_list_idx_2_traj = {}
        Qfiles = glob.glob(os.path.join(self.OPpath, 'Q/*.Q'))
        Gfiles = glob.glob(os.path.join(self.OPpath, 'G/*.G'))
        print(f'Number of Q files found: {len(Qfiles)} | Number of G files found: {len(Gfiles)}')

        # loop through the Qfiles and find matching Gfile
        # then load the Q and G time series into a 2D array
        idx = 0
        for Qf in Qfiles:
            
            traj = Qf.split('_')[-1].replace('.Q','').replace('t','')
            cor_list_idx_2_traj[idx] = int(traj)
            idx += 1

            # get the cooresponding G file
            Gf = [f for f in Gfiles if f't{traj}.G' in f]
            # quality check to ensure that only 1 Gfile was found for each Qfile
            if len(Gf) != 1:
                raise ValueError(f'Zero or more than 1 Gfile found for traj {traj}')
            else:
                Gf = Gf[0]

            # load the G Q data and extract only the time series column
            Qdata = pd.read_csv(Qf)['Q'].values[self.start:self.end:self.stride]
            Gdata = pd.read_csv(Gf)['G'].values[self.start:self.end:self.stride]

            # Quality check that both time series have the same length
            if Qdata.shape != Gdata.shape:
                raise ValueError(f'Q and G do not have the same shape! {Qdata.shape} {Gdata.shape}')
            
            data = np.stack((Qdata, Gdata)).T
            print(data[:3, :])
            cor_list.append(data)

        print(f'Number of trajecotry OP coordinate loaded: {len(cor_list)}')
        self.cor_list = cor_list

        print(f'Mapping of cor_list index to trajID in file names: {cor_list_idx_2_traj}')
        self.cor_list_idx_2_traj = cor_list_idx_2_traj
    #######################################################################################  

    #######################################################################################
    def standardize(self,):
        """
        Standardizes your OP by taking the mean and std Q and G across all traj data and rescaling each trajectorys data by
        Z = (d - mean)/std 
        """
        data_con = self.cor_list[0]
        for i in range(1, len(self.cor_list)):
            data_con = np.vstack((data_con, self.cor_list[i]))
        self.data_mean = np.mean(data_con, axis=0)
        self.data_std = np.std(data_con, axis=0)
        self.standard_cor_list = [(d - self.data_mean) / self.data_std for d in self.cor_list]
    #######################################################################################

    #######################################################################################
    def unstandardize(self, data):
        """
        Unstandardizes your OP by taking the mean and std Q and G across all traj data and rescaling each trajectorys data by
        Z*std + mean = d
        """
        return data*self.data_std + self.data_mean
    #######################################################################################

    #######################################################################################
    def cluster(self,):
        """
        Cluster the GQ data across all trajectories using kmeans. 
        dtrajs contains the resulting kmeans cluster labels for each trajectory time series 
        centers contains the standardized GQ coordinates of the cluster centers
        """
        self.clusters = pem.coordinates.cluster_kmeans(self.standard_cor_list, k=self.n_cluster, max_iter=5000, stride=self.kmean_stride)
        #print(cluster)
        dtrajs = self.clusters.dtrajs
        print(f'dtrajs: {len(dtrajs)} {dtrajs[0].shape}\n{dtrajs[0][:10]}')
        standard_centers = self.clusters.clustercenters
        unstandard_centers = self.unstandardize(standard_centers)
        print(f'unstandard_centers:\n{unstandard_centers} {unstandard_centers.shape}')
    #######################################################################################

    #######################################################################################
    def build_msm(self, lagtime=1):
        print(f'Building MSM model with a lag time of {lagtime}')

        # Get count matrix and connective groups of microstates
        c_matrix = msmtools.estimation.count_matrix(self.clusters.dtrajs, lagtime).toarray()
        print(f'c_matrix:\n{c_matrix} {c_matrix.shape}')
        sub_groups = msmtools.estimation.connected_sets(c_matrix)
        print(f'Total number of sub_groups: {len(sub_groups)}\n{sub_groups}')
        
        # Build the MSM models for any connected sets that have more than 1 microstate
        msm_list = []        
        for sg in sub_groups:
            cm = msmtools.estimation.largest_connected_submatrix(c_matrix, lcc=sg)
            print(f'For sub_group: {sg}')
            if len(cm) == 1:
                msm = None
            else:
                print(f'Building Transition matrix and MSM model')
                T = msmtools.estimation.transition_matrix(cm, reversible=True)
                msm = pem.msm.markov_model(T, dt_model=str(self.dt)+' ns')
            msm_list.append(msm)
        logging.info(f'Number of models: {len(msm_list)}')

        # Coarse grain out the metastable macrostates in the models
        print(f'Coarse grain out the metastable macrostates in the models')
        logging.info(f'Coarse grain out the metastable macrostates in the models')
        meta_dist = []
        meta_set = []
        eigenvalues_list = []
        for idx_msm, msm in enumerate(msm_list):

            # the first model should contain the largest connected state so use the largest number of metastable states
            # for every other subgroup use the smallest
            if idx_msm == 0:
                n_states = self.n_large_states
            else:
                n_states = self.n_small_states

            if msm == None:
                eigenvalues_list.append(None)
                dist = np.zeros(self.n_cluster)
                iidx = sub_groups[idx_msm][0]
                dist[iidx] = 1.0
                meta_dist.append(dist)
                meta_set.append(sub_groups[idx_msm])

            else:
                eigenvalues_list.append(msm.eigenvalues())
                # coarse-graining 
                while n_states > 1:
                    tag_empty = False
                    pcca = msm.pcca(n_states)
                    for ms in msm.metastable_sets:
                        if ms.size == 0:
                            tag_empty = True
                            break
                    if not tag_empty:
                        break
                    else:
                        n_states -= 1
                        print('Reduced number of states to %d for active group %d'%(n_states, idx_msm+1))
                        logging.info('Reduced number of states to %d for active group %d'%(n_states, idx_msm+1))
                if n_states == 1:
                    # use observation prob distribution for non-active set
                    dist = np.zeros(self.n_cluster)
                    for nas in sub_groups[idx_msm]:
                        for dtraj in dtrajs:
                            dist[nas] += np.count_nonzero(dtraj == nas)
                    dist /= np.sum(dist)
                    meta_dist.append(dist)
                    meta_set.append(sub_groups[idx_msm])
                else:
                    for i, md in enumerate(msm.metastable_distributions):
                        dist = np.zeros(self.n_cluster)
                        s = np.sum(md[msm.metastable_sets[i]])
                        set_0 = []
                        for idx in msm.metastable_sets[i]:
                            iidx = sub_groups[idx_msm][idx]
                            dist[iidx] = md[idx]
                            set_0.append(iidx)
                        dist = dist / s
                        meta_dist.append(dist)
                        meta_set.append(set_0)
        meta_dist = np.array(meta_dist)
        meta_set = np.array(meta_set)
   
        ## make microstate to metastable state mapping object
        logging.info(f'\nMetastable state assignment')
        meta_mapping = {}
        for metaID, microstates in enumerate(meta_set):
            logging.info(metaID, microstates)
            for m in microstates:
                if m not in meta_mapping:
                    meta_mapping[m] = metaID
                else:
                    raise ValueError(f'Microstate {m} already in a metastable state!')
        print(f'meta_mapping: {meta_mapping} {len(meta_mapping)}')

        # map those microstate states to the metastable state
        metastable_dtraj = []
        for dtraj_idx, dtraj in enumerate(self.clusters.dtrajs):
            mapped_dtraj = []
            for d in dtraj:
                mapped_dtraj.append(meta_mapping[d])

            #rint(mapped_dtraj)
            metastable_dtraj += [np.asarray(mapped_dtraj)]

        print(f'Metastable state mapping:')
        for dtraj_idx, dtraj in enumerate(metastable_dtraj):
            print(dtraj_idx, self.clusters.dtrajs[dtraj_idx][:10] , dtraj[:10])

        ## get samples of metastable states by most populated microstates
        cluster_indexes = msmtools.dtraj.index_states(self.clusters.dtrajs)
        if len(cluster_indexes) < self.n_cluster:
            cluster_indexes = list(cluster_indexes)
            for i in range(len(cluster_indexes), self.n_cluster):
                cluster_indexes.append(np.array([[]]))
            cluster_indexes = np.array(cluster_indexes)
        samples = msmtools.dtraj.sample_indexes_by_distribution(cluster_indexes, meta_dist, 5)
        print(f'samples: {samples} {len(samples)}')
        
        ## Make the output dataframe that has assignments for each frame of each traj
        df = {'traj':[], 'frame':[], 'microstate':[], 'metastablestate':[], 'Q':[], 'G':[], 'StateSample':[]}
        print(f'Active & inactive metastable state mapping')
        for k,v in enumerate(metastable_dtraj):
            traj = self.cor_list_idx_2_traj[k]
            print(k, traj, v[:10])
            for frame, macrostate in enumerate(v):
                microstate = self.clusters.dtrajs[k][frame]
                Q = self.cor_list[k][frame, 0]
                G = self.cor_list[k][frame, 1]

                if [k, frame] in samples[macrostate].tolist():
                    StateSample = True
                else:
                    StateSample = False
                #print(k, frame, microstate, macrostate, StateSample)
                df['traj'] += [traj]
                df['frame'] += [frame]
                df['microstate'] += [microstate]
                df['metastablestate'] += [macrostate]
                df['Q'] += [Q]
                df['G'] += [G]
                df['StateSample'] += [StateSample]

        df = pd.DataFrame(df)
        print(f'Final MSM mapping DF:\n{df}')
        df_outfile = os.path.join(self.outpath, f'{self.outname}_MSMmapping.csv')
        df.to_csv(df_outfile, index=False)
        logging.info(f'SAVED: {df_outfile}')

        # Plot the metastable state membership and free energy surface
        xall = np.hstack([dtraj[:, 0] for dtraj in self.cor_list])
        yall = np.hstack([dtraj[:, 1] for dtraj in self.cor_list])
        states = np.hstack(metastable_dtraj)
        print(f'xall: {xall} {xall.shape}')
        print(f'yall: {yall} {yall.shape}')
        print(f'states: {states} {states.shape}')

        stateplot_outfile = os.path.join(self.outpath, f'{self.outname}_StateAndFEplot.png')
        self.plot_state_map_and_FE(xall, yall, states, stateplot_outfile)
    #######################################################################################

    #######################################################################################
    def plot_state_map_and_FE(self, x, y, states, outfile, cmap='viridis', point_size=50, alpha=0.85, title='State Map'):
        """
        Plots a state map using x and y values colored by state assignments with labeled colorbar.
        Parameters:
            x (array-like): The x-coordinates of the points.
            y (array-like): The y-coordinates of the points.
            states (array-like): The state assignment for each point.
            cmap (str or Colormap): Colormap for state coloring (default is 'viridis').
            point_size (int): Size of the scatter plot points (default is 50).
            alpha (float): Transparency of the points (default is 0.7).
            title (str): Title of the plot (default is 'State Map with Labels').
        """

        # Create a figure and subplots with 1 row and 2 columns
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # Ensure states are integers
        states = np.asarray(states, dtype=int)
        unique_states = np.unique(states)
        n_states = len(unique_states)

        # Create a colormap with distinct colors for each state
        if isinstance(cmap, str):
            cmap = plt.get_cmap(cmap, n_states)
        else:
            cmap = ListedColormap(cmap.colors[:n_states])

        # Create a normalized colormap with discrete boundaries
        norm = BoundaryNorm(boundaries=np.arange(n_states+1)-0.5, ncolors=n_states)

        # Plot the state map
        scatter = axes[0].scatter(x, y, c=states, cmap=cmap, s=point_size, edgecolor='k', alpha=alpha, norm=norm)

        # Create a colorbar with state labels
        cbar = fig.colorbar(scatter, ax=axes[0], ticks=np.arange(n_states), boundaries=np.arange(n_states + 1) - 0.5)
        cbar.set_label('States')
        cbar.set_ticks(unique_states)
        cbar.set_ticklabels([f'{state}' for state in unique_states])

        axes[0].set_xlabel('Q')
        axes[0].set_ylabel('G')
        axes[0].set_title(title)
        axes[0].grid(True)

        _, axes[1] = pem.plots.plot_free_energy(x, y, cmap='viridis')
        axes[1].set_title('Population map')
        plt.tight_layout()
        plt.savefig(outfile)
        logging.info(f'SAVED: {outfile}')
    #######################################################################################

    #######################################################################################
    def plot_implied_timescales(self,):
        """
        Should be done first before building the model to find a proper lagtime for which the timescales (eignenvalues of the transition matrix) of the model are no longer dependant.
        Look for the point or range where the implied timescales stop changing significantly with increasing lag times. 
        This lag time is generally a good choice for building your MSM, as it suggests the dynamics are being captured without undue dependence on the initial conditions.
        """
        #nits = -1
        lag_times = np.arange(1, 100, 10)  # adjust the range based on your system
        n_states = len(np.unique(self.clusters.dtrajs))  # or a predefined number of states
        its = pem.msm.its(self.clusters.dtrajs, lags=lag_times, errors='bayes')
        pem.plots.plot_implied_timescales(its)
        ITS_outfile = os.path.join(self.outpath, f'{self.outname}_ITS.png')
        plt.savefig(ITS_outfile)
        logging.info(f'SAVED: {ITS_outfile}')
    #######################################################################################

    #######################################################################################
############## MAIN #################
def main():
    
    script_name = f'BuildKineticModel'
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("--OPpath", type=str, required=True, help="Path to directory containing G and Q directories created by GQ.py")
    parser.add_argument("--outname", type=str, required=True, help="base name for output files")
    parser.add_argument("--psf", type=str, required=True, help="Path to CA protein structure file")
    parser.add_argument("--dcds", type=str, required=True, help="Path to trajectories created by ")
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

    # load the G and Q data
    anal.load_OP()

    # apply the standard scalar transformation to the data 
    anal.standardize()

    # cluster the standardized data using kmeans clustering with a stride of 10, change this if necessary
    anal.cluster()

    # genereate the implied timescales plot to check for a suitable lag time
    # Should be done first before building the model to choose a suitable lag time
    #anal.plot_implied_timescales()
    
    # Build the MSM model with the choosen lagtime
    anal.build_msm(lagtime=40)

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    
print(f'NORMAL TERMINATION: {time.time() - start_time}')
logging.info(f'NORMAL TERMINATION: {time.time() - start_time}')