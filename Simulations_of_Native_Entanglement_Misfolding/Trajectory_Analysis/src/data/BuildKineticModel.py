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
#import msmtools
import deeptime
from matplotlib.cm import get_cmap
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.colors as mcolors
import seaborn as sns
from scipy.stats import mode
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


        #self.dcds = args.dcds 
        #logging.info(f'dcds: {self.dcds}')

        self.start = args.start
        self.end = args.end
        self.stride = args.stride
        print(f'START: {self.start} | END: {self.end} | STRIDE: {self.stride}')

        self.n_cluster = 250 # Number of k-means clusters to group. Default is 400.
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
        Gfiles = glob.glob(os.path.join(self.OPpath, 'Cluster_ChangesInEnt/*_clustered.G'))
        print(f'Number of Q files found: {len(Qfiles)} | Number of G files found: {len(Gfiles)}')

        # loop through the Qfiles and find matching Gfile
        # then load the Q and G time series into a 2D array
        idx = 0
        for Qf in Qfiles:
            print(f'Qf: {Qf}')

            traj = Qf.split('_')[-1].replace('.Q','').replace('t','')
            print(f'Traj: {traj}')
            cor_list_idx_2_traj[idx] = int(traj)
            idx += 1

            # get the cooresponding G file
            Gf = [f for f in Gfiles if f't{traj}_clustered.G' in f]
            # quality check to ensure that only 1 Gfile was found for each Qfile
            if len(Gf) != 1:
                raise ValueError(f'Zero or more than 1 Gfile found for traj {traj}')
            else:
                Gf = Gf[0]
            print(f'Gf: {Gf}')

            # load the G Q data and extract only the time series column
            #Qdata = pd.read_csv(Qf)['Q'].values[self.start:self.end:self.stride]
            Qdata = pd.read_csv(Qf)
            Qdata = Qdata[(Qdata['Frame'] >= self.start) & (Qdata['Frame'] <= self.end)]
            Qdata = Qdata['Q'].values

            Gdata = pd.read_csv(Gf)
            Gdata = Gdata[(Gdata['Frame'] >= self.start) & (Gdata['Frame'] <= self.end)]
            Gdata = Gdata['G'].values 
            print(f'Shape of OP: Q {Qdata.shape} G {Gdata.shape}')

            # Quality check that both time series have the same length
            if Qdata.shape != Gdata.shape:
                raise ValueError(f'Q and G do not have the same shape! {Qdata.shape} {Gdata.shape}')
            
            data = np.stack((Qdata, Gdata)).T
            print(data[:3, :], data.shape)
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

        if the number of unique centers found is not equal to self.n_cluster then adjust it to reflect the number found. This can happen if you have data that has a narrow distribution. 
        """
        self.clusters = pem.coordinates.cluster_kmeans(self.standard_cor_list, k=self.n_cluster, max_iter=5000, stride=self.kmean_stride)
        
        # Get the microstate tagged trajectories and their state counts
        self.dtrajs = self.clusters.dtrajs
        print(f'dtrajs: {len(self.dtrajs)} {self.dtrajs[0].shape}\n{self.dtrajs[0][:10]}')
        clusterIDs, counts = np.unique(self.dtrajs, return_counts=True)
        print(f'Number of unique microstate IDs: {len(clusterIDs)} {clusterIDs}')
        
        state_counts = {}
        for i,c in zip(clusterIDs, counts):
            state_counts[i] = c
        print(f'state_counts: {state_counts}')
        
        # Quality check that all microstate ids are assigned
        # If not renumber from 0
        if len(clusterIDs) != self.n_cluster:
            print(f'The number of microstate IDs assigned does not match the number specified: {len(clusterIDs)} != {self.n_cluster}')

            mapping_dict = {}
            for new,old in enumerate(clusterIDs):
                mapping_dict[old] = new
            print(f'mapping_dict: {mapping_dict}')

            # Convert the dictionary to a numpy array for efficient mapping
            max_key = max(mapping_dict.keys())
            mapping_array = np.zeros(max_key + 1, dtype=int)
            for key, value in mapping_dict.items():
                mapping_array[key] = value

            # Map the arrays using the mapping array
            self.dtrajs = [mapping_array[arr] for arr in self.dtrajs]
            
            clusterIDs, counts = np.unique(self.dtrajs, return_counts=True)
            print(f'Number of unique microstate IDs after mapping: {len(clusterIDs)} {clusterIDs}')
            state_counts = {}
            for i,c in zip(clusterIDs, counts):
                state_counts[i] = c
            print(f'state_counts: {state_counts}')

            self.n_cluster = len(clusterIDs)
        

        standard_centers = self.clusters.clustercenters
        unstandard_centers = self.unstandardize(standard_centers)
        print(f'unstandard_centers:\n{unstandard_centers} {unstandard_centers.shape}')
        print(f'self.n_cluster: {self.n_cluster}')
        
    #######################################################################################

    #######################################################################################
    def build_msm(self, lagtime=1):
        print(f'Building MSM model with a lag time of {lagtime}')

        # Get count matrix and connective groups of microstates
        c_matrix = deeptime.markov.tools.estimation.count_matrix(self.dtrajs, lagtime).toarray()
        print(f'c_matrix:\n{c_matrix} {c_matrix.shape}')
        
        sub_groups = deeptime.markov.tools.estimation.connected_sets(c_matrix)
        print(f'Total number of sub_groups: {len(sub_groups)}\n{sub_groups}')
        
        # Build the MSM models for any connected sets that have more than 1 microstate
        msm_list = []        
        for sg in sub_groups:
            cm = deeptime.markov.tools.estimation.largest_connected_submatrix(c_matrix, lcc=sg)
            print(f'For sub_group: {sg}')
            if len(cm) == 1:
                msm = None
            else:
                print(f'Building Transition matrix and MSM model')
                T = deeptime.markov.tools.estimation.transition_matrix(cm, reversible=True)
                msm = pem.msm.markov_model(T, dt_model=str(self.dt)+' ns')
            msm_list.append(msm)
        logging.info(f'Number of models: {len(msm_list)}')
        print(f'Number of models: {len(msm_list)}')

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
        print(f'meta_set: {meta_set}')
        #meta_set = np.array(meta_set)
   
        ## make microstate to metastable state mapping object
        print(f'\nMetastable state assignment')
        meta_mapping = {}
        for metaID, microstates in enumerate(meta_set):
            print(metaID, microstates)
            for m in microstates:
                if m not in meta_mapping:
                    meta_mapping[m] = metaID
                else:
                    raise ValueError(f'Microstate {m} already in a metastable state!')
        print(f'meta_mapping: {meta_mapping} {len(meta_mapping)}')

        # map those microstate states to the metastable state
        metastable_dtraj = []
        for dtraj_idx, dtraj in enumerate(self.dtrajs):
            mapped_dtraj = []
            for d in dtraj:
                mapped_dtraj.append(meta_mapping[d])

            #rint(mapped_dtraj)
            metastable_dtraj += [np.asarray(mapped_dtraj)]

        print(f'Metastable state mapping:')
        for dtraj_idx, dtraj in enumerate(metastable_dtraj):
            print(dtraj_idx, self.dtrajs[dtraj_idx][:10] , dtraj[:10], dtraj.shape)
      

        ## get samples of metastable states by most populated microstates
        print(len(self.dtrajs), self.dtrajs[0].shape)
        cluster_indexes = deeptime.markov.sample.compute_index_states(self.dtrajs)
        print(f'cluster_indexes: {len(cluster_indexes)}')
        print(f'meta_dist: {len(meta_dist)} {meta_dist.shape}')

        samples = deeptime.markov.sample.indices_by_distribution(cluster_indexes, meta_dist, 5)
        print(f'samples: {samples} {len(samples)}')
        
        ## Make the output dataframe that has assignments for each frame of each traj
        df = {'traj':[], 'frame':[], 'microstate':[], 'metastablestate':[], 'Q':[], 'G':[], 'StateSample':[]}
        print(f'Active & inactive metastable state mapping')
        for k,v in enumerate(metastable_dtraj):
            traj = self.cor_list_idx_2_traj[k]
            print(k, traj, v[:10])
            for frame, macrostate in enumerate(v):
                microstate = self.dtrajs[k][frame]
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
        df['frame'] += self.start # correct the frame index to start from the start specified by the user as this frame index starts from 0
        print(f'Final MSM mapping DF:\n{df}')
        df_outfile = os.path.join(self.outpath, f'{self.outname}_MSMmapping.csv')
        df.to_csv(df_outfile, index=False)
        logging.info(f'SAVED: {df_outfile}')
        print(f'SAVED: {df_outfile}')

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
        #############################################################################################
        # Create a figure and subplots with 1 row and 2 columns
        fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

        ### plot FE surface on left plot
        # Define the number of bins for the 2D histogram
        num_bins = 20

        # Calculate the 2D histogram
        hist, xedges, yedges = np.histogram2d(x, y, bins=num_bins, density=True)

        # Calculate the probability as the histogram values
        probability = hist / np.sum(hist)
        #print(f'probability: {probability} {np.unique(probability)}')

        # Compute the free energy as -log10(probability)
        with np.errstate(divide='ignore'):  # Ignore divide-by-zero warnings
            free_energy = -np.log10(probability)
            #free_energy[np.isinf(free_energy)] = np.nan  # Set infinities to NaN for better plotting
            free_energy[np.isinf(free_energy) | np.isnan(free_energy)] = np.nanmax(free_energy[np.isfinite(free_energy)]) #+ 1  # Replace NaN/Inf with a large value
        print(f'free_energy: {free_energy} {free_energy.shape} {np.unique(free_energy)}')

        # Create the meshgrid for the contour plot
        X, Y = np.meshgrid(xedges[:-1], yedges[:-1])
        print(f'X: {X.shape}\nY: {Y.shape}')

        # Create a custom colormap
        #cmap = plt.cm.viridis
        #cmap = plt.cm.magma
        cmap = plt.cm.gist_ncar

        # Plotting the contour plot
        contour = axes[0].contourf(X, Y, free_energy.T, levels=100, cmap=cmap)  # Transpose to align axes
        fig.colorbar(contour, ax=axes[0], label='Free Energy (-log10 Probability)')
        axes[0].set_xlabel('Q')
        axes[0].set_ylabel('G')
        axes[0].set_title('2D Free Energy Contour Plot')
        axes[0].set_xlim(0,1)


        #############################################################################################
        ## Plot state map
        #_, axes[1], _ = pem.plots.plot_state_map(x, y, states)

        # Create a 2D histogram to determine the bin index for each (x, y) pair
        #############################################################################################
        # Step 1: Identify unique states
        unique_states = np.unique(states)
        n_states = len(unique_states)
        print(f'unique_states: {unique_states} {n_states}')

        # Step 2: Create a colormap with one color per unique state
        # You can use any colormap, or define specific colors if desired
        colors = plt.cm.get_cmap('tab20', n_states)  # 'tab10' has up to 10 colors; change if needed
        cmap = ListedColormap([colors(i) for i in range(n_states)])

        # Step 3: Map states to color indices
        state_to_index = {state: i for i, state in enumerate(unique_states)}
        color_indices = np.vectorize(state_to_index.get)(states)

        # # Step 4: Create scatter plot
        scatter = axes[1].scatter(x, y, c=color_indices, cmap=cmap, s=50, edgecolor='k')  # Customize marker size, etc.

        # Step 5: Add a colorbar with labels
        cbar = plt.colorbar(scatter, ax=axes[1], ticks=np.linspace(0.5, n_states - 1.5, num=n_states), label=f'Metastable States')
        cbar.ax.set_yticklabels(unique_states)  # Label colorbar with the unique state values
        #############################################################################################

        axes[1].set_xlabel('Q')
        axes[1].set_ylabel('G')
        axes[1].set_title('2D state map')
        axes[1].set_xlim(0,1)

        #plt.tight_layout()
        plt.savefig(outfile)
        logging.info(f'SAVED: {outfile}')
        print(f'SAVED: {outfile}')
        plt.clf()

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
        n_states = len(np.unique(self.dtrajs))  # or a predefined number of states
        its = pem.msm.its(self.dtrajs, lags=lag_times, errors='bayes')
        pem.plots.plot_implied_timescales(its)
        ITS_outfile = os.path.join(self.outpath, f'{self.outname}_ITS.png')
        plt.savefig(ITS_outfile)
        logging.info(f'SAVED: {ITS_outfile}')
        print(f'SAVED: {ITS_outfile}')
    #######################################################################################

    #######################################################################################
############## MAIN #################
def main():
    
    script_name = f'BuildKineticModel'
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("--OPpath", type=str, required=True, help="Path to directory containing G and Q directories created by GQ.py")
    parser.add_argument("--outname", type=str, required=True, help="base name for output files")
    parser.add_argument("--start", type=int, required=False, help="First frame to analyze 0 indexed", default=0)
    parser.add_argument("--end", type=int, required=False, help="Last frame to analyze 0 indexed", default=9999999)
    parser.add_argument("--stride", type=int, required=False, help="Frame stride", default=1)
    parser.add_argument("--ITS", type=str, required=False, help="Find optimal lag time with ITS", default='False')
    parser.add_argument("--lagtime", type=int, required=False, help="lagtime to build the model", default=1)
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
    if args.ITS == 'True':
        anal.plot_implied_timescales()
        print(f'Analysis terminated since ITS was selected. Check the figure and choose an approrate lagtime')
        quit()
    
    # Build the MSM model with the choosen lagtime
    anal.build_msm(lagtime=args.lagtime)

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    
#print(f'NORMAL TERMINATION: {time.time() - start_time}')
#logging.info(f'NORMAL TERMINATION: {time.time() - start_time}')