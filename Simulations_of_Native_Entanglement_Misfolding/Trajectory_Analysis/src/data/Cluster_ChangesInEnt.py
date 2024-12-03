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
from scipy.stats import mode, permutation_test
from sklearn.datasets import make_blobs
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
from sklearn.cluster import DBSCAN
pd.set_option('display.max_rows', 5000)

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
        ("--EntInfofile", type=str, required=True, help="file containing change in entanglement info from GQ.py")
        ("--outname", type=str, required=True, help="base name for output files")
        """

        # parse the parameters 
        self.EntInfofile = args.EntInfofile
        #logging.info(f'EntInfofile: {self.EntInfofile}')

        self.outpath = args.outpath
        logging.info(f'outpath: {self.outpath}')

        self.outname = args.outname
        logging.info(f'outname: {self.outname}')

        self.start = args.start
        self.end = args.end
        print(f'START: {self.start} | END: {self.end}')

    #######################################################################################

    #######################################################################################
    def load_OP(self,):
        """
        """
        self.EntInfo_df = pd.read_csv(self.EntInfofile)
        #print(f'EntInfo_df:\n{self.EntInfo_df}')     

        # add column with counts of unique crossings for N and C terminus
        self.EntInfo_df['len_crossingsN'] = self.EntInfo_df['crossingsN'].apply(crossings_length)
        self.EntInfo_df['len_crossingsC'] = self.EntInfo_df['crossingsC'].apply(crossings_length)
        print(f'EntInfo_df:\n{self.EntInfo_df}')  
        #print(self.EntInfo_df[self.EntInfo_df['Frame'] == 8154])

        # convert all NaN values in the crossings columns to 0
        self.EntInfo_df['crossingsN'] = self.EntInfo_df['crossingsN'].fillna(0)
        self.EntInfo_df['crossingsC'] = self.EntInfo_df['crossingsC'].fillna(0)
        #print(self.EntInfo_df[self.EntInfo_df['Frame'] == 8154])



    #######################################################################################  


    #######################################################################################
    def Cluster(self,):
        """
        Clustering steps
        Step 1: get only the change row 
        Step 2: DBSCAN cluster the ij array


        input dataframe stucture from .EntInfo file
                Time(ns)  Frame    i    j       gn   Gn  crossingsN       gc   Gc  crossingsC NchangeType CchangeType
        298365  168.900004   2251    3    7  0.00000  0.0         NaN -0.01805 -0.0         NaN    NoChange    NoChange
        298369  168.900004   2251    7   58  0.00000  0.0         NaN  0.10698  0.0         NaN    NoChange  LossChiral
        

        """
        # Define the output dataframe for the clusters
        clustered_dfs = []

        ########################################################################################
        # Loop through frames and cluster the changes in ent in each
        for frame, frame_df in self.EntInfo_df.groupby('Frame'):
            # skip the -1 frame as it is a reference structure
            if frame == -1:
                print(f'Reference:\n{frame_df}')
                frame_df['cID'] = -1
                ref_df = frame_df.copy()
                clustered_dfs += [ref_df]
                continue

            # check if the frame is outside the specified window
            if frame < self.start or frame > self.end:
                #print(f'Frame {frame} is outside the bounds [{self.start, self.end}]')
                continue

            # check if the frame has no change
            NchangeType, CchangeType = frame_df['NchangeType'].value_counts(), frame_df['CchangeType'].value_counts()
            if len(NchangeType) == 1 and len(CchangeType) == 1:
                if NchangeType['NoChange'] and CchangeType['NoChange']:
                    continue

            # check if the frame has no contacts as can be the case in the first frame of quenching sometimes
            if len(frame_df) == 1 and all(frame_df[['i', 'j']].isnull()):
                print(f'No native contacts found this frame {frame} and will skip')
                continue

            ########################################################################################
            ## Step 1: get only the change row
            frame_df = frame_df[(frame_df['NchangeType'] != 'NoChange') | (frame_df['CchangeType'] != 'NoChange')]
            print(frame_df)


            ########################################################################################
            ## Step 2: scan for changes of entanglement that are around the threshold of 0.6
            ref_rows = ref_df.merge(frame_df[['i', 'j']], on=['i', 'j'], how='inner')
            #print(ref_rows)

            if len(frame_df) != len(ref_rows):
                raise ValueError(f'The number of rows found in the reference df does not match those in the frame df {len(frame_df)} != {len(ref_rows)}\nframe_df:\n{frame_df}\nref_rows:\n{ref_rows}')

            rowidx_to_drop = []
            for idx, (rowi, row) in enumerate(frame_df.iterrows()):
                #print(idx, rowi)
                ref_row = ref_rows.iloc[idx]

                #QC to ensure the ij match for the rows being compared
                if row['i'] != ref_row['i'] or row['j'] != ref_row['j']:
                    raise ValueError(f'There is a mismatch in the native contacts between the row in the frame df and the row in the reference df')
                
                # check if there is a N terminal change and if it was around the threshold
                if row['NchangeType'] != 'NoChange':
                    if row['gn'] < 0:
                        lb, ub = -0.7, -0.5
                    elif row['gn'] > 0:
                        lb, ub = 0.5, 0.7

                    if abs(row['gn']) < ub and abs(row['gn']) > lb:
                        if abs(ref_row['gn']) < ub and abs(ref_row['gn']) > lb:
                            #print(f'Found a potential false change due to threshold at rowi {rowi}')
                            rowidx_to_drop += [rowi]

            #print(f'rowidx_to_drop: {rowidx_to_drop}')
            if len(rowidx_to_drop) != 0:
                frame_df = frame_df.drop(rowidx_to_drop)
            
            # check that if after the removal of changes around the threshold that there are no ents left. if so skip this frame
            if len(frame_df) == 0:
                continue

            ########################################################################################
            ## Step 3: DBSCAN cluster the ij array
            ij = frame_df[['i', 'j']].values
            print(f'ij:\n{ij}')
            db = DBSCAN(eps=20, min_samples=1).fit(ij)
            labels = db.labels_
            print(f'labels: {labels} {np.unique(labels)}')
            
            ## Step 4: for each cluster of loop closing contacts with change in entanglement
            # (1) cluster based on the following vector [Gn, NchangeType, Gc, CchangeType]
            # (2) for each of those resulting clusters cluster based on the 
            cID = 0
            for label in np.unique(labels):

                # get the rows in the ij cluster
                label_idx = np.where(labels == label)
                label_df = frame_df.iloc[label_idx]

                # group by [Gn, NchangeType, Gc, CchangeType]
                for change_vec, change_df in label_df.groupby(['Gn', 'NchangeType', 'Gc', 'CchangeType', 'len_crossingsN', 'len_crossingsC']):

                    # group by the number of N and C crossings 
                    #for crossing_vec, crossings_df in change_df.groupby(['len_crossingsN', 'len_crossingsC']):

                    # DBSCAN the crossings to get the final clustered change in ent
                    #crossings = crossings_df[['crossingsN', 'crossingsC']].values
                    crossings = change_df[['crossingsN', 'crossingsC']].values
                    crossings = processes_crossings_array(crossings)
                    cross_db = DBSCAN(eps=10, min_samples=1).fit(crossings)
                    cross_labels = cross_db.labels_
                    #print(f'cross_labels: {cross_labels} {np.unique(cross_labels)}')

                    # for each final group of raw changes in ent prepare the output data
                    for cross_label in np.unique(cross_labels):
                        #cross_label_df = crossings_df.iloc[np.where(cross_labels == cross_label)]
                        cross_label_df = change_df.iloc[np.where(cross_labels == cross_label)]
                        cross_label_df['cID'] = cID
                        cID += 1
                        clustered_dfs += [cross_label_df]

        clustered_dfs = pd.concat(clustered_dfs)
        print(f'clustered_dfs:\n{clustered_dfs}')

        outfile = os.path.join(self.outpath, f'{self.outname}_clustered.EntInfo')
        clustered_dfs.to_csv(outfile, index=False)
        logging.info(f'SAVED: {outfile}')
        print(f'SAVED: {outfile}')

 
    #######################################################################################

#######################################################################################
# Define function to calculate length
def crossings_length(value):
    if pd.isna(value):                # Check for NaN
        return 0
    elif isinstance(value, float):    # Check if it's a single float
        return 1
    else:                             # Otherwise, split on comma and count items
        return len(str(value).split(','))
#######################################################################################

#######################################################################################
def processes_crossings_array(arr):
    """
    for each row determine if any strings are present. if so expand and convert to int
    """
    newarr = []
    for row in arr:
        newrow = []
        for elem in row:
            if isinstance(elem, str):
                elem = elem.split(',')
                elem = [abs(int(e.replace('*', ''))) for e in elem]
                newrow += elem
            else:
                newrow += [abs(int(elem))]
        newarr += [newrow]
    
    return np.asarray(newarr)

#######################################################################################

#######################################################################################
def statistic(x, y, axis):
    return np.mean(x, axis=axis) - np.mean(y, axis=axis)
#######################################################################################

############## MAIN #################
def main():
    """
    This script is designed to take the .EntInfo file generated by GQ.py and cluster the resulting changes in entanglement within each frame. 

    """
    script_name = f'Cluster_ChangesInEnt'
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("--EntInfofile", type=str, required=True, help="file containing change in entanglement info from GQ.py")
    parser.add_argument("--outname", type=str, required=True, help="base name for output files")
    parser.add_argument("--start", type=int, required=False, help="First frame to analyze 0 indexed", default=0)
    parser.add_argument("--end", type=int, required=False, help="Last frame to analyze 0 indexed", default=99999999999)
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

    # Step 0: initialize the simulation object 
    anal = Analysis(args)

    # Step 1: load the G and Q data
    anal.load_OP()

    # Step 2: screen for changes that are likely to be phantom entanglements
    #anal.Screen()
    #quit()

    # Step 3: for each candidate get the native G and Q states
    anal.Cluster()


if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    
print(f'NORMAL TERMINATION: {time.time() - start_time}')
logging.info(f'NORMAL TERMINATION: {time.time() - start_time}')