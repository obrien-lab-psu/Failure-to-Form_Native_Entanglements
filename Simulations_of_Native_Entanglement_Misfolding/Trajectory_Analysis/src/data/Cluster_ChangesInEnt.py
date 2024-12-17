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
from scipy.spatial.distance import pdist, squareform

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
        self.EntInfo_df = pd.read_csv(self.EntInfofile, low_memory=False)
        #print(f'EntInfo_df:\n{self.EntInfo_df}')     

        # convert all NaN values in the crossings columns to 0
        self.EntInfo_df['crossingsN'] = self.EntInfo_df['crossingsN'].fillna(0)
        self.EntInfo_df['crossingsC'] = self.EntInfo_df['crossingsC'].fillna(0)
        #print(self.EntInfo_df[self.EntInfo_df['Frame'] == 8154])

    #######################################################################################  

    #######################################################################################
    def agg_method(self, df, ref_df):
        """
        This is a clustering method that first clusters based on median crossing and then loop overlap. 
        See https://github.com/obrien-lab/cg_simtk_protein_folding/wiki/find_representative_chg_ent.py for more details
        """
        #print(df)
        frame = df['Frame'].values[0]

        df.loc[:, 'crossingsN'] = df['crossingsN'].astype(str)
        df.loc[:, 'crossingsC'] = df['crossingsC'].astype(str)

        # count the number of crossings
        df['len_crossingsN'] = df['crossingsN'].apply(crossings_length)
        df['len_crossingsC'] = df['crossingsC'].apply(crossings_length)
        #print(df)

        cID = 0
        df['cID'] = -1
        for change_vec, change_df in df.groupby(['Gn', 'NchangeType', 'Gc', 'CchangeType', 'len_crossingsN', 'len_crossingsC']):

            ######################################################
            ## get the i, j, and median crossings array for further agglomerative cluststering
            #print(f'change_vec: {change_vec}')
            ijr = []
            for rowi, row in change_df.iterrows():
                i = row['i']
                j = row['j']
                Ncrossings = row['crossingsN']
                Ccrossings = row['crossingsC']

                #print(i, j, Ncrossings, Ccrossings)
                crossings = []
                if row['NchangeType'] != 'NoChange':
                    if isinstance(Ncrossings, float):
                        crossings += [float(Ncrossings)]
                    else:
                        Ncrossings = str(Ncrossings).split(',')
                        crossings += [float(c) for c in Ncrossings]    

                if row['CchangeType'] != 'NoChange':
                    if isinstance(Ccrossings, float):
                        crossings += [float(Ccrossings)]
                    else:  
                        Ccrossings = str(Ccrossings).split(',')
                        crossings += [float(c) for c in Ccrossings]   

                #print(i, j, crossings, np.median(crossings))   
                ijr += [[i, j, np.median(crossings)]]  
            ijr = np.asarray(ijr)
            #print(ijr)

            ######################################################
            ## cluster based on crossing medians
            crossing_dists = ijr[:,2].reshape(-1, 1)
            distances = squareform(pdist(crossing_dists, metric='euclidean'))
            #print(distances)
            db = DBSCAN(eps=10, min_samples=1, metric='precomputed').fit(distances)
            labels = db.labels_
            #print(f'labels: {labels} {np.unique(labels)}')

            ######################################################
            ## within each crossing cluster cluster based on loop distance
            ijr = np.hstack([ijr, labels[:, None]])
            for crossing_cluster in np.unique(labels):
                cluster_ijr = ijr[np.where(ijr[:,-1] == crossing_cluster)]
                #print(cluster_ijr)

                db = DBSCAN(eps=1, min_samples=1, metric=custom_distance).fit(cluster_ijr[:, :2])
                labels = db.labels_
                #print(f'labels: {labels} {np.unique(labels)}')

                cluster_ijr = np.hstack([cluster_ijr, labels[:, None]])
                #print(cluster_ijr)
                
                ## update the final dataframe
                for loop_cluster in np.unique(labels):
                    loop_ijr = cluster_ijr[np.where(cluster_ijr[:,-1] == loop_cluster)]
                    #print(loop_ijr)
           
                    for i,j in loop_ijr[:, :2]:
                        df.loc[(df['i'] == i) & (df['j'] == j), 'cID'] = cID
                    cID += 1

        #########################################################
        ### check for cross contamination
        # have not had to deal with this in my dataset since I cluster on the frame level
        # if someone comes after you can choose to divisivly cluster if necessary
        for cID, cID_df in df.groupby('cID'):
            #print(cID_df)
            ij = cID_df[['i', 'j']].values
            Nchange = cID_df['NchangeType'].values[0]
            Cchange = cID_df['CchangeType'].values[0]
            contamination = []
            if Nchange != 'NoChange':
                for cross in cID_df['crossingsN']:
                    cross = str(cross)
                    for c in cross.split(','):
                        c = abs(float(c))
                        for i,j in ij:
                            if c > i and c < j:
                                d = min([c - i, j - c])/(j-i)
                                if d > 0.1:
                                    contamination += [d]
            if Cchange != 'NoChange':
                for cross in cID_df['crossingsC']:
                    cross = str(cross)
                    for c in cross.split(','):
                        c = abs(float(c))
                        for i,j in ij:
                            if c > i and c < j:
                                d = min([c - i, j - c])/(j-i)
                                if d > 0.1:
                                    contamination += [d]
            if len(contamination) > 0:
                print(df.to_string())
                print(cID_df)
                print(f'contamination: {contamination}')
                #raise ValueError(f'Found contamination in frame: {frame} {cID}')

        return df

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
        frame_N = 0
        newGs = {'Frame':[], 'G':[]}
        for frame, frame_df in self.EntInfo_df.groupby('Frame'):

            frame_df = frame_df.copy()

            # skip the -1 frame as it is a reference structure
            if frame == -1:
                #print(f'Reference:\n{frame_df}')
                frame_df['cID'] = -1
                ref_df = frame_df.copy()
                ref_N = len(ref_df)
                #print(ref_df.to_string())
                #clustered_dfs += [ref_df]
                continue

            # check if the frame is outside the specified window
            if frame < self.start or frame > self.end:
                #print(f'Frame {frame} is outside the bounds [{self.start, self.end}]')
                continue

            # check if the frame has no change
            NchangeType, CchangeType = frame_df['NchangeType'].value_counts(), frame_df['CchangeType'].value_counts()
            if len(NchangeType) == 1 and len(CchangeType) == 1:
                if NchangeType['NoChange'] and CchangeType['NoChange']:
                    newGs['Frame'] += [frame]
                    newGs['G'] += [0]                    
                    continue

            # check if the frame has no contacts as can be the case in the first frame of quenching sometimes
            if len(frame_df) == 1 and all(frame_df[['i', 'j']].isnull()):
                #print(f'No native contacts found this frame {frame} and will skip')
                newGs['Frame'] += [frame]
                newGs['G'] += [0]     
                continue

            ########################################################################################
            ## Step 1: get only the change row
            frame_df = frame_df[(frame_df['NchangeType'] != 'NoChange') | (frame_df['CchangeType'] != 'NoChange')]
            #print(frame_df)

            ## Step 2: check for phantom entanglement changes by calculating the deltagn and deltagc metrics and removing those below a certain threshold (0.35)
            frame_df = self.check_for_phantoms(frame_df, ref_df)
            #print(frame_df.to_string())

            ## Step 3: if any gains do not have a crossing identified by topoly in the frame or any loses do not have a crossing identified by topoly in the ref then remove
            frame_df = self.check_for_crossings(frame_df, ref_df)
            #print(frame_df.to_string())

            ## Step 4: IF there are remaining changes after controling for the phatom changes then cluster
            if len(frame_df) != 0:

                frame_df = self.agg_method(frame_df, ref_df)
                #print(frame_df.to_string())

                frame_N = len(frame_df)
                clustered_dfs += [frame_df]

            newG = frame_N/ref_N
            newGs['Frame'] += [frame]
            newGs['G'] += [newG]

        if len(clustered_dfs) == 0:
            print(f'No changes in ENT found in this traj')
        else:
            clustered_dfs = pd.concat(clustered_dfs)
            #print(f'clustered_dfs:\n{clustered_dfs}')

            outfile = os.path.join(self.outpath, f'{self.outname}_clustered.EntInfo')
            clustered_dfs.to_csv(outfile, index=False)
            logging.info(f'SAVED: {outfile}')
            print(f'SAVED: {outfile}')

        # save the new G file
        outfile = os.path.join(self.outpath, f'{self.outname}_clustered.G')
        newGs = pd.DataFrame(newGs)
        newGs.to_csv(outfile, index=False)
        print(f'SAVED: {outfile}')

    #######################################################################################

    #######################################################################################
    def check_for_phantoms(self, df, ref_df, threshold=0.35):
        """
        For a given EntInfo dataframe determine if the absolute value of the change in linking value is greater than 0.35 and only pass those that are
        """
        delta_gns = []
        delta_gcs = []
        for rowi, row in df.iterrows():
            i, j = row['i'], row['j']
            gn, gc = row['gn'], row['gc']
            Nchange, Cchange = row['NchangeType'], row['CchangeType']

            ref_rows = ref_df[(ref_df['i'] == i) & (ref_df['j'] == j)]
            ref_gn, ref_gc = ref_rows['gn'].values[0], ref_rows['gc'].values[0]

            if Nchange != 'NoChange':
                delta_gns += [abs(gn - ref_gn)]
            else:
                delta_gns += [99]

            if Cchange != 'NoChange':  
                delta_gcs += [abs(gc - ref_gc)]  
            else:
                delta_gcs += [99]

        #df['delta_gn'] = delta_gns
        #df['delta_gc'] = delta_gcs
        df.loc[:, 'delta_gn'] = delta_gns
        df.loc[:, 'delta_gc'] = delta_gcs
        df = df[(df['delta_gn'] > threshold) & (df['delta_gc'] > threshold)]
        return df
    #######################################################################################

    #######################################################################################
    def check_for_crossings(self, df, ref_df):
        """
        This function checks a frame dataframe for topoly crossings. 
        If it is a gain and there is no crossing in the frame then set changeType to NoChange
        If it is a loss and there is no crossing in the ref then set changeType to NoChange
        """

        df.loc[:, 'crossingsN'] = df['crossingsN'].astype(str)
        df.loc[:, 'crossingsC'] = df['crossingsC'].astype(str)
        
        ## if loss is present get native crossings
        for rowi, row in df.iterrows():
            i = row['i']
            j = row['j']
            Ncrossings = str(row['crossingsN'])
            Ccrossings = str(row['crossingsC'])

            Gn_sign = np.sign(row['Gn'])
            if Gn_sign > 0:
                Gn_sign = '+'
            elif Gn_sign < 0:
                Gn_sign = '-'
            else:
                Gn_sign = ''

            Gc_sign = np.sign(row['Gc'])
            if Gc_sign > 0:
                Gc_sign = '+'
            elif Gc_sign < 0:
                Gc_sign = '-'
            else:
                Gc_sign = ''


            ref_ij = ref_df[(ref_df['i'] == i) & (ref_df['j'] == j)]
            ref_Ncrossings = str(ref_ij['crossingsN'].values[0])
            ref_Ccrossings = str(ref_ij['crossingsC'].values[0])

            ref_Gn_sign = np.sign(ref_ij['Gn'].values[0])
            if ref_Gn_sign > 0:
                ref_Gn_sign = '+'
            elif ref_Gn_sign < 0:
                ref_Gn_sign = '-'
            else:
                ref_Gn_sign = ''

            ref_Gc_sign = np.sign(ref_ij['Gc'].values[0])
            if ref_Gc_sign > 0:
                ref_Gc_sign = '+'
            elif ref_Gc_sign < 0:
                ref_Gc_sign = '-'
            else:
                ref_Gc_sign = ''
            #print(i, j, Ncrossings, Ccrossings, [Gn_sign, Gc_sign], ref_Ncrossings, ref_Ccrossings, [ref_Gn_sign, ref_Gc_sign])

            ## if loss was detected then replace crossing with that from reference state
            # if there is a * make sure to add the chirality using the Gn or Gc value as a proxy
            if 'Loss' in row['NchangeType']:
                if ref_Ncrossings == '0':
                    df.loc[rowi, 'NchangeType'] = 'NoChange'
                else:
                    df.loc[rowi, 'crossingsN'] = ref_Ncrossings.replace('*', ref_Gn_sign)
            elif 'Gain' in row['NchangeType']:
                if Ncrossings == '0':
                    df.loc[rowi, 'NchangeType'] = 'NoChange'
                else:
                    df.loc[rowi, 'crossingsN'] = Ncrossings.replace('*', Gn_sign)
            elif 'Pure' in row['NchangeType']:
                df.loc[rowi, 'crossingsN'] = Ncrossings.replace('*', Gn_sign)
            elif 'NoChange' in row['NchangeType']:
                df.loc[rowi, 'crossingsN'] = Ncrossings.replace('*', Gn_sign)

            if 'Loss' in row['CchangeType']:
                if ref_Ccrossings == '0':
                    df.loc[rowi, 'CchangeType'] = 'NoChange'
                else:
                    df.loc[rowi, 'crossingsC'] = ref_Ccrossings.replace('*', ref_Gc_sign)
            elif 'Gain' in row['CchangeType']:
                if Ccrossings == '0':
                    df.loc[rowi, 'CchangeType'] = 'NoChange'
                else:
                    df.loc[rowi, 'crossingsC'] = Ccrossings.replace('*', Gc_sign)
            elif 'Pure' in row['CchangeType']:
                df.loc[rowi, 'crossingsC'] = Ccrossings.replace('*', Gn_sign)
            elif 'NoChange' in row['CchangeType']:
                df.loc[rowi, 'crossingsC'] = Ccrossings.replace('*', Gc_sign)

        df = df[(df['NchangeType'] != 'NoChange') | (df['CchangeType'] != 'NoChange')]
        return df
    #######################################################################################

# Define a custom distance function
def custom_distance(point1, point2):
    # Example: Weighted Euclidean distance
    point1_L = np.arange(point1[0], point1[1] + 1)
    point2_L = np.arange(point2[0], point2[1] + 1)
    combined = np.hstack([point1_L, point2_L])
    d = (max(combined) - min(combined))/(len(point1_L) + len(point2_L))
    #print(point1, point2, d)
    return d

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

    # Step 4: updated G values 


if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    
print(f'NORMAL TERMINATION: {time.time() - start_time}')
logging.info(f'NORMAL TERMINATION: {time.time() - start_time}')