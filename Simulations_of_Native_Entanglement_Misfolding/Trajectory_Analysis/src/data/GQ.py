#!/usr/bin/env python3
import requests, logging, os, sys
import time
import argparse
import pandas as pd
import numpy as np
import glob
import MDAnalysis as mda
from scipy.spatial.distance import pdist, squareform
from topoly import lasso_type  # used pip
import itertools
import concurrent.futures

pd.set_option('display.max_rows', 5000)

class Analysis:
    """
    A class to handel the analyssis of a C-alpha CG trajectory. 
    Current analysis available:
    (1) - Fraction of native contacts (Q)
    (2) - Change of entanglement (G)
    """
    #######################################################################################
    def __init__(self, args):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
        ("--outpath", type=str, required=True, help="Path to output directory")
        ("--outname", type=str, required=True, help="base name for output files")
        ("--log", type=str, required=True, help="Path to logging file")
        ("--psf", type=str, required=True, help="Path to CA protein structure file")
        ("--cor", type=str, required=True, help="Path to CA native coordinates file")
        ("--dcd", type=str, required=True, help="Path to trajectory to analyze")
        ("--sec_elements", type=str, required=True, help="Path to STRIDE secondary structure elements file")
        ("--start", type=int, required=False, help="First frame to analyze 0 indexed", default=0)
        ("--end", type=int, required=False, help="Last frame to analyze 0 indexed", default=-1)
        ("--stride", type=int, required=False, help="Frame stride", default=1)
        """

        # parse the parameters 
        self.outpath = args.outpath
        logging.info(f'outpath: {self.outpath}')

        self.outname = args.outname
        logging.info(f'outname: {self.outname}')

        self.psf = args.psf
        logging.info(f' psf: {self. psf}')

        self.cor = args.cor
        logging.info(f'cor: {self.cor}')

        self.dcd = args.dcd 
        logging.info(f'dcd: {self.dcd}')

        self.sec_elements = args.sec_elements
        logging.info(f'sec_elements: {self.sec_elements}')

        self.ref_universe = mda.Universe(self.psf, self.cor, format='CRD')
        print(f'ref_universe: {self.ref_universe}')

        self.traj_universe = mda.Universe(self.psf, self.dcd, format='DCD')
        print(f'traj_universe: {self.traj_universe}')

        self.start = args.start
        self.end = args.end
        self.stride = args.stride
        print(f'START: {self.start} | END: {self.end} | STRIDE: {self.stride}')
    #######################################################################################

    #######################################################################################
    def Q(self, ):
        print(f'Calculating the fraction of native contacts (Q)')
        logging.info(f'Calculating the fraction of native contacts (Q)')
        """
        Calculate the fraction of native contacts in each frame of the DCD where a native contact is defined between secondary structures 
        and for residues atleast that are atleast 3 residues apart. So if i = 1 then j at a minimum can be 5. 
        For a contact to be present the distance between i and j must be less than 8A in the native structure and in a trajectory frame be less than 1.2*native distance.
        """
        # make directory for Q data if it doesnt exist
        self.Qpath = os.path.join(self.outpath, 'Q')
        if not os.path.exists(self.Qpath):
            os.makedirs(self.Qpath)
            print(f'Made directory: {self.Qpath}')

        # Step 0: load the reference structure and topology
        ref_coor = self.ref_universe.atoms.positions
        #print(f'ref_coor:\n{ref_coor} {ref_coor.shape}')        


        # Step 1: Get the secondary structure information
        # get both those resides in the secondary structures and those not
        print(f'Step 1: Get the secondary structure information')
        resid_in_sec_elements = np.loadtxt(self.sec_elements, dtype=int)
        resid_in_sec_elements = [np.arange(x[1], x[2] + 1) for x in resid_in_sec_elements]
        resid_in_sec_elements = np.hstack(resid_in_sec_elements)
        logging.info(f'resid_in_sec_elements: {resid_in_sec_elements}')

        resid_not_in_sec_elements = np.asarray([r for r in range(1, len(ref_coor) + 1) if r not in resid_in_sec_elements]) # residue ID not in secondary structures
        logging.info(f'resid_not_in_sec_elements: {resid_not_in_sec_elements}')


        # Step 2: Get the native distance map for the native state cordinates
        print(f'Step 2: Get the native distance map for the native state cordinates')
        # Zero the resulting distance map up to the 4th diagonal so only those residues with more than 3 residues between them can be in contact
        # Zero out any secondary structure element residues
        # Zero out any distance not less than 8A
        ref_distances = np.triu(squareform(pdist(ref_coor)), k=4)
        ref_distances[resid_not_in_sec_elements - 1, :] = 0
        ref_distances[:, resid_not_in_sec_elements - 1] = 0
        ref_distances[ref_distances > 8] = 0
        NumNativeContacts = np.count_nonzero(ref_distances)
        print(f'NumNativeContacts: {NumNativeContacts}')
        logging.info(f'NumNativeContacts: {NumNativeContacts}')

        # Step 3: Analyze each frame of the traj_universe and get the distance map
        print(f'Step 3: Analyze each frame of the traj_universe and calc Q')
        # then determine the fraction of native contacts by those distances less than 1.2*native distance
        Qoutput = {'Time(ns)':[], 'Frame':[], 'FrameNumNativeContacts':[], 'Q':[]}
        for ts in self.traj_universe.trajectory[self.start:self.end:self.stride]:
            frame_coor = self.traj_universe.atoms.positions
            frame_distances = np.triu(squareform(pdist(frame_coor)), k=4)
            frame_distances[resid_not_in_sec_elements - 1, :] = 0
            frame_distances[:, resid_not_in_sec_elements - 1] = 0

            cond = (frame_distances <= 1.2*ref_distances) & (ref_distances != 0)

            FrameNumNativeContacts = np.sum(cond)
            #print(f'FrameNumNativeContacts: {FrameNumNativeContacts} for frame {ts.frame}')

            Q = FrameNumNativeContacts/NumNativeContacts
            #print(f'Q: {Q} for frame {ts.frame}')

            frame_time = ts.time/1000
            Qoutput['Frame'] += [ts.frame]
            Qoutput['FrameNumNativeContacts'] += [FrameNumNativeContacts]
            Qoutput['Q'] += [Q]
            Qoutput['Time(ns)'] += [frame_time]
        
        # Step 4: save Q output 
        print(f'Step 4: save Q output')
        Qoutput = pd.DataFrame(Qoutput)
        Qoutfile = os.path.join(self.Qpath, f'{self.outname}.Q')
        Qoutput.to_csv(Qoutfile, index=False)
        print(f'SAVED: {Qoutfile}')
        logging.info(f'SAVED: {Qoutfile}')
    #######################################################################################  

    #######################################################################################
    def G(self,):
        print(f'Calculating the fraction of native contacts with a change in Gauss Linking Number (G)')
        logging.info(f'Calculating the fraction of native contacts with a change in Gauss Linking Number (G)')
        """
        Calculate the fraction of native contacts that have a change in entanglement status. 
        These native contacts are not restricted to secondary structures.
        """
        # make directory for Q data if it doesnt exist
        self.Gpath = os.path.join(self.outpath, 'G')
        if not os.path.exists(self.Gpath):
            os.makedirs(self.Gpath)
            print(f'Made directory: {self.Gpath}')

        # Step 1: load the reference structure and topology and get reference state GaussLinking numbers
        print(f'Step 1: load the reference structure and topology and get reference state GaussLinking numbers')
        ref_coor = self.ref_universe.atoms.positions
        #print(f'ref_coor:\n{ref_coor} {ref_coor.shape}')        
        ref_nc_gdict = self.GaussLink(ref_coor)
        Nnative = len(ref_nc_gdict)
        print(f'Nnative: {Nnative}')
        logging.info(f'Nnative: {Nnative}')

        # Step 1a: Initiate the output dataframe that will contain all the detailed change in entanglement information. 
        # the native structure data will be populated into Frame -1 and Time(ns) -1 for reference 
        EntOutput = {'Time(ns)':[], 'Frame':[], 'i':[], 'j':[], 'gn':[], 'Gn':[], 'crossingsN':[], 'gc':[], 'Gc':[], 'crossingsC':[], 'NchangeType':[], 'CchangeType':[]}
        nc_list = []
        for nc, ref_info in ref_nc_gdict.items():
            i,j = nc[0], nc[1]
            nc_list += [[i, j]] # update this list for the contact_mask when looping through the trajectory data
            EntOutput['Time(ns)'] += [-1]
            EntOutput['Frame'] += [-1]
            EntOutput['i'] += [i]
            EntOutput['j'] += [j]
            EntOutput['gn'] += [ref_info['gn']]
            EntOutput['Gn'] += [ref_info['Gn']]
            EntOutput['crossingsN'] += [ref_info['crossingsN']]
            EntOutput['gc'] += [ref_info['gc']]
            EntOutput['Gc'] += [ref_info['Gc']]
            EntOutput['crossingsC'] += [ref_info['crossingsC']]
            EntOutput['NchangeType'] += ['None-Ref']
            EntOutput['CchangeType'] += ['None-Ref']
            

        # Step 1b: Initiate the output dataframe that will contain the high level stats for changes in entanglement and the G metric used in higherarchical modeling
        Goutput = {'Time(ns)':[], 'Frame':[], 
                    'Nnative':[], 'N-Loss':[], 'N-LossChiral':[], 'N-Gain':[], 'N-GainChiral':[], 'N-PureChiral':[], 
                    'C-Loss':[], 'C-LossChiral':[], 'C-Gain':[], 'C-GainChiral':[], 'C-PureChiral':[], 'G':[]}

        # Step 2: loop through the frames and calculate the changes in entanglement
        print(f'Step 2: loop through the frames and calculate the changes in entanglement')
        cpu_cores = os.cpu_count()
        print(f"Number of available CPU cores: {cpu_cores}")

        # Define the input data for each frame
        frames_data = []
        for ts in self.traj_universe.trajectory[self.start:self.end:self.stride]:
            frame_coor = self.traj_universe.atoms.positions
            frame_data = (
                frame_coor,                # Current frame coordinates
                nc_list,                   # The contact mask (nc_list)
                ref_nc_gdict,              # The reference nc_gdict
                ts.time/1000,              # Time in ns (dividing by 1000)
                ts.frame,                  # Frame index
                self.GaussLink,            # Function GaussLink
                self.GetLinkChanges,         # Function GetLinkChanges
                Nnative                     # number of native contacts       
            )
            frames_data.append(frame_data)

        # Use ProcessPoolExecutor to parallelize the work
        with concurrent.futures.ProcessPoolExecutor(max_workers=cpu_cores) as executor:
            # Map the process_frame function to the frames_data
            results = executor.map(process_frame, frames_data)

            # Collect the results
            for change_info, count_info in results:

                # Adding values of change_info to the master EntOutput dictionary
                for key in EntOutput:
                    EntOutput[key] += change_info[key]

                frame = change_info['Frame'][0]
                frame_time = change_info['Time(ns)'][0]

                # Adding values of count_info to the master Goutput dictionary
                Goutput['Time(ns)'] += [frame_time]
                Goutput['Frame'] += [frame]
                Goutput['Nnative'] += [Nnative]
                for key, value in count_info.items():
                    Goutput[key] += [value]

        # Step 4a: save the dataframe that contains the detailed information about entanglementts
        print(f'Step 4: Save both the detailed entanglment info file .EntInfo and the stats summary file .G')
        EntOutput = pd.DataFrame(EntOutput)
        print(f'EntOutput:\n{EntOutput.to_string()}')
        Entoutfile = os.path.join(self.Gpath, f'{self.outname}.EntInfo')
        EntOutput.to_csv(Entoutfile, index=False)
        print(f'SAVED: {Entoutfile}')
        logging.info(f'SAVED: {Entoutfile}')

        # Step 4b: save the dataframe that contains the high level stats information about entanglement changes in each frame
        # should have 1 line per frame of the trajectory
        Goutput = pd.DataFrame(Goutput)
        #Goutput['G'] /= Nnative
        print(f'Goutput:\n{Goutput.to_string()}')
        Goutfile = os.path.join(self.Gpath, f'{self.outname}.G')
        Goutput.to_csv(Goutfile, index=False)
        print(f'SAVED: {Goutfile}')
        logging.info(f'SAVED: {Goutfile}')
    #######################################################################################

    #######################################################################################
    def GaussLink(self, coor, contact_mask=[]):
        """
        Optimized function to calculate the Gaussian Linking values for each native contact in a structure
        """
        Nterm_thresh = 5
        Cterm_thresh = 5

        # make native contact contact map
        dist_matrix = squareform(pdist(coor))
        native_cmap = np.triu(np.where(dist_matrix <= 8, 1, 0), k=4)
        nc_indexs = np.stack(np.nonzero(native_cmap)).T

        # If a contact mask list is specified then only grab those contacts in that list to analyze
        # will speed up things when computing over trajectory frames as i only care about native contcats
        if contact_mask:
            #print(f'contact_mask found and will be used')
            masked_nc_indexs = []
            for (i,j) in nc_indexs:
                if [i, j] in contact_mask:
                    masked_nc_indexs += [[i, j]]
            nc_indexs = np.asarray(masked_nc_indexs)

        l = len(dist_matrix)
        range_l = np.arange(0, l - 1)
        range_next_l = np.arange(1, l)

        coor = coor.astype(np.float32)
        R = 0.5 * (coor[range_l] + coor[range_next_l])
        dR = coor[range_next_l] - coor[range_l]

        # Efficient calculation of dR_cross using broadcasting
        dR_cross = np.cross(dR[:, np.newaxis], dR)

        # Efficient Runit matrix calculation using broadcasting
        diff = R[:, np.newaxis] - R
        norms = np.linalg.norm(diff, axis=2, keepdims=True)
        
        # Avoid division by zero by replacing zeros in norms with a small value
        norms = np.where(norms == 0, 1e-8, norms)  # Replace 0 with small value to avoid division by zero
        norms = norms ** 3

        Runit = diff / norms

        # Use broadcasting for dot product calculation
        dot_matrix = np.einsum('ijk,ijk->ij', Runit, dR_cross)

        # Initialize nc_gdict and list to store entangled native contacts
        nc_gdict = {}
        ent_nc = []
        
        for i, j in nc_indexs:
            loop_range = np.arange(i, j)
            nterm_range = np.arange(Nterm_thresh, i - 5)
            cterm_range = np.arange(j + 6, l - (Cterm_thresh + 1))

            gn_pairs_array = np.array(np.meshgrid(nterm_range, loop_range)).T.reshape(-1, 2)
            gc_pairs_array = np.array(np.meshgrid(loop_range, cterm_range)).T.reshape(-1, 2)

            # Calculate gn and gc values
            gn_val = np.sum(dot_matrix[gn_pairs_array[:, 0], gn_pairs_array[:, 1]]) / (4.0 * np.pi) if gn_pairs_array.size else 0
            gc_val = np.sum(dot_matrix[gc_pairs_array[:, 0], gc_pairs_array[:, 1]]) / (4.0 * np.pi) if gc_pairs_array.size else 0

            # Store native contacts that exceed the threshold
            if np.abs(gn_val) >= 0.6 or np.abs(gc_val) >= 0.6:
                ent_nc.append((i, j))
            
            nc_gdict[(i, j)] = {
                'gn': round(gn_val, 5), 'gc': round(gc_val, 5),
                'Gn': custom_round(gn_val), 'Gc': custom_round(gc_val),
                'crossingsN': '', 'crossingsC': ''
            }

        # Get crossings for native contacts
        res = lasso_type(coor.tolist(), loop_indices=ent_nc, more_info=True, min_dist=[10, 6, 5])
        for nc, top_info in res.items():
            nc_gdict[nc]['crossingsN'], nc_gdict[nc]['crossingsC'] = ','.join(top_info['crossingsN']), ','.join(top_info['crossingsC'])

        return nc_gdict
    #######################################################################################

    #######################################################################################
    def GetLinkChanges(self, ref_data, frame_data, frame_time, frame, Nnative):
        """
        Calculate the change types in this frame. Each terminal tail (N or C or both) can have one of these change types
        change # any change in the linking number status
        Loss # a loss in the linking number status but no change in chirality (i.e. 2 to 1 | 1 to 0)
        LossChiral # a loss in the linking number status with a change in chirality (i.e. 2 to -1 | -2 to 1)
        Gain # a gain in the linking number status but no change in chirality (i.e. 1 to 2 | 0 to 1)
        GainChiral # a gain in the linking number status with a change in chirality (i.e. 1 to -2 | -1 to 2)
        PureChiral # a pure switch in chirality of the linking number (i.e. -1 to 1 | 2 to -2)
        """
        # this dictionary will track all the detailed information about changes of entanglement observed at each frame
        changes = {'Time(ns)':[], 'Frame':[], 'i':[], 'j':[], 'gn':[], 'Gn':[], 'crossingsN':[], 'gc':[], 'Gc':[], 'crossingsC':[], 'NchangeType':[], 'CchangeType':[]}

        #QC there should be no non-native contacts in the frame_data
        for nc in frame_data:
            if nc not in ref_data:
                raise ValueError(f'There should be no non-native contacts in frame_data but {nc}')

        # This dictionary will track the counts of the changes to output a less detailed text file
        counts = {'N-Loss':0, 'N-LossChiral':0, 'N-Gain':0, 'N-GainChiral':0, 'N-PureChiral':0, 
                    'C-Loss':0, 'C-LossChiral':0, 'C-Gain':0, 'C-GainChiral':0, 'C-PureChiral':0, 'G':0}

        # loop through frame native contacts and look for changes
        for nc, frame_ent_info in frame_data.items():
            ref_Gn = ref_data[nc]['Gn']
            ref_Gc = ref_data[nc]['Gc']
            
            frame_Gn = frame_ent_info['Gn']
            frame_Gc = frame_ent_info['Gc']
            
            # determine N change type
            Nloss = abs(frame_Gn) < abs(ref_Gn)
            Ngain = abs(frame_Gn) > abs(ref_Gn)
            NnoChange = abs(frame_Gn) == abs(ref_Gn)
            NchiralChange = np.sign(frame_Gn) != np.sign(ref_Gn)
            
            if Nloss and not NchiralChange:
                NchangeType = 'Loss' # Loss of entanglement but no chiral switch
                counts['N-Loss'] += 1

            elif Nloss and NchiralChange:
                NchangeType = 'LossChiral' # Loss of entanglement with chiral switch
                counts['N-LossChiral'] += 1

            elif Ngain and not NchiralChange:
                NchangeType = 'Gain' # Gain of entanglement but no chiral switch
                counts['N-Gain'] += 1

            elif Ngain and NchiralChange:
                NchangeType = 'GainChiral' # Gain of entanglement with chiral switch
                counts['N-GainChiral'] += 1

            elif NnoChange and NchiralChange:
                NchangeType = 'PureChiral' # No change in abs(linking value) but a pure chiral switch
                counts['N-PureChiral'] += 1
                
            elif NnoChange and not NchiralChange:
                NchangeType = 'NoChange' # No change in abs(linking value) and no change in chirality

            Closs = abs(frame_Gc) < abs(ref_Gc)
            Cgain = abs(frame_Gc) > abs(ref_Gc)
            CnoChange = abs(frame_Gc) == abs(ref_Gc)
            CchiralChange = np.sign(frame_Gc) != np.sign(ref_Gc)

            if Closs and not CchiralChange:
                CchangeType = 'Loss' # Loss of entanglement but no chiral switch
                counts['C-Loss'] += 1

            elif Closs and CchiralChange:
                CchangeType = 'LossChiral' # Loss of entanglement with chiral switch
                counts['C-LossChiral'] += 1

            elif Cgain and not CchiralChange:
                CchangeType = 'Gain' # Gain of entanglement but no chiral switch
                counts['C-Gain'] += 1

            elif Cgain and CchiralChange:
                CchangeType = 'GainChiral' # Gain of entanglement with chiral switch
                counts['C-GainChiral'] += 1

            elif CnoChange and CchiralChange:
                CchangeType = 'PureChiral' # No change in abs(linking value) but a pure chiral switch
                counts['C-PureChiral'] += 1

            elif CnoChange and not CchiralChange:
                CchangeType = 'NoChange' # No change in abs(linking value) and no change in chirality

            # count if any change occured to track G metric
            if NchangeType != 'NoChange' or CchangeType != 'NoChange':
                counts['G'] += 1

            #print(nc, ref_data[nc], frame_ent_info, NchangeType, CchangeType)
            
            i,j = nc[0], nc[1]
            changes['Time(ns)'] += [frame_time]
            changes['Frame'] += [frame]
            changes['i'] += [i]
            changes['j'] += [j]
            changes['gn'] += [frame_ent_info['gn']]
            changes['Gn'] += [frame_ent_info['Gn']]
            changes['crossingsN'] += [frame_ent_info['crossingsN']]
            changes['gc'] += [frame_ent_info['gc']]
            changes['Gc'] += [frame_ent_info['Gc']]
            changes['crossingsC'] += [frame_ent_info['crossingsC']]
            changes['NchangeType'] += [NchangeType]
            changes['CchangeType'] += [CchangeType]

        counts["G"] /= Nnative
        logging.info(f'Time(ns): {frame_time} Frame: {frame} with G: {counts}')
        return changes, counts
    #######################################################################################


## Round GaussLink values
def custom_round(number):
    if number >= 0:
        # For positive numbers, round up if fractional part >= 0.6
        return np.ceil(number) if number % 1 >= 0.6 else np.floor(number)
    else:
        # For negative numbers, round down if the absolute fractional part >= 0.6
        # need to take the abs of the number first else the modulus does work right for negative numbers?
        return np.floor(number) if abs(abs(number) % 1) >= 0.6 else np.ceil(number)

def process_frame(frame_data):
    frame_coor, nc_list, ref_nc_gdict, time, frame, GaussLink, GetLinkChanges, Nnative = frame_data
    # Call GaussLink function
    frame_nc_gdict = GaussLink(frame_coor, contact_mask=nc_list)
    # Call GetLinkChanges function
    change_info, count_info = GetLinkChanges(ref_nc_gdict, frame_nc_gdict, time, frame, Nnative)
    return change_info, count_info
 
############## MAIN #################
def main():
    
    script_name = f'GQ'
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("--outname", type=str, required=True, help="base name for output files")
    parser.add_argument("--psf", type=str, required=True, help="Path to CA protein structure file")
    parser.add_argument("--cor", type=str, required=True, help="Path to CA native coordinates file")
    parser.add_argument("--dcd", type=str, required=True, help="Path to trajectory to analyze")
    parser.add_argument("--sec_elements", type=str, required=True, help="Path to STRIDE secondary structure elements file")
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

    # Get the fraction of native contacts Q
    anal.Q()

    # Get the Change in entanglement parameter G and its associated info
    anal.G()


if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    
print(f'NORMAL TERMINATION: {end_time - start_time}')
logging.info(f'NORMAL TERMINATION: {end_time - start_time}')