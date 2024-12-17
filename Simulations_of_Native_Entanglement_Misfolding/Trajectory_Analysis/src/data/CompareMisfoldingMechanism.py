import pandas as pd
import glob
from scipy.stats import bootstrap, permutation_test
import numpy as np
import time
import argparse
import sys, os, logging
from multiprocessing import Pool, cpu_count

class Analysis:
    """
    A class to perform analysis on Molecular dynamics trajecotry data and plot the results. 
    """
    #######################################################################################
    def __init__(self, args):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
        ("--outpath", type=str, required=True, help="Path to output directory")
        ("--candidates", type=str, required=True, help="A file containing two columns. The candidate tag and the groupID")
        ("--CollectedOPpath", type=str, required=True, help="path to the CollectAndProcessesOP DATA/ folder")
        ("--toplevel", type=str, required=True, help="file containing relative paths to either native state GQ files or the MSM file for various proteins")
        ("--outname", type=str, required=True, help="base name for output files")
        ("--setID", type=int, required=True, help="setID to analyze (2 or 3)")
        """

        # parse the parameters 
        self.setID = args.setID
        self.candidates = pd.read_csv(args.candidates)
        self.candidates = self.candidates[self.candidates['set'] == self.setID]
        logging.info(f'candidates:\n{self.candidates}')
        print(f'candidates:\n{self.candidates}')

        self.outpath = args.outpath
        logging.info(f'outpath: {self.outpath}')
        print(f'outpath: {self.outpath}')

        self.outname = args.outname
        logging.info(f'outname: {self.outname}')
        print(f'outname: {self.outname}')

        self.CollectedOPpath = args.CollectedOPpath
        logging.info(f'CollectedOPpath: {self.CollectedOPpath}')
        print(f'CollectedOPpath: {self.CollectedOPpath}')

        self.toplevel = args.toplevel
        logging.info(f'toplevel: {self.toplevel}')
        print(f'toplevel: {self.toplevel}')

        ## make Plots dir
        self.plot_path = os.path.join(self.outpath, 'Plots')
        if not os.path.exists(self.plot_path):
            os.makedirs(self.plot_path)
            print(f'Made directory: {self.plot_path}')  
        print(f'plot_path: {self.plot_path}') 
        
        ## make data dir
        self.data_path = os.path.join(self.outpath, 'DATA')
        if not os.path.exists(self.data_path):
            os.makedirs(self.data_path)
            print(f'Made directory: {self.data_path}')
        print(f'data_path: {self.data_path}')    
    #######################################################################################

    #######################################################################################
    def load_OP(self,):
        """
        Load the following files
        1. Load the mirror file containing those hand annotated mirror trajectories
        2. Load the collected GQK dataframe
        3. 
        """

        ###############################################################################################
        ## Load the collected GQ data
        #set2_OPfile: ../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CollectAndProcessOP/setID2/DATA/Quench_Collected_GQK.csv
        GQdata_f = os.path.join(self.CollectedOPpath, f'setID3/DATA/Quench_Collected_GQK.csv')
        print(GQdata_f)
        GQdata = pd.read_csv(GQdata_f)
        GQdata = GQdata[GQdata['gene'].isin(self.candidates['gene'])]
        GQdata = GQdata[GQdata['Mirror'] == False]
        print(GQdata)
        print(GQdata['gene'].value_counts())
        ###############################################################################################


        ###############################################################################################
        ## Get the thresholding metrics for later 
        self.threshold_metrics = {'gene':[], 'traj':[], '<Q>':[], 'MisfoldingProp':[], 'num_frames':[], 'num_misfolded_frames':[]}
        for gene, gene_df in GQdata.groupby('gene'):
            for traj, traj_df in gene_df.groupby('traj'):
                #print(traj_df)
                avgQ = np.mean(traj_df['Q'].values)
                MisfoldingProp = 1 - np.mean(traj_df['NativeByRef'].values)
                print(f'{gene} {traj} <Q> = {avgQ} and MisfoldingProp = {MisfoldingProp}')
                self.threshold_metrics['gene'] += [gene]
                self.threshold_metrics['traj'] += [traj]
                self.threshold_metrics['<Q>'] += [avgQ]
                self.threshold_metrics['MisfoldingProp'] += [MisfoldingProp]
                self.threshold_metrics['num_frames'] += [len(traj_df)]
                self.threshold_metrics['num_misfolded_frames'] += [np.sum(~traj_df['NativeByRef'])]

        threshold_metrics_outfile = os.path.join(self.data_path, f'Combined_and_Processed_threshold_metrics_setID{self.setID}.csv')
        self.threshold_metrics = pd.DataFrame(self.threshold_metrics)
        self.threshold_metrics = self.threshold_metrics.sort_values(by=['MisfoldingProp'])
        print(self.threshold_metrics.to_string())
        self.threshold_metrics.to_csv(threshold_metrics_outfile, index=False)
        print(f'SAVED: {threshold_metrics_outfile}')               
        ###############################################################################################

        ###############################################################################################
        ### Get the combine EntInfo file
        ## only get those frames that are from non Mirror traj and considered misfolded by Ref
        combined_EntInfo_f = os.path.join(self.CollectedOPpath, f'setID3/DATA/Collected_EntInfo.csv')
        print(combined_EntInfo_f)
        self.combined_EntInfo = pd.read_csv(combined_EntInfo_f, low_memory=False)
        print(self.combined_EntInfo)
        new_dfs = []
        for gene, gene_df in GQdata.groupby('gene'):
            for traj, traj_df in gene_df.groupby('traj'):
                #print(traj_df)
                misfolded_frames = GQdata[(GQdata['NativeByRef'] == False)]['frame'].values
                #print(f'misfolded_frames: {misfolded_frames} {misfolded_frames.shape}')

                EntInfo_df = self.combined_EntInfo[(self.combined_EntInfo['gene'] == gene) & (self.combined_EntInfo['traj'] == traj)]
                #EntInfo_df = EntInfo_df[EntInfo_df['Frame'].isin(misfolded_frames)]
                #print(EntInfo_df)
                new_dfs += [EntInfo_df]

        self.combined_EntInfo = pd.concat(new_dfs, ignore_index=True)
        print(f'combined_EntInfo:\n{self.combined_EntInfo}')
        ###############################################################################################

        self.combined_EntInfo['crossingsN'] = self.combined_EntInfo['crossingsN'].astype(str)
        self.combined_EntInfo['crossingsC'] = self.combined_EntInfo['crossingsC'].astype(str)
    #######################################################################################  

    #######################################################################################
    def Mechanism_stats(self,):

        PerTraj_Mechanism_summary_outfile = os.path.join(self.data_path, f'PerTraj_Mechanism_summary_setID{self.setID}.csv')
        if not os.path.exists(PerTraj_Mechanism_summary_outfile):
            # quick and dirty counts and stats for Loss and Gain fractions
            Loss_dist = []
            Gain_dist = []
            summary_df = {'gene':[], 'traj':[], 'num_frames':[], 'Loss':[], 'Gain':[], 'Total':[], 'LossOverlap':[], 'GainOverlap':[], 'Total_Overlap':[]}

            #for gene, gene_df in self.combined_EntInfo.groupby('gene'):
            #    for traj, traj_df in gene_df.groupby('traj'):
            #        for frame, frame_df in traj_df.groupby('Frame'):
            #            print(frame_df)
            #            Loss, Gain, LossOverlap, GainOverlap = find_overlap(frame_df)
            #            print(f'{gene} {traj} {frame} || Loss: {Loss} | Gain: {Gain} | LossOverlap: {LossOverlap} | GainOverlap: {GainOverlap}')
            #        quit()

            results = parallel_process_genes(self.combined_EntInfo, self.toplevel)
            for gene, traj, num_frames, traj_Loss, traj_Gain, Total, traj_LossOverlap, traj_GainOverlap, Total_Overlap in results:

                print(f'Gene: {gene} | Traj: {traj} | num_frames: {num_frames} | traj_Loss: {traj_Loss} | traj_Gain: {traj_Gain} | Total: {Total} | traj_LossOverlap: {traj_LossOverlap} | traj_GainOverlap: {traj_GainOverlap} | Total_Overlap: {Total_Overlap}')

                summary_df['gene'] += [gene]
                summary_df['traj'] += [traj]
                summary_df['num_frames'] += [num_frames]
                summary_df['Loss'] += [traj_Loss]
                summary_df['Gain'] += [traj_Gain]
                summary_df['Total'] += [Total]
                summary_df['LossOverlap'] += [traj_LossOverlap]
                summary_df['GainOverlap'] += [traj_GainOverlap]
                summary_df['Total_Overlap'] += [Total_Overlap]

            for k,v in summary_df.items():
                print(k, len(v))
            summary_df = pd.DataFrame(summary_df)
            print(summary_df.to_string())
            summary_df.to_csv(PerTraj_Mechanism_summary_outfile)
            print(f'SAVED: {PerTraj_Mechanism_summary_outfile}')
        else:
            summary_df = pd.read_csv(PerTraj_Mechanism_summary_outfile)
            print(f'LOADED: {PerTraj_Mechanism_summary_outfile}')
            print(summary_df.to_string())
        ##########################################################################################
 

        ##########################################################################################
        ## filter for those that meet the threshold for sturcutre content and misfolding propensity
        print(f'Filter for <Q> >= 0.6 and MisfoldingProp > 0.8')
        newdf = []
        for rowi, row in self.threshold_metrics.iterrows():
            gene, traj, avgQ, MisfoldingProp = row['gene'], row['traj'], row['<Q>'], row['MisfoldingProp']
            if avgQ >= 0.6:
                #if MisfoldingProp >= 0.8:
                if MisfoldingProp >= 0.8:
                    newdf += [summary_df[(summary_df['gene'] == gene) & (summary_df['traj'].astype(str) == str(traj))]]
        summary_df = pd.concat(newdf)
        print(summary_df.to_string())
        ##########################################################################################


        ##########################################################################################
        stats_df = {'type':[], 'mean':[], 'lb':[], 'ub':[], 'perm_stat':[], 'pvalue':[]}
        ## make the fraction of Loss and Gain columns
        summary_df['FracLoss'] = summary_df['Loss']/summary_df['Total']
        summary_df['FracGain'] = summary_df['Gain']/summary_df['Total']

        ## get the statistics for the Loss and Gain
        median, mean, std, lb, ub = get_stats(summary_df['FracLoss'].values)
        print(f'FracLoss stats: {mean} ({lb}, {ub})')
        stats_df['type'] += ['FracLoss']
        stats_df['mean'] += [mean]
        stats_df['lb'] += [lb]
        stats_df['ub'] += [ub]

        median, mean, std, lb, ub = get_stats(summary_df['FracGain'].values)
        print(f'FracGain stats: {mean} ({lb}, {ub})')
        stats_df['type'] += ['FracGain']
        stats_df['mean'] += [mean]
        stats_df['lb'] += [lb]
        stats_df['ub'] += [ub]

        perm_res = permutation_test((summary_df['FracLoss'].values, summary_df['FracGain'].values), statistic)
        pvalue = perm_res.pvalue
        perm_stat = perm_res.statistic
        print(f'Permutation test: stat = {perm_stat} with p-value {pvalue}')
        stats_df['perm_stat'] += [perm_stat, perm_stat]
        stats_df['pvalue'] += [pvalue, pvalue]


        ## Get the stats for the probaility of overlap and the conditional probability of overlap
        summary_df['FracOverlap'] = summary_df['Total_Overlap']/summary_df['Total']
        median, mean, std, lb, ub = get_stats(summary_df['FracOverlap'].values)
        print(f'FracOverlap stats: {mean} ({lb}, {ub})')
        stats_df['type'] += ['FracOverlap']
        stats_df['mean'] += [mean]
        stats_df['lb'] += [lb]
        stats_df['ub'] += [ub]
        stats_df['perm_stat'] += [np.nan]
        stats_df['pvalue'] += [np.nan]

        summary_df['CondLossOverlap'] = summary_df['LossOverlap']/summary_df['Loss']
        CondLossOverlap = summary_df['CondLossOverlap'].values
        median, mean, std, lb, ub = get_stats(CondLossOverlap[~np.isnan(CondLossOverlap)])
        print(f'CondLossOverlap stats: {mean} ({lb}, {ub})')
        stats_df['type'] += ['CondLossOverlap']
        stats_df['mean'] += [mean]
        stats_df['lb'] += [lb]
        stats_df['ub'] += [ub]
        stats_df['perm_stat'] += [np.nan]
        stats_df['pvalue'] += [np.nan]

        summary_df['CondGainOverlap'] = summary_df['GainOverlap']/summary_df['Gain']
        median, mean, std, lb, ub = get_stats(summary_df['CondGainOverlap'].values)
        print(f'CondGainOverlap stats: {mean} ({lb}, {ub})')
        stats_df['type'] += ['CondGainOverlap']
        stats_df['mean'] += [mean]
        stats_df['lb'] += [lb]
        stats_df['ub'] += [ub]
        stats_df['perm_stat'] += [np.nan]
        stats_df['pvalue'] += [np.nan]

        print(summary_df.to_string())

        stats_df = pd.DataFrame(stats_df)
        print(f'stats_df:\n{stats_df}')

        Mechanism_stats_summary_outfile = os.path.join(self.data_path, f'Mechanism_stats_summary_setID{self.setID}.csv')
        stats_df.to_csv(Mechanism_stats_summary_outfile)
        print(f'SAVED: {Mechanism_stats_summary_outfile}')
        #######################################################################################

#######################################################################################
def find_overlap(df):
    #print(df)
    #df['crossingsN'] = df['crossingsN'].astype(str)
    #df['crossingsC'] = df['crossingsC'].astype(str)

    ## if there is only a single unique ent then determine if its a loss or gain and return
    Loss = 0
    Gain = 0
    LossOverlap = 0
    GainOverlap = 0
    if len(df['cID'].unique()) == 1:

        if df['NchangeType'].str.contains('Loss').any():
            Loss += 1
        if df['CchangeType'].str.contains('Loss').any():
            Loss += 1
            
        if df['NchangeType'].str.contains('Gain').any():
            Gain += 1
        if df['CchangeType'].str.contains('Gain').any():
            Gain += 1
            

    else:
        ## else find overlaps
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
            
            cID_res[cID] = {'ijres':ijres, 'Nres':Nres, 'Ntype':NchangeType, 'Cres':Cres, 'Ctype':CchangeType}

        ################################
        # determine if there is overlap
        for cID, cID_info in cID_res.items():
            #print('\n', cID, cID_info)

            # count if the change is a loss or a gain
            if 'Loss' in cID_info['Ntype']:
                Loss += 1
            if 'Loss' in cID_info['Ctype']:
                Loss += 1
                
            if 'Gain' in cID_info['Ntype']:
                Gain += 1
            if 'Gain' in cID_info['Ctype']:
                Gain += 1

            ### Check if the N term has overlap with anything 
            if cID_info['Ntype'] != 'NoChange':

                for Comp_cID, Comp_cID_info in cID_res.items():

                    if Comp_cID == cID:
                        continue

                    # is there overlap
                    if cID_info[f'Nres'].intersection(Comp_cID_info[f'Nres']):
                        if ('Loss' in cID_info[f'Ntype'] and 'Gain' in Comp_cID_info[f'Ntype']):
                            LossOverlap += 1
                            break
                        
                        if ('Gain' in cID_info[f'Ntype'] and 'Loss' in Comp_cID_info[f'Ntype']):
                            GainOverlap += 1   
                            break


            if cID_info['Ctype'] != 'NoChange':

                for Comp_cID, Comp_cID_info in cID_res.items():

                    if Comp_cID == cID:
                        continue

                    # is there overlap
                    if cID_info[f'Cres'].intersection(Comp_cID_info[f'Cres']):
                        if ('Loss' in cID_info[f'Ctype'] and 'Gain' in Comp_cID_info[f'Ctype']):
                            LossOverlap += 1
                            break
                        
                        if ('Gain' in cID_info[f'Ctype'] and 'Loss' in Comp_cID_info[f'Ctype']):
                            GainOverlap += 1   
                            break

    return  Loss, Gain, LossOverlap, GainOverlap
#######################################################################################

#######################################################################################
def process_gene(args):
    gene, gene_df, toplevel = args
    results = []

    for traj, traj_df in gene_df.groupby('traj'):
        traj_Loss = 0
        traj_Gain = 0
        traj_LossOverlap = 0
        traj_GainOverlap = 0
        num_frames = 0

        for frame, frame_df in traj_df.groupby('Frame'):
            num_frames += 1

            ## determine the overlap stats for the frame
            Loss, Gain, LossOverlap, GainOverlap = find_overlap(frame_df)

            traj_Loss += Loss
            traj_Gain += Gain
            traj_LossOverlap += LossOverlap
            traj_GainOverlap += GainOverlap

        Total = traj_Loss + traj_Gain
        Total_Overlap = traj_LossOverlap + traj_GainOverlap
        if Total > 0:
            results.append((gene, traj, num_frames, traj_Loss, traj_Gain, Total, traj_LossOverlap, traj_GainOverlap, Total_Overlap))

    return results

    
####################################################
def parallel_process_genes(combined_EntInfo, toplevel):
    # Prepare arguments for parallel processing
    args = [(gene, gene_df, toplevel) for gene, gene_df in combined_EntInfo.groupby('gene')]

    # Use multiprocessing to parallelize the processing
    num_cpu = cpu_count()
    #num_cpu = 1
    print(f'cpu_count: {num_cpu}')
    with Pool(num_cpu) as pool:
        results = pool.map(process_gene, args)

    # Flatten the list of results
    flattened_results = [item for sublist in results for item in sublist]
    return flattened_results
#######################################################################################

#######################################################################################
def get_stats(arr):
    """
    Get the <> and 95% ci for a stats array
    """
    mean = np.mean(arr)
    std = np.std(arr)
    res = bootstrap((arr,), np.mean)
    lb = res.confidence_interval.low
    ub = res.confidence_interval.high
    median = np.median(arr)
    return (median, mean, std, lb, ub)
#######################################################################################

#######################################################################################
def statistic(x, y, axis):
    return np.mean(x, axis=axis) - np.mean(y, axis=axis)
#######################################################################################

############## MAIN #################
def main():
    
    script_name = f'CompareMisfoldingMechanism'
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("--candidates", type=str, required=True, help="A file containing two columns. The candidate tag and the groupID")
    parser.add_argument("--CollectedOPpath", type=str, required=True, help="path to the CollectAndProcessesOP DATA/ folder")
    parser.add_argument("--toplevel", type=str, required=True, help="file containing relative paths to either native state GQ files or the MSM file for various proteins")
    parser.add_argument("--outname", type=str, required=True, help="base name for output files")
    parser.add_argument("--setID", type=int, required=True, help="setID to analyze (2 or 3)")
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


    anal = Analysis(args)
    anal.load_OP()
    anal.Mechanism_stats()


if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()

print(f'NORMAL TERMINATION: {time.time() - start_time}')
logging.info(f'NORMAL TERMINATION: {time.time() - start_time}')