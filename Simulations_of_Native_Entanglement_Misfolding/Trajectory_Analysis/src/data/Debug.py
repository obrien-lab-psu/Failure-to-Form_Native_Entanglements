import pandas as pd
import glob
from scipy.stats import permutation_test, mode
import numpy as np
import time
import argparse
import sys, os, logging
from multiprocessing import Pool, cpu_count
import matplotlib.pyplot as plt
start_time = time.time()

combined_EntInfo_f = f'../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis/CollectAndProcessOP_v2.1/setID3/DATA/Collected_EntInfo.csv'
print(combined_EntInfo_f)
combined_EntInfo = pd.read_csv(combined_EntInfo_f, low_memory=False)
print(f'combined_EntInfo:\n{combined_EntInfo}')

gene = 'P0A6B4'
traj = 16

gene_combined_EntInfo = combined_EntInfo[(combined_EntInfo['gene'] == gene) & (combined_EntInfo['traj'] == traj)]
print(f'gene_combined_EntInfo:\n{gene_combined_EntInfo}')

## count the number of unique loss, unique gain, total changes
Loss = 0
Gain = 0
for frame, frame_df in gene_combined_EntInfo.groupby('Frame'):
    print(f'{"#"*100}\nFRAME: {frame}\n',frame_df.to_string())
    if frame == 24029:
        quit()
    frame_Loss = 0
    frame_Gain = 0
    for cID, cID_df in frame_df.groupby('cID'):
        print(cID_df.to_string())
        if cID_df['NchangeType'].str.contains('Loss').any():
            frame_Loss += 1
        if cID_df['CchangeType'].str.contains('Loss').any():
            frame_Loss += 1
        if cID_df['NchangeType'].str.contains('Gain').any():
            frame_Gain += 1
        if cID_df['CchangeType'].str.contains('Gain').any():
            frame_Gain += 1 
        print(f'frame_Loss: {frame_Loss}')
        print(f'frame_Gain: {frame_Gain}')
    Loss += frame_Loss
    Gain += frame_Gain
    print(f'frame_Loss: {frame_Loss}')
    print(f'frame_Gain: {frame_Gain}')
total = Loss + Gain
print(Loss, Gain, total)        

print(f'NORMAL TERMINATION: {time.time() - start_time}')
logging.info(f'NORMAL TERMINATION: {time.time() - start_time}')