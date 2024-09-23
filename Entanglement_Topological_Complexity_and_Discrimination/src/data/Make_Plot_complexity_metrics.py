import sys,os
import numpy as np
#python src/data/Plot_complexity_metrics.py 
#-m ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Entanglement_Topological_Complexity_and_Discrimination/EXP/merged_stats_uent_data_\*_0.csv 
#-l ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Entanglement_Topological_Complexity_and_Discrimination/logs/Plot_complexity_metrics_EXP_0.log 
#-o ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Entanglement_Topological_Complexity_and_Discrimination/Plots/EXP/ -s mean

path_to_slug = '../../../git_slugs/Failure-to-Form_Native_Entanglements_slug'
for set_type in ['EXP', 'AF']:
    for spa in np.arange(0,100,10):
        for metric in ['mean', 'median']:

            script = f'python src/data/Plot_complexity_metrics.py'
            outpath = f'-o {path_to_slug}/Entanglement_Topological_Complexity_and_Discrimination/Plots/{set_type}/'
            log = f'-l {path_to_slug}/Entanglement_Topological_Complexity_and_Discrimination/logs/Plot_complexity_metrics_{set_type}_Rall_spa{spa}_LiPMScov50.log'
            merged = f'-m ../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Entanglement_Topological_Complexity_and_Discrimination/{set_type}/merged_stats_uent_data_\*_{spa}.csv'
            metric = f'-s {metric} -t {set_type}_{spa}'


            cmd = ' '.join([script, outpath, merged, metric, log])
            print(cmd)

        



