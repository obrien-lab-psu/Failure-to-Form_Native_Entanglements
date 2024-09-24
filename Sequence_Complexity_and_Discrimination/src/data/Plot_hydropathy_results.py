import os
import numpy as np
import glob
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import MaxNLocator

class DataLoader:
    """
    This class is responsible for loading and consolidating regression data from multiple files.
    """
    def __init__(self, file_pattern, outpath):
        """
        Initializes the data loader with a glob pattern for file matching.

        :param file_pattern: A glob pattern to find input files.
        """
        self.file_pattern = file_pattern
        self.data = pd.DataFrame()
        self.outpath = outpath

    def load_data(self):
        """
        Loads and concatenates regression data from files matching the input pattern.
        """
        files = glob.glob(os.path.join(self.file_pattern, '*'))
        print(f'Files found: {len(files)}')
        #compare_hydropathy_FLiPPR_SPAwCOV_PROD/EXP/CD_0/HydropathyAnalyzerOutpath/stats.csv
        #compare_hydropathy_outpath/CD_0/HydropathyAnalyzerOutpath/stats.csv
        for file in files:
            print(file)
            temp_df = pd.read_csv(file, sep='|')
            buff = file.split('/')[-1].split('_')[1]
            spa = file.split('/')[-1].split('_')[2]
            spa = spa.replace('spa','')
            print(file, buff, spa)

            temp_df['buff'] = buff
            temp_df['spa'] = spa

            if self.data.empty:
                self.data = temp_df
            else:
                additional_data = pd.read_csv(file, sep='|')
                self.data = pd.concat([self.data, temp_df], ignore_index=True)
        print(f'Data loaded and consolidated. {len(files)} files processed.')

        ## Save total loaded df used to plot data
        outfile = os.path.join(self.outpath, f'hydropathy_results.csv')
        self.data.to_csv(outfile, index=False)
        print(f'SAVED: {outfile}')

class Plotter:
    """
    This class is responsible for plotting regression data.
    """
    def __init__(self, data):
        """
        Initializes the plotter with regression data.

        :param data: A DataFrame containing regression data.
        """
        self.data = data
        print(f'Plotter_data:\n{self.data}')

    def plot_data(self, outpath):
        """
        Generates and saves plots for the specified regression variable.

        :param outpath: Directory where the plots will be saved.
        :param tag: Tag to include in the filename for identifying the output.
        :param variable: The variable on which regression is performed.
        """
        outfile = os.path.join(outpath, f'hydropathy_results.png')
        fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(12, 8))
        for Hclass_i, Hclass in enumerate(['SHphob', 'WHphob', 'Hphil', 'Arom']):

            Hclass_df = self.data[self.data['class'] == Hclass]
            min_y = np.min(Hclass_df[['lower_ci_Ess', 'lower_ci_NonEss']].values) - 0.02
            max_y = np.max(Hclass_df[['upper_ci_Ess', 'upper_ci_NonEss']].values) + 0.02

            rung = 0.025
            min_y_floor = np.floor(min_y / rung) * rung
            max_y_ceil = np.ceil(max_y / rung) * rung
            print(f'Hclass: {Hclass} | ({min_y}, {max_y}) | ({min_y_floor}, {max_y_ceil})')
            
            for buff_i, buff in enumerate(['C', 'CD', 'CG']):
                plot_df = self.data[(self.data['buff'] == buff) & (self.data['class'] == Hclass)]
                plot_df['upper_ci_max'] = plot_df[[f'upper_ci_Ess', f'upper_ci_NonEss']].max(axis=1)
                print(plot_df)

                for tag, color in {'Ess':'blue', 'NonEss':'red'}.items():
                    x = plot_df['spa']
                    y = plot_df[f'mean_{tag}']
                    y_lb_delta = y - plot_df[f'lower_ci_{tag}']
                    y_ub_delta = plot_df[f'upper_ci_{tag}'] - y

                    buff_tag = {'C': 'cyto-serum', 'CD': 'cyto-serum + DnaK', 'CG': 'cyto-serum + GroEL'}[buff]
                    axes[Hclass_i, buff_i].errorbar(x, y, yerr=[y_lb_delta, y_ub_delta], label=tag, marker='o', ls='none', fillstyle='none', color=color, capsize=3)
                    axes[Hclass_i, buff_i].set_title(buff_tag)

                    if buff_i == 0:
                        axes[Hclass_i, buff_i].set_ylabel('fract loop forming NC')
                    if Hclass_i == 3:
                        axes[Hclass_i, buff_i].set_xlabel('SPA threshold')

                    axes[Hclass_i, buff_i].set_title(f'{buff_tag} | {Hclass}')
                    axes[Hclass_i, buff_i].set_ylim(min_y_floor, max_y_ceil)
                    axes[Hclass_i, buff_i].yaxis.set_major_locator(MultipleLocator(rung))

                    if buff_i == 2:
                        axes[Hclass_i, buff_i].legend(loc='upper left', bbox_to_anchor=(1, 1))


                for i, (xi, yi, upper, pval) in enumerate(zip(x, y, plot_df[f'upper_ci_max'], plot_df['pvalue'])):
                    if pval < 0.001:
                        axes[Hclass_i, buff_i].text(xi, upper + 0.01, '***', ha='center', color='k')
                    elif pval < 0.01:
                        axes[Hclass_i, buff_i].text(xi, upper + 0.01, '**', ha='center', color='k')
                    elif pval < 0.05:
                        axes[Hclass_i, buff_i].text(xi, upper + 0.01, '*', ha='center', color='k')
                    print(i, xi, yi, upper, pval)

        #plt.suptitle(f'')
        plt.tight_layout()
        #plt.show()
        plt.savefig(outfile)
        plt.close()
        print(f'SAVED: {outfile}')


def main():
    parser = argparse.ArgumentParser(description="Process regression data and generate plots.")
    parser.add_argument("-f", "--inp_files", type=str, required=True, help="Input file pattern for hydropathy data.")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory.")
    args = parser.parse_args()

    if not os.path.exists(args.outpath):
        os.makedirs(args.outpath)
        print(f'MADE: {args.outpath}')

    loader = DataLoader(args.inp_files, args.outpath)
    loader.load_data()

    plotter = Plotter(loader.data)
    plotter.plot_data(args.outpath)

    print('NORMAL TERMINATION')

if __name__ == "__main__":
    main()

