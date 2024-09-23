import os
import glob
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.gridspec as gridspec


class DataLoader:
    """
    This class is responsible for loading and consolidating regression data from multiple files.
    """
    def __init__(self, file_pattern, tag, outpath):
        """
        Initializes the data loader with a glob pattern for file matching.

        :param file_pattern: A glob pattern to find input files.
        """
        self.file_pattern = file_pattern
        self.data = pd.DataFrame()
        self.tag = tag
        self.outpath = outpath

    def load_data(self):
        """
        Loads and concatenates regression data from files matching the input pattern.
        """
        #EnrichedContactSignificance_v2.0/EnrichedContactSignificanceOutput/Contingency_analysis_CD_90_EXP_CD_spa90_group7_ent_genes_numLC1.csv
        files = glob.glob(self.file_pattern)
        print(files)
        dfs = []
        for file in files:
            print(file)
            df = pd.read_csv(file, sep='|')
            if 'EXP' in file:
                df['label'] = 'EXP'
            if 'AF' in file:
                df['label'] = 'AF'
            dfs += [df]
        self.data = pd.concat(dfs)
        print(f'Data loaded and consolidated. {len(files)} files processed.')
        print(self.data)
        outfile_csv = os.path.join(self.outpath, f'Assoc_Client_n_Essential_plot_data_{self.tag}.csv')
        self.data.to_csv(outfile_csv, sep='|', index=False)
        print(f'SAVED: {outfile_csv}')

class Plotter:
    """
    This class is responsible for plotting regression data.
    """
    def __init__(self, data, tag, outpath):
        """
        Initializes the plotter with regression data.

        :param data: A DataFrame containing regression data.
        """
        self.data = data
        print(self.data)
        self.tag = tag
        self.outpath = outpath

    def plot_data(self,):
        """
        Generates and saves plots for the specified regression variable.

        :param outpath: Directory where the plots will be saved.
        :param tag: Tag to include in the filename for identifying the output.
        :param variable: The variable on which regression is performed.
        """
        buff_tag = {'C': 'cyto-serum', 'CD': '+DnaK', 'CG': '+GroEL', 'Any':'All Rep genes'}
        outfile = os.path.join(self.outpath, f'Essentiality_and_client_Fisher_results_{self.tag}.png')

        fig, axes = plt.subplots(1, 3, figsize=(12, 4))

        for label, label_df in self.data.groupby('label'):
            for buff_i, buff in enumerate(['C']):
                buff_df = label_df[label_df['buff'] == buff]
                print(buff_i, buff,  buff_df.to_string())

                x = buff_df['spa']
                OR = buff_df['OR']
                OR_lb = buff_df['OR_lb']
                OR_lb_delta = OR - OR_lb
                OR_ub = buff_df['OR_ub']
                OR_ub_delta = OR_ub - OR
                pvalues = buff_df['pvalue']
                n = buff_df['n']


                axes[0].errorbar(x, OR, yerr=[OR_lb_delta, OR_ub_delta], marker='o', ls='none', fillstyle='none', capsize=3, label=label)
                axes[0].set_title(buff_tag[buff])
                axes[0].axhline(y=1, color='red', linestyle='--')

                axes[1].plot(x, pvalues, label=label, marker='o', ls='none', fillstyle='none')
                axes[1].axhline(y=0.05, color='red', linestyle='--')

                axes[2].plot(x, n, label=label, marker='o', ls='none', fillstyle='none')
                axes[2].legend(loc='upper left', bbox_to_anchor=(1,1))
                self._configure_axes(axes)

        plt.suptitle(f'{self.tag}', fontsize=7)
        plt.tight_layout()
        #plt.show()
        plt.savefig(outfile)
        plt.close()
        print(f'SAVED: {outfile}')

    def _configure_axes(self, axes):
        """
        Configures the axes properties for the plot.

        :param axes: The axes array from the subplot.
        :param buff_i: Index of the current buffer in enumeration.
        """
        OR_ticks = [0.5, 1.0, 1.5, 2.0, 2.5]
        n_ticks = [0, 250, 500, 750, 1000]
        for i in range(3):
            ax = axes[i]
            if i == 0:
                ax.set_ylabel('Odds Ratio')
                #ax.set_ylim(0, self.max_y + 0.1)
                #ax.set_yscale('log')
                ax.set_ylim(0.1, 3)
                ax.set_xlim(-10, 100)
                ax.set_xlabel('SPA threshold')
                #ax.yaxis.set_major_locator(MultipleLocator(1))
                #ax.yaxis.set_minor_locator(MultipleLocator(0.5))
                ax.grid(True, which='both', linestyle='--', linewidth=0.5)

            elif i == 1:
                ax.set_yscale('log')
                ax.set_ylabel('pvalue')
                ax.set_xlabel('SPA threshold')
                ax.grid(True, which='both', linestyle='--', linewidth=0.5)

            elif i == 2:
                ax.set_ylabel('n')
                ax.set_xlabel('SPA threshold')
                ax.set_yscale('log')
                ax.set_ylim(1, 1300)
                ax.grid(True, which='both', linestyle='--', linewidth=0.5)
         


def main():
    parser = argparse.ArgumentParser(description="Process regression data and generate plots.")
    parser.add_argument("-f", "--inp_files", type=str, required=True, help="Input file pattern for regression data.")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory.")
    parser.add_argument("-t", "--tag", type=str, required=True, help="tag for final output image")
    args = parser.parse_args()

    if not os.path.exists(args.outpath):
        os.makedirs(args.outpath)
        print(f'MADE: {args.outpath}')

    loader = DataLoader(args.inp_files, args.tag, args.outpath)
    loader.load_data()

    plotter = Plotter(loader.data, args.tag, args.outpath)
    plotter.plot_data()

    print('NORMAL TERMINATION')

if __name__ == "__main__":
    main()

