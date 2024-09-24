import os
from scipy.stats import false_discovery_control
import glob
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

class RegressionDataLoader:
    """
    This class is responsible for loading and consolidating regression data from multiple files.
    """
    def __init__(self, file_pattern):
        """
        Initializes the data loader with a glob pattern for file matching.

        :param file_pattern: A glob pattern to find input files.
        """
        self.file_pattern = file_pattern
        self.data = pd.DataFrame()

    def load_data(self):
        """
        Loads and concatenates regression data from files matching the input pattern.
        """
        files = glob.glob(self.file_pattern)
        for file in files:
            print(file)
            if self.data.empty:
                self.data = pd.read_csv(file, sep='|')
            else:
                additional_data = pd.read_csv(file, sep='|')
                self.data = pd.concat([self.data, additional_data], ignore_index=True)
        print(f'Data loaded and consolidated. {len(files)} files processed.')

class RegressionPlotter:
    """
    This class is responsible for plotting regression data.
    """
    def __init__(self, data):
        """
        Initializes the plotter with regression data.

        :param data: A DataFrame containing regression data.
        """
        self.data = data

    def plot_data(self, outpath, tag, variable):
        """
        Generates and saves plots for the specified regression variable.

        :param outpath: Directory where the plots will be saved.
        :param tag: Tag to include in the filename for identifying the output.
        :param variable: The variable on which regression is performed.
        """
        buff_tag = {'C': 'cyto-serum', 'CD': '+DnaK', 'CG': '+GroEL'}

        for cov in [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]:
            outfile = os.path.join(outpath, f'{tag}_binomial_regression_results_var-{variable}_LiPMScov{cov}.png')
            outfile_csv = os.path.join(outpath, f'{tag}_binomial_regression_results_var-{variable}_LiPMScov{cov}.csv')
            fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12, 4))

            save_df = {'buff':[], 'spa':[], 'cov':[], 'OR':[], 'OR_lb':[], 'OR_ub':[], 'pvalues':[], 'n':[]}
            for buff_i, buff in enumerate(['C']):
                plot_df = self.data[(self.data['buff'] == buff) & (self.data['var'] == variable) & (self.data['cov'] == cov)]
                print(plot_df.to_string())

                x = plot_df['spa']
                OR = plot_df['OR']
                OR_lb_delta = OR - plot_df['[0.025']
                OR_ub_delta = plot_df['0.975]'] - OR
                pvalues = plot_df['P>|z|'].replace(0, 0.00000000000001)
                ## Apply FDR
                #pvalues = false_discovery_control(pvalues.values)
                n = plot_df['n']

                axes[0].errorbar(x, OR, yerr=[OR_lb_delta, OR_ub_delta], label='OR', marker='o', ls='none', fillstyle='none', color='k', capsize=3)
                axes[0].set_title(buff_tag[buff])
                axes[0].axhline(y=1, color='red', linestyle='--')

                axes[1].plot(x, pvalues, label='pvalue', marker='o', ls='none', fillstyle='none', color='k')
                axes[1].axhline(y=0.05, color='red', linestyle='--')

                axes[2].plot(x, n, label='n', marker='o', ls='none', fillstyle='none', color='k')
                self._configure_axes(axes)

                #add data to save_df
                save_df['buff'] += [buff]*len(x)
                save_df['spa'] += [x for x in x]
                save_df['cov'] += [c for c in plot_df['cov']]
                save_df['OR'] += [o for o in OR]
                save_df['OR_lb'] += [o for o in plot_df['[0.025']]
                save_df['OR_ub'] += [o for o in plot_df['0.975]']]
                save_df['pvalues'] += [p for p in pvalues]
                save_df['n'] += [n for n in n]

                
            save_df = pd.DataFrame(save_df)
            print(f'save_df:\n{save_df}')
            save_df.to_csv(outfile_csv, index=False)
            print(f'SAVED: {outfile_csv}')

            plt.suptitle(f'Binomial Regression | {tag} | var {variable} | LiPMScov{cov}')
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
                ax.set_ylim(0.5, 2)
                ax.yaxis.set_major_locator(MultipleLocator(1))
                ax.yaxis.set_minor_locator(MultipleLocator(0.5))
                ax.set_xlabel('SPA threshold')
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
                ax.set_ylim(10, 1000)
                ax.grid(True, which='both', linestyle='--', linewidth=0.5)




def main():
    parser = argparse.ArgumentParser(description="Process regression data and generate plots.")
    parser.add_argument("-f", "--inp_files", type=str, required=True, help="Input file pattern for regression data.")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory.")
    parser.add_argument("-t", "--tag", type=str, required=True, help="Tag for output filenames.")
    parser.add_argument("-r", "--regression_var", type=str, required=True, help="regression variable you wish to plot")
    args = parser.parse_args()
    print(args)

    if not os.path.exists(args.outpath):
        os.makedirs(args.outpath)
        print(f'MADE: {args.outpath}')

    loader = RegressionDataLoader(args.inp_files)
    loader.load_data()

    plotter = RegressionPlotter(loader.data)
    plotter.plot_data(args.outpath, args.tag, args.regression_var)

    print('NORMAL TERMINATION')

if __name__ == "__main__":
    main()

