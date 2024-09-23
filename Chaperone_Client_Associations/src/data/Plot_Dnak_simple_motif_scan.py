import os
import ast
import glob
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.stats import false_discovery_control

class DnaKBinderDataLoader:
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

        # Make outpath if it doesnt exists
        self.Outpath = outpath
        if not os.path.exists(self.Outpath):
            os.makedirs(self.Outpath)
            print(f'Made directory: {self.Outpath}')

    def load_data(self, tag):
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
        self.data.to_csv(f'{self.Outpath}consolidated_Dnak_binding_data_{tag}.csv', index=False)
        print(f'SAVED: {self.Outpath}consolidated_Dnak_binding_data_{tag}.csv')
        self.data = self.data[self.data['OnlyEnt'] == True]

    def correct_pvalues(self):
        """
        for each D_type considered correct the pvalues for a single SPA threshold
        """
        dfs = []
        for spa, spa_df in self.data.groupby('spa'):
            for D_type, D_type_df in spa_df.groupby('D_type'):
                print(D_type_df)
                D_type_df['q_value'] = false_discovery_control(D_type_df['p_value'].values)
                print(D_type_df)
                dfs += [D_type_df]
        self.data = pd.concat(dfs)
        print(self.data)


class DnakBinderDataPlotter:
    """
    This class is responsible for plotting regression data.
    """
    def __init__(self, data, outpath):
        """
        Initializes the plotter with regression data.

        :param data: A DataFrame containing regression data.
        """
        self.data = data
        self.Outpath = outpath

    def plot_data(self, tag, variable):
        """
        Generates and saves plots for the specified regression variable.

        :param outpath: Directory where the plots will be saved.
        :param tag: Tag to include in the filename for identifying the output.
        :param variable: The variable on which regression is performed.
        """
        pool_tag = {'D_All':'Whole Protein', 'D_Thread':'Entanglement Threads', 'D_Loop':'Entanglement Loops'}
        motifs = ['Bukau4', 'Bukau4RF', 'Bukau4LF', 'Bukau5', 'Bukau5RF', 'Bukau5LF', 'Schymkowitz', 'Emperical4', 'Emperical4RF', 'Emperical4LF', 'Emperical5', 'Emperical5RF', 'Emperical5LF']
        for motif_i, motif in enumerate(motifs):

            #for EntFlag in [True, False]:
            for EntFlag in [True]:

                outfile = os.path.join(self.Outpath, f'{motif}_ENTonly-{EntFlag}_{tag}.png')
                print(f'outfile: {outfile}')
                fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(12, 3.5))

                for pool_i, pool in enumerate(['D_All', 'D_Loop', 'D_Thread']):
                    plot_df = self.data[(self.data['motif'] == motif) & (self.data['D_type'] == pool) & (self.data['OnlyEnt'] == EntFlag)]
                    print(motif, EntFlag, pool)
                    print(plot_df.to_string())

                    x = plot_df['spa']
                    Ess = plot_df['Ess_mean']
                    Ess_lb_delta = plot_df['Ess_mean'] - plot_df['Ess_lower_ci']
                    Ess_ub_delta = plot_df['Ess_upper_ci'] - plot_df['Ess_mean']
                    NonEss = plot_df['NonEss_mean']
                    NonEss_lb_delta = plot_df['NonEss_mean'] - plot_df['NonEss_lower_ci']
                    NonEss_ub_delta = plot_df['NonEss_upper_ci'] - plot_df['NonEss_mean']
                    pvalues = plot_df['q_value']
                    
                    axes[0, pool_i].errorbar(x, Ess, yerr=[Ess_lb_delta, Ess_ub_delta], label='Ess', marker='o', ls='none', fillstyle='none', color='blue', capsize=3)
                    axes[0, pool_i].errorbar(x, NonEss, yerr=[NonEss_lb_delta, NonEss_ub_delta], label='NonEss', marker='o', ls='none', fillstyle='none', color='red', capsize=3)
                    #axes[0, pool_i].axhline(y=1, color='red', linestyle='--')
                    axes[0, pool_i].set_title(pool_tag[pool])
                    axes[0, pool_i].set_ylabel(r'$\frac{\# \text{ Motif Hits}}{\text{protein size}}$')
                    axes[0, pool_i].set_xlabel('SPA threshold')
                    #axes[0, pool_i].yaxis.set_major_locator(MultipleLocator(0.00025))
                    #axes[0, pool_i].set_yscale('log')
                    #axes[0, pool_i].set_ylim(0.0001, 0.001)

                    if pool_i == 2:
                        axes[0, pool_i].legend(loc='upper left', bbox_to_anchor=(1, 1))


                    axes[1, pool_i].plot(x, pvalues, label='pvalue', marker='o', ls='none', fillstyle='none', color='k')
                    axes[1, pool_i].axhline(y=0.05, color='red', linestyle='--')
                    axes[1, pool_i].set_yscale('log')
                    axes[1, pool_i].set_ylabel('pvalue')
                    axes[1, pool_i].grid(True, which='both', linestyle='--', linewidth=0.5)
                    axes[1, pool_i].set_xlabel('SPA threshold')

                    
                plt.suptitle(f'{motif} DnaK binding motif scan (EntOnly-{EntFlag})')
                plt.tight_layout()
                #plt.show()
                #quit()
                plt.savefig(outfile)
                plt.close()
                print(f'SAVED: {outfile}')


def main():
    parser = argparse.ArgumentParser(description="Process regression data and generate plots.")
    parser.add_argument("-f", "--inp_files", type=str, required=True, help="Input file pattern for DnaK binder data")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory.")
    parser.add_argument("-t", "--tag", type=str, required=True, help="Tag for output filenames.")
    args = parser.parse_args()

    if not os.path.exists(args.outpath):
        os.makedirs(args.outpath)
        print(f'MADE: {args.outpath}')

    loader = DnaKBinderDataLoader(args.inp_files, args.outpath)
    loader.load_data(args.tag)
    loader.correct_pvalues()

    plotter = DnakBinderDataPlotter(loader.data, loader.Outpath)
    plotter.plot_data(args.tag, 'region')

    print('NORMAL TERMINATION')

if __name__ == "__main__":
    main()

