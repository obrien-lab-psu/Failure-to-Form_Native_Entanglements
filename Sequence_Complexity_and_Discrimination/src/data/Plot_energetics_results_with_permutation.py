import os
import numpy as np
import glob
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import MaxNLocator
import seaborn as sns

class Plotter:
    """
    This class is responsible for loading and consolidating regression data from multiple files.
    Then plotting the Enorm, deltaE files, deltaDeltaE files
    """
    #def __init__(self, outpath, Enorm_file, Ess_deltaE_file, NonEss_deltaE_file, deltaDeltaE_file, pvalue_file, tag):
    def __init__(self, **kwargs):
        """
        outpath = outpath, 
        Enorm_file = Enorm_file, 

        Ess_deltaE_file = Ess_deltaE_file, 
        Ess_deltaE_pvalues_file = Ess_deltaE_pvalues_file,

        NonEss_deltaE_file = NonEss_deltaE_file,
        NonEss_deltaE_pvalues_files = NonEss_deltaE_pvalues_files,

        deltaDeltaE_file = deltaDeltaE_file, 
        deltaDeltaE_pvalues_file = deltaDeltaE_pvalues_file,

        tag = tag)
        """
        self.data = {}
        ########################################
        # Get outpath
        self.__dict__.update(kwargs)
        print(self.__dict__)
        self.outpath = self.__dict__['outpath']
        print(f'self.outpath: {self.outpath}')

        ########################################
        # Load Enorm df into class
        self.Enorm_file = self.__dict__['Enorm_file']
        print(f'self.Enorm_file: {self.Enorm_file}')
        self.Enorm = pd.read_csv(self.Enorm_file, sep='|')
        self.Enorm.set_index('AA', inplace=True)
        print(f'self.Enorm:\n{self.Enorm.to_string()}')
        self.data['Enorm'] = {'df':self.Enorm, 'pvalues':None}

        ########################################
        # Load Ess deltaE and its pvalue files into class
        self.Ess_deltaE_file = self.__dict__['Ess_deltaE_file']
        print(f'self.Ess_deltaE_file: {self.Ess_deltaE_file}')
        self.Ess_deltaE = pd.read_csv(self.Ess_deltaE_file, sep='|')
        self.Ess_deltaE.set_index('AA', inplace=True)
        print(f'self.Ess_deltaE:\n{self.Ess_deltaE.to_string()}')

        self.Ess_deltaE_pvalues_file = self.__dict__['Ess_deltaE_pvalues_file']
        print(f'self.Ess_deltaE_pvalues_file: {self.Ess_deltaE_pvalues_file}')
        if os.path.exists(self.Ess_deltaE_pvalues_file):
            self.Ess_deltaE_pvalues = pd.read_csv(self.Ess_deltaE_pvalues_file, sep='|')
            print(self.Ess_deltaE_pvalues)
            self.Ess_deltaE_pvalues.set_index('AA', inplace=True)
            print(self.Ess_deltaE_pvalues)
            a = self.Ess_deltaE_pvalues.values
            symmetric_a = np.maximum(a, a.T)
            self.Ess_deltaE_pvalues = pd.DataFrame(symmetric_a, index=self.Ess_deltaE_pvalues.index, columns=self.Ess_deltaE_pvalues.columns)
            self.Ess_deltaE_pvalues.at['C', 'C'] = 1
            print(f'self.Ess_deltaE_pvalues:\n{self.Ess_deltaE_pvalues.to_string()}')
            self.data['Ess_deltaE'] = {'df':self.Ess_deltaE, 'pvalues':self.Ess_deltaE_pvalues}
        else:
            self.data['Ess_deltaE'] = {'df':self.Ess_deltaE, 'pvalues':None}

        ########################################
        # Load NonEss deltaE and its pvalue files into class
        self.NonEss_deltaE_file = self.__dict__['NonEss_deltaE_file']
        print(f'self.NonEss_deltaE_file: {self.NonEss_deltaE_file}')
        self.NonEss_deltaE = pd.read_csv(self.NonEss_deltaE_file, sep='|')
        self.NonEss_deltaE.set_index('AA', inplace=True)
        print(f'self.NonEss_deltaE:\n{self.NonEss_deltaE.to_string()}')

        self.NonEss_deltaE_pvalues_file = self.__dict__['NonEss_deltaE_pvalues_file']
        print(f'self.NonEss_deltaE_pvalues_file: {self.NonEss_deltaE_pvalues_file}')
        if os.path.exists(self.NonEss_deltaE_pvalues_file):
            self.NonEss_deltaE_pvalues = pd.read_csv(self.NonEss_deltaE_pvalues_file, sep='|')
            self.NonEss_deltaE_pvalues.set_index('AA', inplace=True)
            a = self.NonEss_deltaE_pvalues.values
            symmetric_a = np.maximum(a, a.T)
            self.NonEss_deltaE_pvalues = pd.DataFrame(symmetric_a, index=self.NonEss_deltaE_pvalues.index, columns=self.NonEss_deltaE_pvalues.columns)
            self.NonEss_deltaE_pvalues.at['C', 'C'] = 1
            print(f'self.NonEss_deltaE_pvalues:\n{self.NonEss_deltaE_pvalues.to_string()}')
            self.data['NonEss_deltaE'] = {'df':self.NonEss_deltaE, 'pvalues':self.NonEss_deltaE_pvalues}
        else:
            self.data['NonEss_deltaE'] = {'df':self.NonEss_deltaE, 'pvalues':None}


        ########################################
        # Load Ess deltaE and its pvalue files into class
        self.deltaDeltaE_file = self.__dict__['deltaDeltaE_file']
        print(f'self.deltaDeltaE_file: {self.deltaDeltaE_file}')
        self.deltaDeltaE = pd.read_csv(self.deltaDeltaE_file, sep='|')
        self.deltaDeltaE.set_index('AA', inplace=True)
        print(f'self.deltaDeltaE:\n{self.deltaDeltaE.to_string()}')

        self.deltaDeltaE_pvalues_file = self.__dict__['deltaDeltaE_pvalues_file']
        print(f'self.deltaDeltaE_pvalues_file: {self.deltaDeltaE_pvalues_file}')
        if os.path.exists(self.deltaDeltaE_pvalues_file):
            self.deltaDeltaE_pvalues = pd.read_csv(self.deltaDeltaE_pvalues_file, sep='|')
            self.deltaDeltaE_pvalues.set_index('AA', inplace=True)
            a = self.deltaDeltaE_pvalues.values
            symmetric_a = np.maximum(a, a.T)
            self.deltaDeltaE_pvalues = pd.DataFrame(symmetric_a, index=self.deltaDeltaE_pvalues.index, columns=self.deltaDeltaE_pvalues.columns)
            self.deltaDeltaE_pvalues.at['C', 'C'] = 1
            print(f'self.deltaDeltaE_pvalues:\n{self.deltaDeltaE_pvalues.to_string()}')
            self.data['deltaDeltaE'] = {'df':self.deltaDeltaE, 'pvalues':self.deltaDeltaE_pvalues}
        else:
            self.data['deltaDeltaE'] = {'df':self.deltaDeltaE, 'pvalues':None}


        # define a data dictionary to loop through when plottin as the plotting method will be the same for all these
        #self.data = {'Ess_deltaE': {'df':self.Ess_deltaE, 'pvalues':self.Ess_deltaE_pvalues}, 
        #        'NonEss_deltaE':{'df':self.NonEss_deltaE, 'pvalues':self.NonEss_deltaE_pvalues}, 
        #        'deltaDeltaE':{'df':self.deltaDeltaE, 'pvalues':self.deltaDeltaE_pvalues},
        #        'Enorm':{'df':self.Enorm, 'pvalues':None}}

        self.tag = self.__dict__['tag']
        print(f'self.tag: {self.tag}')
        

    def plot_data(self,):
        """
        Generates and saves plots for the specified regression variable.
        """
        print('\nGenerates and saves plots for the specified regression variable.')

        ########################################################################################
        # Create the heatmaps for Enorm and deltaE(Ess) and deltaE(NonEss)
        for info, data_dict in self.data.items():
            print(info)
            #print(info, data_dict)

            df = data_dict['df']
            pvalues_df = data_dict['pvalues']

            df = df.replace([np.inf, -np.inf], np.nan).fillna(0)

            plt.figure(figsize=(8, 6))
            heatmap = sns.heatmap(df, annot=df.round().astype(int), fmt="d", cmap="coolwarm", vmin=-60, vmax=60, linewidths=.5, cbar_kws={"shrink": .8, "aspect": 30})
            
            # if there is a valid pvalue dataframe color those with values below 0.05 in a yellow highlight
            if isinstance(pvalues_df, pd.DataFrame):
                print('Pvalues found')
                # Highlight cells based on df2 values
                for i in range(pvalues_df.shape[0]):
                    for j in range(pvalues_df.shape[1]):
                        if pvalues_df.iloc[i, j] < 0.05:
                            heatmap.add_patch(plt.Rectangle((j, i), 1, 1, fill=False, edgecolor='yellow', lw=3))

            # Set labels
            heatmap.set_yticklabels(df.index, rotation=0)
            heatmap.set_xticklabels(df.columns, rotation=90)

            plt.title(info)
            plt.xlabel('Amino Acids')
            plt.ylabel('Amino Acids')
            #plt.show()
            plt.tight_layout()
            outfile = os.path.join(self.outpath, f'{info}_{self.tag}.png')
            plt.savefig(outfile)
            plt.close()
            print(f'SAVED: {outfile}')


def main():
    parser = argparse.ArgumentParser(description="Process regression data and generate plots.")
    parser.add_argument("-En", "--Enorm", type=str, required=True, help="path to Enorm file to plot")
    parser.add_argument("-EsdE", "--Ess_deltaE", type=str, required=True, help="path to Essential gene GT deltaE file to plot")
    parser.add_argument("-EsdEp", "--Ess_deltaE_pvalues", type=str, required=True, help="path to fdr pvalues file")
    parser.add_argument("-NEsdE", "--NonEss_deltaE", type=str, required=True, help="path to NonEssential gene GT deltaE file to plot")
    parser.add_argument("-NEsdEp", "--NonEss_deltaE_pvalues", type=str, required=True, help="path to fdr pvalues file")
    parser.add_argument("-dDE", "--deltaDeltaE", type=str, required=True, help="path to deltaDeltaE file to plot")
    parser.add_argument("-dDEp", "--deltaDeltaE_pvalues", type=str, required=True, help="path to fdr pvalues file")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory.")
    parser.add_argument("-t", "--tag", type=str, required=True, help="tag for outfile: deltaE, Enorm,e ct....")
    args = parser.parse_args()

    Enorm_file = args.Enorm
    Ess_deltaE_file = args.Ess_deltaE
    Ess_deltaE_pvalues_file = args.Ess_deltaE_pvalues
    NonEss_deltaE_file = args.NonEss_deltaE
    NonEss_deltaE_pvalues_file = args.NonEss_deltaE_pvalues
    deltaDeltaE_file = args.deltaDeltaE
    deltaDeltaE_pvalues_file = args.deltaDeltaE_pvalues
    outpath = args.outpath
    tag = args.tag

    if not os.path.exists(outpath):
        os.makedirs(outpath)
        print(f'MADE: {outpath}')

    plotter = Plotter(
            outpath = outpath, 
            Enorm_file = Enorm_file, 
            Ess_deltaE_file = Ess_deltaE_file, 
            Ess_deltaE_pvalues_file = Ess_deltaE_pvalues_file,
            NonEss_deltaE_file = NonEss_deltaE_file,
            NonEss_deltaE_pvalues_file = NonEss_deltaE_pvalues_file,
            deltaDeltaE_file = deltaDeltaE_file, 
            deltaDeltaE_pvalues_file = deltaDeltaE_pvalues_file,
            tag = tag)

    plotter.plot_data()

    print('NORMAL TERMINATION')

if __name__ == "__main__":
    main()
##########################


