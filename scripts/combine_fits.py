import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.rcParams['axes.formatter.useoffset'] = False

import matplotlib.pyplot as plt
import seaborn as sns


if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="remove time info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--rate-table', type=str, required=True, help="input data")
    parser.add_argument('--output-plot', type=str, help="plot file")
    args = parser.parse_args()

    rates = pd.read_csv(args.rate_table, sep='\t', index_col='clade')

    fig, axs = plt.subplots(1,3, figsize=(15,6))
    for ax, mut_type, ax_label in zip(axs,['nuc', 'aa', 'syn'],['total divergence', 'aa divergence', 'syn divergence']):
        ax.set_ylabel(ax_label)
        for clade, row in rates.iterrows():
            dt = np.linspace(0, 0.75,2)
            t = dt + row[f'{mut_type}_origin']
            slope = row[f'{mut_type}_rate']
            ax.plot(t, dt*slope + row[f'{mut_type}_div'], label=clade)
        ax.set_xlim([2019.8, 2022.5])
        ax.set_ylim(0)
        if mut_type == 'syn': ax.legend(ncol=2)


    if args.output_plot:
        plt.savefig(args.output_plot)
    else:
        plt.show()
