import argparse,json
from datetime import datetime
import numpy as np
import matplotlib as mpl
mpl.rcParams['axes.formatter.useoffset'] = False

import matplotlib.pyplot as plt
import seaborn as sns


if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="remove time info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--counts', nargs='+', type=str, required=True, help="input data")
    parser.add_argument('--output-plot', type=str, help="figure file")
    args = parser.parse_args()

    fig, axs = plt.subplots(1,1, figsize = (8,6))
    ax = axs
    ls = ['-', '-.', '--']
    for fi,fname in enumerate(args.counts):
        with open(fname) as fh:
            counts = json.load(fh)

        clade = fname.split('/')[-1].split('_')[0]
        dates = [datetime.fromordinal(x) for x in counts['bins'][:-1]]

        ax.plot(sorted(counts["mutation_spectrum"].values()), np.linspace(1,0,len(counts["mutation_spectrum"])+1)[:-1],
                       label=clade, ls=ls[fi//10])

    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('mutation frequency')
    ax.set_ylabel('fraction above')
    x = np.logspace(-4,-0.3,21)
    ax.plot([1e-5,1e-1], [1e-0,1e-4], c='k', lw=3, alpha=0.5, label='1/x')
    # ax.plot(x, np.log(x)/np.log(x[0]), c='k', lw=3, alpha=0.5, label='~log(x)')
    plt.legend(ncol=3)
    plt.savefig(args.output_plot)
