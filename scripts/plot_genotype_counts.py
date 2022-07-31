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

    parser.add_argument('--counts', type=str, required=True, help="input data")
    parser.add_argument('--output-plot', type=str, help="figure file")
    args = parser.parse_args()

    with open(args.counts) as fh:
        counts = json.load(fh)

    plt.figure()
    nmuts = [m for m in counts['mutation_number'] if m[-1]!='+']
    for week in range(5,len(counts['mutation_number']['0'])):
        c = np.array([counts['mutation_number'][str(n)][week] for n in nmuts])
        if c.sum()>100:
            plt.plot(range(len(nmuts)), c/c.sum(), '-')


    dates = np.array([datetime.fromordinal(x) for x in counts['bins'][:-1]])
    fig, axs = plt.subplots(1,3, figsize = (18,6))

    ax = axs[2]
    # ax.plot(dates, counts['all_samples'], lw=3, c='k', alpha=0.3)
    all_samples = np.array(counts['all_samples'])
    for m in ['0', '1', '2', '3', '4', '5+']:
        muts = np.array(counts['mutation_number'][m])
        ind = muts>0
        ax.plot(dates[ind], muts[ind]/all_samples[ind], '-o', label=f'{m} mutations')

    ax.set_yscale('log')
    ax.legend()

    ax = axs[1]
    ax.plot(dates, counts['all_samples'], lw=3, c='k', alpha=0.3)
    for m in sorted(counts['mutations'].keys(), key=lambda x:int(x[1:-1])):
        ax.plot(dates, counts['mutations'][m], '-o', label=f'{m}')

    ax.set_yscale('log')
    ax.legend()

    ax = axs[0]
    ax.plot(dates, counts['all_samples'], lw=3, c='k', alpha=0.3)
    for m in sorted(counts['genotypes'].keys(), key=lambda x: len(x)):
        ax.plot(dates, counts['genotypes'][m], '-o', label=f'{m}' if m else "founder")

    ax.set_yscale('log')
    ax.legend()
    fig.autofmt_xdate()

    # ax = axs[0,1]
    # ax.plot(sorted(counts["mutation_spectrum"].values()), np.linspace(1,0,len(counts["mutation_spectrum"])))
    # ax.set_yscale('log')
    # ax.set_xscale('log')
    plt.savefig(args.output_plot)
