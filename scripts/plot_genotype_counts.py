import argparse,json
from datetime import datetime
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


    dates = [datetime.fromordinal(x) for x in counts['bins'][:-1]]
    fig, axs = plt.subplots(3,1, figsize = (8,20), sharex=True)

    axs[0].plot(dates, counts['all_samples'], lw=3, c='k', alpha=0.3)
    for m in counts['mutation_number']:
        axs[0].plot(dates, counts['mutation_number'][m], '-o', label=f'{m} mutations')

    axs[0].set_yscale('log')
    axs[0].legend()

    axs[1].plot(dates, counts['all_samples'], lw=3, c='k', alpha=0.3)
    for m in sorted(counts['mutations'].keys(), key=lambda w:int(m[1:-1])):
        axs[1].plot(dates, counts['mutations'][m], '-o', label=f'{m}')

    axs[1].set_yscale('log')
    axs[1].legend()

    axs[2].plot(dates, counts['all_samples'], lw=3, c='k', alpha=0.3)
    for m in sorted(counts['genotypes'].keys(), key=lambda x: len(x)):
        axs[2].plot(dates, counts['genotypes'][m], '-o', label=f'{m}' if m else "founder")

    axs[2].set_yscale('log')
    axs[2].legend()
    fig.autofmt_xdate()

    plt.savefig(args.output_plot)
