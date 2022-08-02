import argparse,json
from datetime import datetime
import numpy as np
import matplotlib as mpl
mpl.rcParams['axes.formatter.useoffset'] = False

import matplotlib.pyplot as plt
import pandas as pd

from plot_genotype_counts import fit_poisson


if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="remove time info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--counts', nargs='+', type=str, required=True, help="input data")
    parser.add_argument('--output-plot', type=str, help="figure file")
    parser.add_argument('--output-rates', type=str, help="table")
    args = parser.parse_args()
    counts = {}
    for fi,fname in enumerate(args.counts):
        clade = fname.split('/')[-1].split('_')[0]
        with open(fname) as fh:
            counts[clade] = json.load(fh)

    fig, axs = plt.subplots(1,1, figsize = (8,6))
    ax = axs
    ls = ['-', '-.', '--']
    for fi,(clade,d) in enumerate(counts.items()):
        ax.plot(sorted(d["mutation_spectrum"].values()), np.linspace(1,0,len(d["mutation_spectrum"])+1)[:-1],
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


    ls = ['-', '-.', '--']
    data = []
    for fi,(clade,d) in enumerate(counts.items()):
        dates =   np.array([x for x in d['bins'][:-1]])
        datetimes =   np.array([datetime.fromordinal(x) for x in d['bins'][:-1]])
        t0 = dates[0]
        dates -= t0
        total =   np.array(d["all_samples"])
        k_vecs = { i: np.array(d["mutation_number"][f'{i}']) for i in ['0', '1', '2', '3', '4', '5+']}

        ind = total>0
        res = fit_poisson(total[ind], {x:k_vecs[x][ind] for x in k_vecs}, dates[ind])
        start_date = datetime.fromordinal(int(t0 + res['offset'])).strftime("%Y-%m-%d")

        res_fixed_rate = fit_poisson(total[ind], {x:k_vecs[x][ind] for x in k_vecs}, dates[ind], rate=15/365)
        start_date_fixed = datetime.fromordinal(int(t0 + res_fixed_rate['offset'])).strftime("%Y-%m-%d")

        data.append({'clade':clade, 'rate':res['rate']*365, 'origin':start_date, 'origin_fixed_rate':start_date_fixed})

    pd.DataFrame(data).to_csv(args.output_rates, sep='\t')
