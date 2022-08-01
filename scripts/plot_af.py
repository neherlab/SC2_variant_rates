import argparse,json
from datetime import datetime
import numpy as np
import matplotlib as mpl
mpl.rcParams['axes.formatter.useoffset'] = False

import matplotlib.pyplot as plt
import seaborn as sns


def fit_exp(n,k,t, order=0):
    from scipy.optimize import minimize
    eps = 1e-16
    def binom(x,n,k,t, order=0):
        p = x[0]**2*np.exp(-x[1]**2*t)
        if order>0:
            p *= x[1]**2*t
        if order>1:
            p *= x[1]**2*t/2
        p = np.minimum(1.0,p)
        return -np.sum((n-k)*np.log(1-p+eps) + k*np.log(p+eps))

    sol = minimize(binom, (1,0.05), args=(n,k,t, order), method="Powell")

    return {'offset':sol['x'][0]**2, 'rate':sol['x'][1]**2}



if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="remove time info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--counts', nargs='+', type=str, required=True, help="input data")
    parser.add_argument('--output-plot', type=str, help="figure file")
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
    # plt.savefig(args.output_plot)


    fig, axs = plt.subplots(1,1, figsize = (8,6))
    ax = axs
    ls = ['-', '-.', '--']
    for fi,(clade,d) in enumerate(counts.items()):
        dates =   np.array([x for x in d['bins'][:-1]])
        dates -= dates[0]
        total =   np.array(d["all_samples"])
        founder = np.array(d["mutation_number"]['0'])
        ind = founder>0
        if ind.sum()>10:
            print(clade, "founder", fit_exp(total[ind], founder[ind], dates[ind]))
            ax.plot(dates[ind], founder[ind]/total[ind], label=clade, ls=ls[fi//10], c='C0')

        single = np.array(d["mutation_number"]['1'])
        ind = single>0
        if ind.sum()>10:
            print(clade, "single", fit_exp(total[ind], single[ind], dates[ind], order=1))
            ax.plot(dates[ind], single[ind]/total[ind], label=clade, ls=ls[fi//10], c='C1')

        double = np.array(d["mutation_number"]['2'])
        ind = double>0
        if ind.sum()>10:
            print(clade, "double", fit_exp(total[ind], double[ind], dates[ind], order=2))
            ax.plot(dates[ind], double[ind]/total[ind], label=clade, ls=ls[fi//10], c='C2')

    ax.set_yscale('log')
    ax.set_xlabel('mutation frequency')
    ax.set_ylabel('fraction above')
    # x = np.logspace(-4,-0.3,21)
    # ax.plot([1e-5,1e-1], [1e-0,1e-4], c='k', lw=3, alpha=0.5, label='1/x')
    # ax.plot(x, np.log(x)/np.log(x[0]), c='k', lw=3, alpha=0.5, label='~log(x)')
    plt.legend(ncol=3)

