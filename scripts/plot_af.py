import argparse,json
from datetime import datetime
import numpy as np
import matplotlib as mpl
mpl.rcParams['axes.formatter.useoffset'] = False

import matplotlib.pyplot as plt
import pandas as pd


def fit_exp(n,k_list,t):
    from scipy.optimize import minimize
    eps = 1e-16
    def binom(x,n,k_list,t):
        res = 0
        tu = np.maximum(0,x[1]**2*(t-x[0]))
        for order, k in enumerate(k_list):
            p = np.exp(-tu)
            if order>0:
                p *= tu
            if order>1:
                p *= tu/2
            if order>2:
                p *= tu/3
            if order>3:
                p *= tu/4
            p = np.minimum(1.0,p)
            res -= np.sum((n-k)*np.log(1-p+eps) + k*np.log(p+eps))
        return res

    sol = minimize(binom, (-10,0.25), args=(n,k_list,t), method="Powell")

    return {'offset':sol['x'][0], 'rate':sol['x'][1]**2}



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
        fig = plt.figure()
        dates =   np.array([x for x in d['bins'][:-1]])
        datetimes =   np.array([datetime.fromordinal(x) for x in d['bins'][:-1]])
        t0 = dates[0]
        dates -= t0
        total =   np.array(d["all_samples"])
        klist = [ np.array(d["mutation_number"]['0']),
                  np.array(d["mutation_number"]['1']),
                  np.array(d["mutation_number"]['2']),
                  np.array(d["mutation_number"]['3']),
                  np.array(d["mutation_number"]['4'])]
        ind = total>0
        res = fit_exp(total[ind], [x[ind] for x in klist], dates[ind])
        print(clade, res)
        plt.plot(datetimes[ind], klist[0][ind]/total[ind], 'o',c='C0')
        plt.plot(datetimes[ind], klist[1][ind]/total[ind], 'o', c='C1')
        plt.plot(datetimes[ind], klist[2][ind]/total[ind], 'o', c='C2')
        plt.plot(datetimes[ind], klist[3][ind]/total[ind], 'o', c='C3')
        plt.plot(datetimes[ind], klist[4][ind]/total[ind], 'o', c='C4')

        tu = (dates[ind]-res['offset'])*res['rate']
        plt.plot(datetimes[ind], np.exp(-tu), label=clade,         ls='-', c='C0')
        plt.plot(datetimes[ind], np.exp(-tu)*tu ,label=clade,      ls='-', c='C1')
        plt.plot(datetimes[ind], np.exp(-tu)*tu**2/2, label=clade, ls='-', c='C2')
        plt.plot(datetimes[ind], np.exp(-tu)*tu**3/6 ,label=clade, ls='-', c='C3')
        plt.plot(datetimes[ind], np.exp(-tu)*tu**4/24, label=clade, ls='-', c='C4')

        plt.yscale('log')
        start_date = datetime.fromordinal(int(t0 + res['offset'])).strftime("%Y-%m-%d")
        data.append({'clade':clade, 'rate':res['rate']*365, 'origin':start_date})
        plt.title(f"{clade}: {start_date}, {res['rate']*365:1.3}/year")
        fig.autofmt_xdate()
        plt.savefig(f'figures/{clade}_poisson.pdf')

    pd.DataFrame(data).to_csv(args.output_rates, sep='\t')
