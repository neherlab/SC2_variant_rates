import argparse,json
from datetime import datetime
from sys import exec_prefix
import numpy as np
import matplotlib as mpl
mpl.rcParams['axes.formatter.useoffset'] = False

import matplotlib.pyplot as plt
import seaborn as sns
from root_to_tip import add_panel_label
fs=12
def poisson(lam, order):
    p = np.exp(-lam)
    try:
        for i in range(0,int(order)):
            p *= lam/(i+1)
    except:
        p_cum = np.copy(p)
        for i in range(0,int(order[:-1])-1):
            p *= lam/(i+1)
            p_cum += p
        p = 1 - p_cum

    return np.maximum(0.0, np.minimum(1.0,p))

def fit_poisson(n,k_vecs,t, rate=None):
    from scipy.optimize import minimize
    eps = 1e-16
    def binom(x,n,k_vecs,t, rate):
        res = 0
        mu = rate or x[1]**2
        tu = np.maximum(0,mu*(t-x[0]))
        for order, k in k_vecs.items():
            p = poisson(tu, order)
            res -= np.sum(((n-k)*np.log(1-p+eps) + k*np.log(p+eps)))

        return res

    if rate is None:
        sol = minimize(binom, (-10,0.25), args=(n,k_vecs,t,rate), method="Powell")
        res = {'offset':sol['x'][0], 'rate':sol['x'][1]**2}
    else:
        sol = minimize(binom, (-10,), args=(n,k_vecs,t,rate), method="Powell")
        res = {'offset':sol['x'][0], 'rate':rate}

    return res


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
    t0 = counts['bins'][0]
    rel_date = np.array([x-t0 for x in counts['bins'][:-1]])
    fig, axs = plt.subplots(1,3, figsize = (18,6))

    ax = axs[2]
    # ax.plot(dates, counts['all_samples'], lw=3, c='k', alpha=0.3)
    total = np.array(counts['all_samples'])
    ind = total>0
    k_vecs = {m: np.array(counts['mutation_number'][m]) for m in ['0', '1', '2', '3', '4', '5+']}
    res = fit_poisson(total[ind], {x:k_vecs[x][ind] for x in k_vecs}, rel_date[ind])
    tu = (rel_date[ind]-res['offset'])*res['rate']
    for i, (m,k) in enumerate(k_vecs.items()):
        ax.plot(dates[ind], k[ind]/total[ind], 'o', label=f'{m} mutations',c=f"C{i}")
        plt.plot(dates[ind], poisson(tu,m), ls='-', c=f'C{i}')

    ax.set_ylim(8e-5, 2)

    ax = axs[1]
    ax.plot(dates, total, lw=3, c='k', alpha=0.3)
    for m in sorted(counts['mutations'].keys(), key=lambda x:int(x[1:-1])):
        ax.plot(dates, counts['mutations'][m], '-o', label=f'{m}')

    ax = axs[0]
    ax.plot(dates, total, lw=3, c='k', alpha=0.3)
    for m in sorted(counts['genotypes'].keys(), key=lambda x: len(x)):
        ax.plot(dates, counts['genotypes'][m], '-o', label=f'{m}' if m else "founder")

    fig.autofmt_xdate()
    for ax, label in zip(axs, 'DEF'):
        ax.set_yscale('log')
        ax.legend()
        add_panel_label(ax, label, fs=fs*1.8)


    plt.savefig(args.output_plot)
